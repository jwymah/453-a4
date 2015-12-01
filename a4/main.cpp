#include <QImage>
#include <QColor>
#include <math.h>
#include <vector>
#include "algebra.h"

using namespace std;

vector<Shape*> shapes;
vector<Sphere*> lights;
Vector3D bgColour = Vector3D();
int maxDepth = 2;
float ambient = 0.35;

Vector3D reflect(Vector3D incident, Vector3D normal)
{
    incident.normalize();
    normal.normalize();

    Vector3D reflected = incident - 2*((float) incident.dot(normal)*normal);
    reflected.normalize();
    return reflected;
}

Vector3D phong(Point3D p, Vector3D normal, Vector3D pToViewer, Shape *C)
{
    // ambient component
    Vector3D ret = ambient * C->surfaceColour;

    for (int i=0; i<lights.size(); i++)
    {
        Vector3D rayToLightSource = lights[i]->center - p;
        rayToLightSource.normalize();
        for (int j=0; j<shapes.size(); j++)
        {
            float t0, t1;
            if (shapes[j]->intersect(p + 0.0005*normal, rayToLightSource, t0, t1))
            {
                return ret;
            }
        }
        // diffuse component
        Vector3D Ld = Vector3D(1,1,1);
        Vector3D kd = Vector3D(0.75,0.60,0.22);
        Vector3D Id = max(normal.dot(rayToLightSource), 0.0) * kd * Ld;

        Vector3D reflectionVector = 2 * (normal.dot(rayToLightSource)) * normal - rayToLightSource;
        Vector3D viewerVector = pToViewer;
        Vector3D Ls = Vector3D(1,1,1);
        float shininess = 0.4; //TODO shininess should be in shapeclass
        // specular componenet
        Vector3D Is = C->ks * pow(max(reflectionVector.dot(viewerVector), 0.0), shininess) * Ls;

        ret = ret + (Id + Is) * C->surfaceColour;
//        ret = (Is) * C->surfaceColour;
//        ret = (Id) * C->surfaceColour;
    }
    ret[0] = min(ret[0],255.0);
    ret[0] = max(ret[0],0.0);
    ret[1] = min(ret[1],255.0);
    ret[1] = max(ret[1],0.0);
    ret[2] = min(ret[2],255.0);
    ret[2] = max(ret[2],0.0);
    return ret;
}

Vector3D trace(Point3D origin, Vector3D direction, int depth)
{
    Vector3D local, reflected, transmitted;
    Point3D q; //intersection point
    Vector3D n, r , t; // normal, reflection, transmission;

    if (depth > maxDepth) return bgColour;

    Shape *closest = NULL;
    float tnear = INFINITY;

    for (int i =0; i< shapes.size(); i++)
    {
        float t0 = INFINITY, t1 = INFINITY;
        if (shapes[i]->intersect(origin, direction, t0, t1))
        {
            if (t0 < 0)
            {
                t0 = t1;
            }
            if (t0 < tnear)
            {
                tnear = t0;
                closest = shapes[i];
            }
        }
    }
    if (!closest) return Vector3D();    // no intersection

    q = origin + tnear*direction;
    n = closest->normalAt(q);   // assuming it DOES intersect
    r = reflect(direction, n);
//    t = transmit(q, n); //refraction.. is this vector going into shape or leaving it? I think going into
    Vector3D pToViewer = -direction;
    local = phong(q, n, pToViewer, closest);
    reflected = trace(q, n, depth+1);
//    transmitted = trace(q, t, depth+1);
    return local;
//    return closest->surfaceColour;
    //    return local + reflected + transmitted;
}

int main(int argc, char *argv[])
{
	// currently unused parameters
	Q_UNUSED(argc);
	Q_UNUSED(argv);

	// image width and height
	// TODO: prompt user on command line for dimensions
    int width = 640;
    int height = 480;

	// create new image
	QImage image(width, height, QImage::Format_RGB32);

    // Spherere center, radius, colour, emissionColour; ks;
//    shapes.push_back(new Sphere(Vector3D(0,0,-5),   3,      Vector3D(1,1,1),    Vector3D(1,1,1))); //big one at the bottom
//    shapes.push_back(new Sphere(Vector3D(0,0,-6),   3,      Vector3D(1,0,1),    Vector3D(1,1,1))); //big one at the bottom
//    shapes.push_back(new Sphere(Vector3D(0,0,-7),   3,      Vector3D(1,0,0),    Vector3D(1,1,1))); //big one at the bottom
//    shapes.push_back(new Sphere(Vector3D(0,5,-5),   0.5,    Vector3D(0.8,0.75,1),   Vector3D(1,1,1))); // medium one at the top

//    shapes.push_back(new Sphere(Vector3D(-1,3,-5),   0.5,    Vector3D(0.8,0.75,0.1),   Vector3D(1,1,1)));
//    shapes.push_back(new Sphere(Vector3D(-2,4,-5),   0.5,    Vector3D(0.8,0.95,0),   Vector3D(1,1,1)));

//    shapes.push_back(new Sphere(Vector3D(0,4,-5),   0.25,    Vector3D(0.5,0,0),   Vector3D(1,1,1)));

    // spheres from that sample
    shapes.push_back(new Sphere(Point3D(0,-10004,-20), 10000, Vector3D(0.2,0.2,0.2), 0.5));
    shapes.push_back(new Sphere(Point3D(0.0,0.0,-20), 4, Vector3D(1,0.32,0.36), 0.5));
    shapes.push_back(new Sphere(Point3D(5,-1,-15), 2, Vector3D(0.9,0.76,0.46), 0.8));
    shapes.push_back(new Sphere(Point3D(5,0,-25), 3, Vector3D(0.65,0.77,0.97), 0.5));
    shapes.push_back(new Sphere(Point3D(-5.5,-0,-15), 3, Vector3D(0.9,0.90,0.90),0.3));

    lights.push_back(new Sphere(Point3D(0,20,-5), 3, Vector3D(), 1));

	// iterate over the pixels & set colour values
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
        {
            float fov = 45; // 45 degrees in either direction
            float angle = tan(M_PI * fov / 180);
            float invWidth = 1/(float) width;
            float invHeight = 1/(float) height;
            float aspectratio = width / float(height);

            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio; //center image, adjust to fov and ratio
            float yy = (1 - 2 * ((y + 0.5) * invHeight) * angle);

            Vector3D ray = Vector3D(xx, yy, -1);
            ray.normalize();
            Vector3D colour = trace(Point3D(), ray, 0);

            double r = colour[0];
            double g = colour[1];
            double b = colour[2];

//            // try to generate circle at (0,0) diameter half height)
//            if (pow(x-width/2,2) + pow(y-height/2,2) <= pow(height/10,2))
//            {
//                r = 0.5;
//                g = 0.5;
//                b = 0.5;
//            }

			// set pixel value
			image.setPixel(x, y, 
				qRgb(r * 255, g * 255, b * 255));
		}
	}



	// save to file
	// TODO: prompt user on command line for output name
	image.save("output.png");
	
	// application successfully returned
	return 0;
}
