#include <QImage>
#include <QColor>
#include <math.h>
#include <vector>
#include "algebra.h"

using namespace std;

vector<Shape*> shapes;
vector<Sphere*> lights;
int maxDepth = 4;
Vector3D ambient = Vector3D(0.1);

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
            if (shapes[j]->intersect(p + 0.01*normal, rayToLightSource, t0, t1))
            {
                return ret;
            }
        }
        // diffuse component
        Vector3D Ld = lights[i]->ld;
        Vector3D kd = C->kd;
        Vector3D Id = max(normal.dot(rayToLightSource), 0.0) * kd * Ld;

        float shininess = 30.0;
        // specular componenet
        Vector3D reflectionVector = 2 * (max(normal.dot(rayToLightSource),0.0)) * normal - rayToLightSource;
        reflectionVector.normalize();
        Vector3D viewerVector = pToViewer;
        Vector3D Ls = lights[i]->ls;
        Vector3D ks = C->ks;
        Vector3D Is = pow(max(reflectionVector.dot(viewerVector), 0.0), shininess) * ks * Ls;

        ret = ret + (Id + Is) * C->surfaceColour;
    }
    return ret;
}

Vector3D trace(Point3D origin, Vector3D direction, int depth)
{
    Vector3D local, reflected;
    Point3D q; //intersection point
    Vector3D n, r; // normal, reflection

    if (depth >= maxDepth) return Vector3D();

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
    Vector3D pToViewer = -direction;
    local = phong(q, n, pToViewer, closest);
    reflected = trace(q + 0.01*n, n, depth+1);

    Vector3D ret =  local + reflected;
    ret[0] = min(ret[0],1.0);
    ret[0] = max(ret[0],0.0);
    ret[1] = min(ret[1],1.0);
    ret[1] = max(ret[1],0.0);
    ret[2] = min(ret[2],1.0);
    ret[2] = max(ret[2],0.0);
    return ret;
}

int main(int argc, char *argv[])
{
	// currently unused parameters
	Q_UNUSED(argc);
	Q_UNUSED(argv);

    // image width and height
    if (argc < 3)
    {
        cout << "usage: " << argv[0] << " [width] [height]" << endl;
        return 1;
    }
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);

	// create new image
	QImage image(width, height, QImage::Format_RGB32);

    // Spherere center, radius, surfaceColour, kd, ks;

    // spheres from that sample
    shapes.push_back(new Sphere(Point3D(0.0,0.0,-20), 4, Vector3D(1,0.32,0.36), Vector3D(0.5), Vector3D(0.5)));
    shapes.push_back(new Sphere(Point3D(5,-1,-15), 2, Vector3D(0.9,0.76,0.46), Vector3D(0.2), Vector3D(0.2)));
    shapes.push_back(new Sphere(Point3D(5,0,-25), 3, Vector3D(0.65,0.77,0.97), Vector3D(0.5), Vector3D(0.5)));
    shapes.push_back(new Sphere(Point3D(-5.5,-0,-15), 3, Vector3D(0.9,0.90,0.90), Vector3D(0.75,0.60,0.22), Vector3D(0.5)));
    shapes.push_back(new Sphere(Point3D(0, 25,-40), 15, Vector3D(0.3), Vector3D(0.25), Vector3D(0.25)));

    //plane(plane(Point, normal, colour, kd, ks, max coord, min coord)
    shapes.push_back(new Quad(Point3D(0,0,-53), Vector3D(0.3,0.9,1), Vector3D(0.5), Vector3D(0.2), Vector3D(0.2), Vector3D(30,30,-60), Vector3D(-30,-40,-30)));

    //triangle(Point a, Point b, Point c, colour, kd, ks)
/*    shapes.push_back(new MyTriangle(Point3D(0,-5,0), Point3D(1,-4,-1), Point3D(-1,-6,-1),
                                    Vector3D(0.5), Vector3D(0.5), Vector3D(0.5)))*/;

    // Lights are spheres with 2 optional additional argument for ld and ls
    lights.push_back(new Sphere(Point3D(0,20,-20), 0.01, Vector3D(0.3), Vector3D(1), Vector3D(1), Vector3D(0.7), Vector3D(0.7)));
    lights.push_back(new Sphere(Point3D(-50,50,-50), 0.01, Vector3D(0.3), Vector3D(1), Vector3D(1), Vector3D(1), Vector3D(1)));

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

			// set pixel value
			image.setPixel(x, y, 
                qRgb(r*255, g*255, b*255));
		}
	}

    // save to file
	image.save("output.png");
	
	// application successfully returned
	return 0;
}
