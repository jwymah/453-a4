#include <QImage>
#include <QColor>
#include <math.h>
#include <vector>
#include "algebra.h"

using namespace std;

vector<Shape*> shapes;
//bgColour = Vector3D()
//Vector3D normal(Point3D q)

Vector3D trace(Point3D origin, Vector3D direction, int depth)
{
    Vector3D local, reflected, transmitted;
    Point3D q; //intersection point
    Vector3D n, r , t; // normal, reflection, transmission;

//    if (depth > max) return bgColour;

//    n = normal(q);
//    r = reflect(q, n);
//    t = transmit(q, n);
//    local = phong(1, n, r);;  // and shadows

//    reflected = trace(q, r, depth+1);
//    transmitted = trace(q, t, depth+1);
//    return local + reflected + transmitted;
    return Vector3D();
}

int main(int argc, char *argv[])
{
	// currently unused parameters
	Q_UNUSED(argc);
	Q_UNUSED(argv);

	// image width and height
	// TODO: prompt user on command line for dimensions
    int width = 400;
    int height = 400;

	// create new image
	QImage image(width, height, QImage::Format_RGB32);

    // Spherere center, radius, colour, Vector3D surfaceColour, emissionColour; float transparency, reflectivity;
    shapes.push_back((new Sphere(Vector3D(0,0,-5), 3, Vector3D(1,1,1), Vector3D(1,1,1)))); //big one at the bottom
    shapes.push_back((new Sphere(Vector3D(0,5,-5), 0.5, Vector3D(0.8,0.75,1), Vector3D(1,1,1)))); // medium one at the top
//    shapes.push_back(new Sphere(Vector3D(0,0,-20), 4, Vector3D(1,0.32,0.36), Vector3D(1,1,1)));


	// iterate over the pixels & set colour values
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			double r = 0.0;
			double g = 0.0;
			double b = 0.0;

			// compute rgb values
			// TODO: replace with values from ray tracing
			if (x > width / 2 && y < height / 2) r = 1.0;
			if (x < width / 2 && y > height / 2) g = 1.0;
			if (x > width / 2 && y > height / 2) b = 1.0;

            float aspectratio = width / float(height);
            float fov = 30;
            float angle = tan(M_PI * 0.5 * fov / 180);

            float xx = (2 * ((x + 0.5) * 1/float(width)) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * (1/float(height))) * angle);
            Vector3D ray = Vector3D(xx, yy, -1);
            ray.normalize();

            for (int i =0; i< shapes.size(); i++)
            {
                if (shapes[i]->intersect(Vector3D(), ray))
                {
                    r = 1;
                    g = 1;
                    b = 1;
                }
            }
            // try to generate circle at (0,0) diameter half height)
            if (pow(x-width/2,2) + pow(y-height/2,2) <= pow(height/10,2))
            {
                r = 0.5;
                g = 0.5;
                b = 0.5;
            }

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
