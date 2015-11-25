#include <QImage>
#include <QColor>

int main(int argc, char *argv[])
{
	// currently unused parameters
	Q_UNUSED(argc);
	Q_UNUSED(argv);

	// image width and height
	// TODO: prompt user on command line for dimensions
	int width = 200;
	int height = 200;

	// create new image
	QImage image(width, height, QImage::Format_RGB32);

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
