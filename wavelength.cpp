#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vectors.h"
#include "image.h"
#include <vector>
#include <algorithm>
#include "photon.h"

Vec3f wavelengthToRGB(double wavelength){
	double r = 0, g = 0, b = 0;
	if(wavelength < 410 && wavelength >= 380){
		r = .6 - .41 * (wavelength - 380) / 30;
		b = .39 + .6 * (wavelength - 380) / 30;
	}
	else if(wavelength < 440 && wavelength >= 410){
		r = .19 - .19 * (wavelength - 410) / 30;
		b = 1;
	}
	else if(wavelength < 490 && wavelength >= 440){
		g = 1 - (490 - wavelength) / 50;
		b = 1;
	}
	else if(wavelength < 510 && wavelength >= 490){
		g = 1;
		b = (510 - wavelength) / 20;
	}
	else if(wavelength < 580 && wavelength >= 510){
		r = 1 - (580 - wavelength) / 70;
		g = 1;
	}
	else if(wavelength < 640 && wavelength >= 580){
		r = 1;
		g = (640 - wavelength) / 60;
	}
	else if(wavelength < 700 && wavelength >= 640){
		r = 1;
	}
	else if(wavelength <= 780 && wavelength >= 700){
		r = 1 - .65 * (wavelength - 700) / 80;

	}
	return Vec3f(r,g,b);
}
Vec3f RGBtoXYZ(Vec3f rgb){
	Vec3f xyz;
	Matrix transform;
	transform.clear();
	transform.setToIdentity();
	
	transform.set(0,0,0.67);
	transform.set(0,1,0.21);
	transform.set(0,2,0.14);

	transform.set(1,0,0.33);
	transform.set(1,1,0.71);
	transform.set(1,2,0.08);

	transform.set(2,0,0.00);
	transform.set(2,1,0.08);
	transform.set(2,2,0.78);
	return transform * rgb;
}
Vec3f XYZtoRGB(Vec3f xyz){
	Matrix transform;
	transform.clear();
	transform.setToIdentity();
	
	transform.set(0,0,0.67);
	transform.set(0,1,0.21);
	transform.set(0,2,0.14);

	transform.set(1,0,0.33);
	transform.set(1,1,0.71);
	transform.set(1,2,0.08);

	transform.set(2,0,0.00);
	transform.set(2,1,0.08);
	transform.set(2,2,0.78);
	transform.Inverse();

	return transform * xyz;
}
//Given a list of wavelengths, converts to xyY coordinates to add
//chromaticity and luminance, then converts an RGB coordinate system
//for rendering.
//RGB and XYZ coordinates for a specific wavelength is generated by
//determining the color that most closely represents the dominant
//wavelength in RGB representation.
const Vec3f& mixColors(std::vector<Photon> wavelengths){
	if(wavelengths.size() == 0)
		return Vec3f(0,0,0);
	double r = 0, g = 0, b = 0;
	for(int i = 0; i < wavelengths.size(); i++){
		Vec3f rgb = wavelengthToRGB(wavelengths[i].getWavelength());
		r += rgb.r();
		g += rgb.g();
		b += rgb.b();
	}
	r = (r / wavelengths.size()) / 255;
	g = (g / wavelengths.size()) / 255;
	b = (b / wavelengths.size()) / 255;
	double mult = 1;// std::min( std::min(1.0 / r, 1.0 / g), 1.0 / b );
	return Vec3f(r * mult, g * mult, b * mult);
} 
/*
int main() {
	Image full_spectrum, red, green, blue;
	full_spectrum.Allocate(400,10);
	red.Allocate(400,10);
	green.Allocate(400,10);
	blue.Allocate(400,10);
	for(double lambda = 410; lambda < 700; lambda++){
		Vec3f rgb = wavelengthToRGB(lambda);
		int i = (int)(lambda - 380);
		for(int j = 0; j < 10; j++){
			full_spectrum.SetPixel(i, j, Color((int)rgb.r(), (int)rgb.g(), (int)rgb.b()));
			red.SetPixel(i, j, Color((int)rgb.r(), 0, 0));
			green.SetPixel(i, j, Color(0, (int)rgb.g(), 0));
			blue.SetPixel(i, j, Color(0, 0, (int)rgb.b()));
		}
	}
	full_spectrum.Save("full.ppm");
	red.Save("red.ppm");
	green.Save("green.ppm");
	blue.Save("blue.ppm");
	Image cie2;
	cie2.Allocate(400,400);
	cie2.SetAllPixels(Color(0,0,0));
	for(int i = 0; i < 255; i++){
		for(int j = 0; j < 255; j++){
			for(int k = 0; k < 255; k++){
				Vec3f XYZ = RGBtoXYZ(Vec3f(i,j,k));
				double x = XYZ.x() / (XYZ.x() + XYZ.y() + XYZ.z());
				double y = XYZ.y() / (XYZ.x() + XYZ.y() + XYZ.z());
				if(x >= 0 && y >= 0 && x <=1 && y <= 1)
					cie2.SetPixel((int)(x * 400),(int)(y * 400), Color(i,j,k));
			}
		}
	}
	Image cie;
	cie.Allocate(400,400);
	cie.SetAllPixels(Color(0,0,0));
	std::vector<double> wavelengths;
	for(double i = 380.0; i < 780.0; i += 0.01){
		wavelengths.push_back(i);
		Vec3f XYZ = RGBtoXYZ(wavelengthToRGB(i));
		double X = XYZ.x();
		double Y = XYZ.y();
		double Z = XYZ.z();
		Vec3f RGB = XYZtoRGB(XYZ);
		double x = X / (X + Y + Z);
		double y = Y / (X + Y + Z);
		int k = x * 400;
		int j = y * 400;
		cie.SetPixel(k, j, Color(wavelengthToRGB(i).r(), wavelengthToRGB(i).g(), wavelengthToRGB(i).b()));
	}
	std::cout << "Expecting White Light: " << mixColors(wavelengths) << std::endl;
	cie.Save("CIE.ppm");
	cie2.Save("CIE2.ppm");
}
*/

