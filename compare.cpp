#include <stdlib.h>
#include <cassert>
#include <math.h>
#include "vectors.h"
#include "image.h"

Image calculateImageDifference(Image& a, Image& c){
	assert(a.Width() == c.Width() && a.Height() == c.Height());
	Image difference;
	difference.Allocate(a.Width(), a.Height());
	for(int i = 0; i < a.Width(); i++){
		for(int j = 0; j < a.Height(); j++){
			int r = abs(a.GetPixel(i,j).r - c.GetPixel(i,j).r);
			int g = abs(a.GetPixel(i,j).g - c.GetPixel(i,j).g);
			int b = abs(a.GetPixel(i,j).b - c.GetPixel(i,j).b);
			difference.SetPixel(i, j, Color(r, g, b));
		}
	}
	return difference;
}

int main() {
	Image full_spectrum;
	Image full_bad;
	full_spectrum.Load("debug/full.ppm");
	full_bad.Load("full.ppm");
	Image diff = calculateImageDifference(full_spectrum, full_bad);
	diff.Save("diff.ppm");
}
