#ifndef _PHOTON_H_
#define _PHOTON_H_

#include <iostream>
#include "vectors.h"

// ===========================================================
// Class to store the information when a photon hits a surface

class Photon {
 public:
	Photon() {}
	
	// CONSTRUCTOR
	Photon(const Vec3f &p, const Vec3f &d, const float &w, int b) :
		position(p),direction_from(d),bounce(b),wavelength(w){}

	// ACCESSORS
	const Vec3f& getPosition() const { return position; }
	const Vec3f& getDirectionFrom() const { return direction_from; }
	//const Vec3f& getEnergy() const { return energy; }
	int whichBounce() const { return bounce; }

	float getWavelength() const { return wavelength; }

	double getDistance() const { return position.Length(); }

	//Vec3f getNormal() const { return normal; }

    void setPosition(const Vec3f& p) { position = Vec3f(p); }   

/* //I don't think that this function is needed because photons are not actually passed into TracePhoton
	const float getCurrentNVal() const {

		//if there is no current material, use refractive index of air
		//TODO: update current refractive index to be air again when leaving a material (hit will still be of that material's type)
		if (current_material == NULL){
			return 1.000293;
		}
		else {
			return current_material->getRefractiveIndex();
		}
	}
*/

 private:
	// REPRESENTATION
	Vec3f position;
	Vec3f direction_from;
	//Vec3f normal;
	//Vec3f energy;
	int bounce;


	float wavelength;

	//Material* current_material;

};

inline bool compareLengths (const Photon & i, const Photon & j) { return i.getDistance() < j.getDistance(); }

#endif
