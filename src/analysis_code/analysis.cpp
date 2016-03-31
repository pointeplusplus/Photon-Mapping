
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#define PI 3.1415926535897

class Vector{
public:
	Vector(){
		x = 0;
		y = 0;
		z = 0;
	}
	Vector(double a, double b, double c){
		x = a;
		y = b;
		z = c;
	}

	double x;
	double y; 
	double z;

};

class Face{
public:
	Face(){
		area = 0;
		internal_bounces = 0;
		reflected_rays = 0;
		rays_entering = 0;
		rays_leaving = 0;
	}
	Face(std::string m, Vector n, double a, unsigned int i, unsigned int r, unsigned int e, unsigned int l){
		material = m;
		normal = n;
		area = a;
		internal_bounces = i;
		reflected_rays = r;
		rays_entering = e;
		rays_leaving = l;
	}

	std::string material;
	Vector normal;
	double area;
	unsigned int internal_bounces;
	unsigned int reflected_rays;
	unsigned int rays_entering;
	unsigned int rays_leaving;
};

class ModelInfo{
public:
	ModelInfo(){
		num_less_than_45 = 0;
		num_between_45_and_90 = 0;
		num_between_90_and_135 = 0;
		num_greater_than_135 = 0;
		total_leaving = 0;
		num_leaving = 0;

		for(int d = 0; d < 180; d++){
			degree_totals.push_back(0);
		}
		for(int d = 0; d < 180; d++){
			degree_percentages.push_back(0);
		}
	}
	
	void AddFace(std::string material, Vector normal, double area, unsigned int internal_bounces, unsigned int reflected_rays, unsigned int rays_entering, unsigned int rays_leaving){
		faceInfo.push_back(Face(material, normal, area, internal_bounces, reflected_rays, rays_entering, rays_leaving));
	}

	void IncrementNumLeaving(){ num_leaving++; }

	void IncrementDegree(int location){
		if(location < 0 || location > 179){
			std::cout << "ERROR: trying to increment bad location" << std::endl;
			return;
		}
		degree_totals[location]++;
	}

	int GedDegreeInfo(int location){
		if(location < 0 || location > 179){
			std::cout << "ERROR: trying to get information on bad location" << std::endl;
			return 0;
		}
		return degree_totals[location];
	}

	double GedDegreePercentage(int location){
		if(location < 0 || location > 179){
			std::cout << "ERROR: trying to get information on bad location" << std::endl;
			return 0;
		}
		return degree_percentages[location];
	}

	void CalculatePercentages(){
		for(unsigned int d = 0; d < degree_totals.size(); d++){
			degree_percentages[d] = double(degree_totals[d])/double(num_leaving);
			//std::cout << "percentage: " << degree_percentages[d] << " degree num: " << degree_totals[d] << std::endl;
		}
	}
	
	const Face& getFace(unsigned int a){
		if( a < faceInfo.size()){
			return faceInfo[a];
		}
		else{
			//this shouldn't happen
			return faceInfo[0];
		}
	}
	unsigned int numFaces(){
		return faceInfo.size();
	}

	unsigned int num_less_than_45;
	unsigned int num_between_45_and_90;
	unsigned int num_between_90_and_135;
	unsigned int num_greater_than_135;
	unsigned int total_leaving;

	//std::vector<double> angles_leaving;
	unsigned int num_leaving; //same function as total leaving, but for degree totals
	std::vector<int> degree_totals;
	std::vector<double> degree_percentages;


private:
	std::vector<Face> faceInfo;
};

double AngleFromUP(const Vector& vec){

		double len = sqrt((vec.x*vec.x) + (vec.y*vec.y) + (vec.z*vec.z));

		//std::cout << "angle " << acos(vec.x/len) << std::endl;

		return acos(vec.y/len); 
	}

double AngleFromUPDegrees(const Vector& vec){

	return AngleFromUP(vec)*180/PI;
}

Vector ParseNormal(std::string normal){
	
	normal = normal.substr(1,normal.length() - 2); //remove {}s

	std::istringstream ss(normal);
	std::string token;

	std::getline(ss, token, ','); //x
	double x = std::stod(token);
	std::getline(ss, token, ','); //y
	double y = std::stod(token);
	std::getline(ss, token, ','); //z
	double z = std::stod(token);

	return Vector(x, y, z);

}

void CollectModelInfo(char* filename, ModelInfo& model){
	std::ifstream file(filename);
	std::string token;

	std::string material;
	Vector normal;
	double area;
	unsigned int internal_bounces;
	unsigned int reflected_rays;
	unsigned int rays_entering;
	unsigned int rays_leaving;

	bool data_time = false;
	//the first thing read in is the number of the face
	while (file >> token){
		//std::cout << token << ": ";
		//keep going until you hit the first face line
		if(token != "0:" && data_time == false) {
			continue;
		}

		//skip the line with the totals
		if(token[token.length()-1] != ':'){
			continue;
		}

		if(token != "leaving_directions:"){
			//Gathering per face photon leaving data
			data_time = true;

			file >> material;

			file >> token;
			normal = ParseNormal(token);
			
			file >> token; 
			area = std::stod(token.c_str());
			
			file >> token;
			internal_bounces = atoi(token.c_str());

			file >> token;
			reflected_rays = atoi(token.c_str());

			file >> token;
			rays_entering = atoi(token.c_str());

			file >> token;
			rays_leaving = atoi(token.c_str());

			model.AddFace(material, normal, area, internal_bounces, reflected_rays, rays_entering, rays_leaving);
		}
		

		//Getting all leaving angles
		if(token == "leaving_directions:"){
			double angle_from_top = 0;
			while(file >> token){
				normal = ParseNormal(token);
				angle_from_top = AngleFromUPDegrees(normal);
				//std::cout << angle_from_top << std::endl;
				//increment all of the numbers less than the angle for cumulative sum
				int degree = 180;
				int location = 179;
				model.IncrementNumLeaving();
				for(; angle_from_top < degree; degree--, location--){
					model.IncrementDegree(location);
				}
			}
			model.CalculatePercentages();
		}
	}

}

void PrintPerDegreeData(std::vector<ModelInfo>& modelInfo){
	for(int i = 0; i < 180; i++){
		std::cout << i +1 << " ";

		for(unsigned int m = 0; m < modelInfo.size(); m++){
			std::cout << modelInfo[m].GedDegreePercentage(i)  << " ";
		}
		std::cout << std::endl; 
	}
}

void Print4AngleData(std::vector<ModelInfo>& modelInfo){
	for(unsigned int m = 0; m < modelInfo.size(); m++){
		std::cout << "Model: " << m << std::endl;
		for(int f = 0; f < modelInfo[m].numFaces()-1; f ++){
			const Face face = modelInfo[m].getFace(f);
		//	std::cout << "    " << f << " " << face.material << " " << face.normal.x << " " << face.normal.y << " " 
		//						<< face.normal.z << " " << face.area << " " << face.internal_bounces << " " << face.reflected_rays 
		//						<< " " << face.rays_entering << " " << face.rays_leaving << std::endl;
		}

		for(int f = 0; f < modelInfo[m].numFaces()-1; f ++){
			if(modelInfo[m].getFace(f).material == "diamond"){
				if(AngleFromUP(modelInfo[m].getFace(f).normal) < PI/4.0){
					modelInfo[m].num_less_than_45 += modelInfo[m].getFace(f).rays_leaving;
				}
				else if (AngleFromUP(modelInfo[m].getFace(f).normal) < PI/2.0){
					modelInfo[m].num_between_45_and_90 += modelInfo[m].getFace(f).rays_leaving;
				}
				else if(AngleFromUP(modelInfo[m].getFace(f).normal) < (3.0 * PI/4.0)){
					modelInfo[m].num_between_90_and_135 += modelInfo[m].getFace(f).rays_leaving;
				}
				else{
					modelInfo[m].num_greater_than_135 += modelInfo[m].getFace(f).rays_leaving;
				}
			}
		}
		modelInfo[m].total_leaving = modelInfo[m].num_less_than_45 + modelInfo[m].num_between_45_and_90 + modelInfo[m].num_between_90_and_135 + modelInfo[m].num_greater_than_135; 

		std::cout << modelInfo[m].num_less_than_45/(double)modelInfo[m].total_leaving << " " << modelInfo[m].num_between_45_and_90/(double)modelInfo[m].total_leaving << " "
					<< modelInfo[m].num_between_90_and_135/(double)modelInfo[m].total_leaving << " " << modelInfo[m].num_greater_than_135/(double)modelInfo[m].total_leaving << std::endl;
	}
}

int main(int argc, char* argv[]){

	if(argc < 2){
		std::cout << "Please enter file names" << std::endl;
	}

	std::vector<ModelInfo> modelInfo;

	int f = 1;

	//determine which information to print
	bool print_per_degree_data = false;
	bool print_4_angle_data = false;

	while (argv[f][0] == '-'){
		std::string arg = argv[f];
		if (arg == "-degree_data"){
			print_per_degree_data = true;
		}
		if (arg == "-face_data"){
			print_4_angle_data = true;
		}		
		f++;
	}

	//Get model info for each file
	for(; f < argc; f++){
		ModelInfo model;
		CollectModelInfo(argv[f], model);
		modelInfo.push_back(model);
	}

	if(print_per_degree_data){
		PrintPerDegreeData(modelInfo);
	}

	if(print_4_angle_data){
		Print4AngleData(modelInfo);
	}


	return 0;
}