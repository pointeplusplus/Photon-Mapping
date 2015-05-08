#include "kdtree.h"

#define MAX_PHOTONS_BEFORE_SPLIT 100
#define MAX_DEPTH 15

// Special constructor to create a balanced KDTree.
KDTree::KDTree(
		const BoundingBox &_bbox, 
		int _depth, 
		const std::vector<Photon> & _photons) 
{
	bbox = _bbox;
	depth = _depth;
	photons = _photons;
	child1 = NULL;
	child2 = NULL;
	
	if(depth < MAX_DEPTH && photons.size() > MAX_PHOTONS_BEFORE_SPLIT)
	{
		BoundingBox b1;
		BoundingBox b2;
		findBalancedChildBoundingBox(b1, b2);
		std::vector<Photon> v1;
		v1.insert(v1.end(), photons.begin(), photons.begin() + photons.size() /2);
		child1 = new KDTree(b1, depth + 1, v1);
		
		std::vector<Photon> v2;
		v2.insert(v2.end(), photons.begin() + photons.size() /2, photons.end());
		child2 = new KDTree(b2, depth + 1, v2);
		photons.clear();
	}
}




// ==================================================================
// DESTRUCTOR
// ==================================================================
KDTree::~KDTree() 
{
	if (!isLeaf()) 
	{
		delete child1;
		delete child2;
	} 
	else 
	{
		// delete all the photons (this is done automatically since they
		// are stored directly in an STL vector, not using pointers)
	}
}


// ==================================================================
// HELPER FUNCTIONS

bool KDTree::PhotonInCell(const Photon &p, const BoundingBox & bb) 
{
  const Vec3f& min = bb.getMin();
  const Vec3f& max = bb.getMax();
  const Vec3f &position = p.getPosition();
  if (position.x() > min.x() - EPSILON &&
		position.y() > min.y() - EPSILON &&
		position.z() > min.z() - EPSILON &&
		position.x() < max.x() + EPSILON &&
		position.y() < max.y() + EPSILON &&
		position.z() < max.z() + EPSILON) 
	{
		return true;
	}
  return false;
}

bool KDTree::overlaps(const BoundingBox &bb) const 
{
	const Vec3f& bb_min = bb.getMin();
	const Vec3f& bb_max = bb.getMax();
	const Vec3f& tmp_min = bbox.getMin();
	const Vec3f& tmp_max = bbox.getMax();
	if (bb_min.x() > tmp_max.x()) 
	{
		return false;
	}
	if (tmp_min.x() > bb_max.x()) 
	{
		return false;
	}
	if (bb_min.y() > tmp_max.y()) 
	{
		return false;
	}
	if (tmp_min.y() > bb_max.y()) 
	{
		return false;
	}
	if (bb_min.z() > tmp_max.z()) 
	{
		return false;
	}
	if (tmp_min.z() > bb_max.z()) 
	{
		return false;
	}
	return true;
}


// ==================================================================
void KDTree::AddPhoton(const Photon &p) 
{
	const Vec3f &position = p.getPosition();
	if(!PhotonInCell(p, bbox)) 
	{
		return; 
	}
	if (isLeaf()) 
	{
		// this cell is a leaf node
		photons.push_back(p);
		if (photons.size() > MAX_PHOTONS_BEFORE_SPLIT && depth < MAX_DEPTH) {
		  SplitCell();
		}
	} 
	else 
	{
		// this cell is not a leaf node
		// decide which subnode to recurse into
		if (split_axis == 0)
		{
			if (position.x() < split_value)
			{
				child1->AddPhoton(p);
			}
			else
			{
				child2->AddPhoton(p);
			}
		} 
		else if (split_axis == 1) 
		{
			if (position.y() < split_value)
			{
				child1->AddPhoton(p);
			}
			else
			{
				child2->AddPhoton(p);
			}
		} 
		else 
		{
			assert (split_axis == 2);
			if (position.z() < split_value)
			{
				child1->AddPhoton(p);
			}
		  else
		  {
			child2->AddPhoton(p);
		  }
		}
	}
}


// ==================================================================
void KDTree::CollectPhotonsInBox(
	const BoundingBox &bb, 
	std::vector<Photon> &photons) const 
{
	// explicitly store the queue of cells that must be checked (rather
	// than write a recursive function)
	std::vector<const KDTree*> todo;  
	todo.push_back(this);
	while (!todo.empty()) 
	{
		const KDTree *node = todo.back();
		todo.pop_back(); 
		if (!node->overlaps(bb))
		{
			continue;
		}
		if (node->isLeaf()) 
		{
			// if this cell overlaps & is a leaf, add all of the photons into the master list
			// NOTE: these photons may not be inside of the query bounding box
			const std::vector<Photon> &photons2 = node->getPhotons();
			int num_photons = photons2.size();
			for (int i = 0; i < num_photons; i++) 
			{
				photons.push_back(photons2[i]);
			}
		} 
		else 
		{
			// if this cell is not a leaf, explore both children
			todo.push_back(node->getChild1());
			todo.push_back(node->getChild2());
		} 
	}
}

void KDTree::findChildBoundingBox(BoundingBox & child1, BoundingBox & child2)
{
	const Vec3f& min = bbox.getMin();
	const Vec3f& max = bbox.getMax();
	double dx = max.x()-min.x();
	double dy = max.y()-min.y();
	double dz = max.z()-min.z();
	// split this cell in the middle of the longest axis
	Vec3f min1,min2,max1,max2;
	if (dx >= dy && dx >= dz) 
	{
		split_axis = 0;
		split_value = min.x()+dx/2.0;
		min1 = Vec3f(min.x(), min.y(), min.z());
		max1 = Vec3f(split_value, max.y(), max.z());
		min2 = Vec3f(split_value, min.y(), min.z());
		max2 = Vec3f(max.x(), max.y(), max.z());
	} 
	else if (dy >= dx && dy >= dz) 
	{
		split_axis = 1;
		split_value = min.y()+dy/2.0;
		min1 = Vec3f(min.x(),min.y()    ,min.z());
		max1 = Vec3f(max.x(),split_value,max.z());
		min2 = Vec3f(min.x(),split_value,min.z());
		max2 = Vec3f(max.x(),max.y()    ,max.z());
	} 
	else 
	{
		assert (dz >= dx && dz >= dy);
		split_axis = 2;
		split_value = min.z()+dz/2.0;
		min1 = Vec3f(min.x(),min.y(),min.z()    );
		max1 = Vec3f(max.x(),max.y(),split_value);
		min2 = Vec3f(min.x(),min.y(),split_value);
		max2 = Vec3f(max.x(),max.y(),max.z()    );
	}
	
	child1 = BoundingBox(min1,max1);
	child2 = BoundingBox(min2,max2);
}

bool sortByX(const Photon & p1, const Photon & p2)
{
	return p1.getPosition().x() < p2.getPosition().x();
}

bool sortByY(const Photon & p1, const Photon & p2)
{
	return p1.getPosition().y() < p2.getPosition().y();
}

bool sortByZ(const Photon & p1, const Photon & p2)
{
	return p1.getPosition().z() < p2.getPosition().z();
}

void KDTree::findBalancedChildBoundingBox(BoundingBox & child1, BoundingBox & child2)
{
	const Vec3f& min = bbox.getMin();
	const Vec3f& max = bbox.getMax();
	double dx = max.x()-min.x();
	double dy = max.y()-min.y();
	double dz = max.z()-min.z();
	// split this cell in the middle of the longest axis
	Vec3f min1,min2,max1,max2;
	if (dx >= dy && dx >= dz) 
	{
		split_axis = 0;
		
		//Calculate split value
		std::sort(photons.begin(), photons.end(), sortByX);
		split_value = (photons[photons.size()/2]).getPosition().x();
		
		min1 = Vec3f(min.x(), min.y(), min.z());
		max1 = Vec3f(split_value, max.y(), max.z());
		min2 = Vec3f(split_value, min.y(), min.z());
		max2 = Vec3f(max.x(), max.y(), max.z());
	} 
	else if (dy >= dx && dy >= dz) 
	{
		
		split_axis = 1;
		
		//Calculate split value
		std::sort(photons.begin(), photons.end(), sortByY);
		split_value = (photons[photons.size()/2]).getPosition().y();
		
		min1 = Vec3f(min.x(),min.y()    ,min.z());
		max1 = Vec3f(max.x(),split_value,max.z());
		min2 = Vec3f(min.x(),split_value,min.z());
		max2 = Vec3f(max.x(),max.y()    ,max.z());
	} 
	else 
	{
		assert (dz >= dx && dz >= dy);
		split_axis = 2;
		//Calculate split value
		std::sort(photons.begin(), photons.end(), sortByZ);
		split_value = (photons[photons.size()/2]).getPosition().z();
		
		min1 = Vec3f(min.x(),min.y(),min.z()    );
		max1 = Vec3f(max.x(),max.y(),split_value);
		min2 = Vec3f(min.x(),min.y(),split_value);
		max2 = Vec3f(max.x(),max.y(),max.z()    );
	}
	
	child1 = BoundingBox(min1,max1);
	child2 = BoundingBox(min2,max2);
}

// ==================================================================
void KDTree::SplitCell() 
{
	BoundingBox c1;
	BoundingBox c2;
	findChildBoundingBox(c1, c2);
	
	// create two new children
	child1 = new KDTree(c1, depth+1);
	child2 = new KDTree(c2, depth+1);
	int num_photons = photons.size();
	std::vector<Photon> tmp = photons;
	photons.clear();
	// add all the photons to one of those children
	for (int i = 0; i < num_photons; i++) 
	{
		const Photon &p = tmp[i];
		this->AddPhoton(p);
	}
}

// ==================================================================

void KDTree::printLeaves(int level, std::vector<int> & nums) const
{
	if(isLeaf())
	{
		std::cout << "Leaf is at level " << level << " and has " 
			<< photons.size() << " photons.\n";
		nums.push_back(photons.size());
		return;
	}
	else
	{
		if(child1 != NULL)
		{
			child1->printLeaves(level + 1, nums);
		}
		if(child2 != NULL)
		{
			child2->printLeaves(level + 1, nums);
		}
	}
}

