#ifndef __CELL_H__
#define __CELL_H__

#include "commondefs.h"

using namespace std;

//!A cell in the expanding cell mesh
/*!
\version 1.0
\author Scott Pratt
\date March 2011

In the CMSU_Boltzmann model, the model space is expressed as a mesh grid of cells that expand as time propogates. The mesh is populated by cells, which are CMSU_BoltzmannCell objects. This class keeps track of the particles populating it (in a particle map), as well as its spatial dimensions and neighbors. The neighbors are especially relevant, as actions such as collisions are scheduled by checking against particles inside its the current cell, as well as all neighboring cells.
*/

class CMSU_BoltzmannCell{
public:
	class CMSU_BoltzmannCell *neighbor[3][3][3];
	int ix,iy,ieta;
	class CMSU_BoltzmannCell *creflection;
	int ireflection;
	double xmin,xmax,ymin,ymax,etamin,etamax;
	CMSUPartMap partmap;
	void PrintPartMap(CMSUPartMap *partmap);
	void KillAllParts();
	void ReKeyAllParts();
	void Print();
	vector<double> dens;
	CMSU_BoltzmannCell(double xmin,double xmax,double ymin,double ymax,double etamin,double etamax);
	static CMSU_Boltzmann *boltzmann;
	static char *message;
};

//!A cell in the expanding cell mesh
/*!
\version 1.0
\author Scott Pratt
\date March 2011

In the CMSU_Boltzmann model, the model space is expressed as a mesh grid of cells that expand as time propogates. The mesh is populated by cells, which are CMSU_BoltzmannCell objects. This class keeps track of the particles populating it (in a particle map), as well as its spatial dimensions and neighbors. The neighbors are especially relevant, as actions such as collisions are scheduled by checking against particles inside its the current cell, as well as all neighboring cells.
*/

#endif