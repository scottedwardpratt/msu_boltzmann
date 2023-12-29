#ifndef __MSU_DECAY_H__
#define __MSU_DECAY_H__

#include "msu_boltzmann.h"

using namespace std;

namespace NMSUPratt{

class CMSU_Decay{
public:
	vector<double> mass,probmax;
	vector<FourVector> p;
	FourVector u,pprime;
	char message[200];
	CDecay_NBody *decay_nbody;
	Crandy *randy;
	CMSU_Decay(Crandy *randy);
	void GetMassesForDecay(vector<double> &mass,int nbodies,array<CMSUPart *,5> &daughter);
	void Decay(CMSUPart *mother,int &nbodies,array<CMSUPart *,5> &daughter);
	CMSU_Boltzmann *boltzmann;
};

}

#endif
