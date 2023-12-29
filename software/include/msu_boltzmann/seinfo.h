#ifndef __SEINFO_H__
#define __SEINFO_H__

#include "msu_commonutils/commondefs.h"

using namespace std;
namespace NMSUPratt;
using namespace NMSUPratt;

class CMSU_Boltzmann;

class CSEInfo{
public:
	CSEInfo(CMSU_Boltzmann *boltzmannset);
	// vectors hold information for different times
	vector<double> Pbar,Tzz,epsilon,nhadrons,K0,F0,uperpbar;
	// this is information for latest time
	void Zero(); // sets arrays to zero
	void SECalc();
	void Print();
	double R,TAU0,TAUMAX,DELTAU,ETAOVERS,RMAX;
	int NTAU;
	int NETEVENTS;
	CMSU_Boltzmann *boltzmann;
	char message[500];
};

#endif