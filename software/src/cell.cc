#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/cell.h"

using namespace std;
using namespace NMSUPratt;

CMSU_Boltzmann *CMSU_BoltzmannCell::boltzmann=nullptr;
char *CMSU_BoltzmannCell::message=new char[CLog::CHARLENGTH];

CMSU_BoltzmannCell::CMSU_BoltzmannCell(double xminset,double xmaxset,double yminset,double ymaxset,double etaminset,double etamaxset){
	xmin=xminset; xmax=xmaxset; ymin=yminset; ymax=ymaxset; etamin=etaminset; etamax=etamaxset;	
	ireflection=0;
	creflection=nullptr;
}

void CMSU_BoltzmannCell::Print(){
	snprintf(message,CLog::CHARLENGTH,"___ CELL INFO _____\n");
	snprintf(message,CLog::CHARLENGTH,"%six=%d, iy=%d, ieta=%d, xmin=%g, xmax=%g, ymin=%g, ymax=%g, etamin=%g, etamax=%g\n",message, ix,iy,ieta,xmin,xmax,ymin,ymax,etamin,etamax);
	snprintf(message,CLog::CHARLENGTH,"%s%d parts in cell\n",message,int(partmap.size()));
	snprintf(message,CLog::CHARLENGTH,"%s---------------------\n",message);
	CLog::Info(message);
}
