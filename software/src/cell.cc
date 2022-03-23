#include "boltzmann.h"
#include "cell.h"

CMSU_Boltzmann *CMSU_BoltzmannCell::boltzmann=NULL;
char *CMSU_BoltzmannCell::message=new char[500];

CMSU_BoltzmannCell::CMSU_BoltzmannCell(double xminset,double xmaxset,double yminset,double ymaxset,double etaminset,double etamaxset){
	xmin=xminset; xmax=xmaxset; ymin=yminset; ymax=ymaxset; etamin=etaminset; etamax=etamaxset;	
	ireflection=0;
	creflection=NULL;
}

void CMSU_BoltzmannCell::Print(){
	sprintf(message,"___ CELL INFO _____\n");
	sprintf(message,"%six=%d, iy=%d, ieta=%d, xmin=%g, xmax=%g, ymin=%g, ymax=%g, etamin=%g, etamax=%g\n",message, ix,iy,ieta,xmin,xmax,ymin,ymax,etamin,etamax);
	sprintf(message,"%s%d parts in cell\n",message,int(partmap.size()));
	sprintf(message,"%s---------------------\n",message);
	CLog::Info(message);
}
