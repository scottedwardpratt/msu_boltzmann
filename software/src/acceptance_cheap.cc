//
//  acceptance_CHEAP.cc
//
#include "msu_boltzmann/acceptance.h"
#include "msu_boltzmann/msupart.h"
#include "msu_sampler/resonances.h"
#include "msu_commonutils/parametermap.h"

using namespace std;

CAcceptance_CHEAP::CAcceptance_CHEAP(CparameterMap *parmapin) : CAcceptance(){
	ETAMAX=parmapin->getD("MSU_BOLTZMANN_ETAMAX",6.0);
	ETAMIN=-ETAMAX;
	ptmin=0.2;
	ptmax=20000000.0;
}

void CAcceptance_CHEAP::CalcAcceptance(bool &accept,double &efficiency,CMSUPart *part){
	double pt;
	double eta,pmag;
	
	accept=false;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	
	if(pt<10000000 && fabs(eta)<ETAMAX){
		accept=true;
		efficiency=1.0;
	}
	else{
		accept=false;
		efficiency=0.0;
	}
}

void CAcceptance_CHEAP::CalcAcceptanceNoID(bool &accept,double &efficiency,CMSUPart *part){
	double pt,eta,pmag;
	accept=false;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	
	if(pt<10000000 && fabs(eta)<ETAMAX){
		accept=true;
		efficiency=1.0;
	}
	else{
		accept=false;
		efficiency=0.0;
	}
	
	
}
