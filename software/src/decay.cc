#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/randy.h"
#include "msu_eos/resonances.h"
#include "msu_boltzmann/msupart.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"
#include "msu_boltzmann/decay.h"
#include "msu_commonutils/decay_nbody.h"

using namespace std;
using namespace NMSUPratt;

CMSU_Decay::CMSU_Decay(Crandy *randyset){
	randy=randyset;
	mass.resize(6);
	p.resize(5);
	probmax.resize(5);
	decay_nbody=new CDecay_NBody(randy);
}

void CMSU_Decay::GetMassesForDecay(vector<double> &mass,int nbodies,array<CMSUPart *,5> &daughter){
	double minmass_daughters=0.0,mass_daughters=0.0,width_daughters=0.0;
	double mtot,E0,netprob,Estarmax;
	int imass,ibody;
	int ntry;
	for(ibody=0;ibody<nbodies;ibody++){
		minmass_daughters+=daughter[ibody]->resinfo->minmass;
		mass_daughters+=daughter[ibody]->resinfo->mass;
		width_daughters+=daughter[ibody]->resinfo->width;
	}
	if(mass_daughters<minmass_daughters || mass[0]<minmass_daughters){
		snprintf(message,CLog::CHARLENGTH,"mothermass=%g, daughter masses=(",mass[0]);
		for(int i=1;i<=nbodies;i++){
			snprintf(message,CLog::CHARLENGTH,"%s%g, ",message,mass[i]);
		}
		snprintf(message,CLog::CHARLENGTH,"%s)\n",message);
		CLog::Info(message);
		CLog::Fatal("masses out of whack in decay\n");
	}
	if(mass[0]>minmass_daughters+0.5*(mass_daughters-minmass_daughters) && mass[0]>mass_daughters-width_daughters){
		ntry=0;
		do{
			mtot=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				netprob=randy->ran();
				mass[ibody+1]=daughter[ibody]->resinfo->GenerateMassFromSF(netprob);
				mtot+=mass[ibody+1];
			}
			ntry+=1;
			if(ntry>20){
				mass[0]+=0.005;
				if(ntry>200)
					CLog::Info("In CMSU_Decay::GetMassesForDecay, ntry="+to_string(ntry)+"\n");
			}
		}while(mtot>mass[0]);
	}
	else{
		Estarmax=mass[0]-minmass_daughters;
		for(ibody=0;ibody<nbodies;ibody++){
			if(daughter[ibody]->resinfo->decay){
				probmax[ibody]=0.0;
				E0=daughter[ibody]->resinfo->SpectEVec[0];
				imass=0;
				while(daughter[ibody]->resinfo->SpectEVec[imass]-E0<Estarmax && imass<daughter[ibody]->resinfo->NSPECTRAL){
					probmax[ibody]+=daughter[ibody]->resinfo->SpectVec[imass];
					imass+=1;
				}
			}
		}
		ntry=0;
		do{
			mtot=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				if(daughter[ibody]->resinfo->decay){
					do{
						netprob=randy->ran()*probmax[ibody];
					}while(netprob>1.0);
					mass[ibody+1]=daughter[ibody]->resinfo->GenerateMassFromSF(netprob);
				}
				else
					mass[ibody+1]=daughter[ibody]->resinfo->mass;
				mtot+=mass[ibody+1];
			}
			ntry+=1;
			if(ntry>20){
				mass[0]+=0.005;
				if(ntry>200)
					CLog::Info("ntry="+to_string(ntry)+"\n");
			}
		}while(mtot>mass[0]); 
	}
	double mcheck=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		mcheck+=mass[ibody+1];
	}
	if(mcheck>mass[0]){
		snprintf(message,CLog::CHARLENGTH,"CMSU_Boltzmann::GetMassesForDecay -- masses too big %g < %g\n",mass[0],mcheck);
		CLog::Fatal(message);
	}
}

void CMSU_Decay::Decay(CMSUPart *mother,int &nbodies,array<CMSUPart *,5> &daughter){
	int ibody,alpha;
	CMSUPart *dptr;
	double tau=boltzmann->tau;
	mass[0]=sqrt(mother->msquared);

	GetMassesForDecay(mass,nbodies,daughter);
	
	if(nbodies==2){
		decay_nbody->SetMasses2(mass[0],mass[1],mass[2]);
		decay_nbody->GenerateMomenta2(p[0],p[1]);
	}
	else if(nbodies==3){
		decay_nbody->SetMasses3(mass[0],mass[1],mass[2],mass[3]);
		decay_nbody->GenerateMomenta3(p[0],p[1],p[2]);
	}
	else if(nbodies==4){
		decay_nbody->SetMasses4(mass[0],mass[1],mass[2],mass[3],mass[4]);
		decay_nbody->GenerateMomenta4(p[0],p[1],p[2],p[3]);
	}
	else{
		decay_nbody->SetMasses(nbodies,mass);
		decay_nbody->GenerateMomenta(p);
	}
	
	/* Boost the new particles */
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=mother->p[alpha]/mass[0];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::Boost(u,p[ibody],pprime);
		dptr->active=true;
		dptr->balanceID=mother->balanceID;
		for(alpha=0;alpha<4;alpha++)
			dptr->p[alpha]=pprime[alpha];
		//dptr->CopyPositionInfo(mother);
		dptr->r[1]=mother->r[1];
		dptr->r[2]=mother->r[2];
		dptr->eta=mother->eta;
		dptr->tau0=tau;
		dptr->r[0]=tau*cosh(dptr->eta);
		dptr->r[3]=tau*sinh(dptr->eta);
		//dptr->msquared=pow(dptr->resinfo->mass,2);
		dptr->msquared=mass[ibody+1]*mass[ibody+1];
		dptr->Setp0();
		//dptr->Setp0();
		dptr->SetY();
		//dptr->eta0=mother->eta0;
		dptr->phi0=mother->phi0;
	}
	p.clear();
	mass.clear();
	
}
