#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/randy.h"
#include "msu_sampler/resonances.h"
#include "msu_boltzmann/msupart.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"

void CMSU_Boltzmann::GetMassesForDecay(vector<double> &mass,int nbodies,array<CMSUPart *,5> &daughter){
	double minmass_daughters=0.0,mass_daughters=0.0,width_daughters=0.0;
	vector<double> probmax;
	double mtot,E0,netprob,Estarmax;
	int imass,ibody;
	int ntry;
	for(ibody=0;ibody<nbodies;ibody++){
		minmass_daughters+=daughter[ibody]->resinfo->minmass;
		mass_daughters+=daughter[ibody]->resinfo->mass;
		width_daughters+=daughter[ibody]->resinfo->width;
	}
	if(mass[0]>minmass_daughters+0.02 && mass[0]>mass_daughters-1.5*width_daughters){
		ntry=0;
		do{
			mtot=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				netprob=randy->ran();
				mass[ibody+1]=daughter[ibody]->resinfo->GenerateMassFromSF(netprob);
				mtot+=mass[ibody+1];
			}
			ntry+=1;
		}while(mtot>mass[0]);
		if(ntry>10000){
			printf("A: ntry=%d\n",ntry);
			printf("mass[0]=%g, minmass_daughters=%g, mass_daughters=%g\n",mass[0],minmass_daughters, mass_daughters);
		}
	}
	else{
		probmax.resize(nbodies);
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
		do{
			mtot=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				if(daughter[ibody]->resinfo->decay){
					netprob=randy->ran()*probmax[ibody];
					mass[ibody+1]=daughter[ibody]->resinfo->GenerateMassFromSF(netprob);
				}
				else
					mass[ibody+1]=daughter[ibody]->resinfo->mass;
				mtot+=mass[ibody+1];
			}
		}while(mtot>mass[0]);
	}
	double mcheck=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		mcheck+=mass[ibody+1];
	}
	if(mcheck>mass[0]){
		printf("CMSU_Boltzmann::GetMassesForDecay -- masses too big %g < %g\n",mass[0],mcheck);
		exit(1);
	}
}

void CMSU_Boltzmann::Decay(CMSUPart *mother,int &nbodies,array<CMSUPart *,5> &daughter){
	int ibody,alpha;
	CMSUPart *dptr;
	vector<double> mass(6);
	vector<FourVector> p(5);
	FourVector u,pprime;

	mass[0]=mother->GetMass();

	GetMassesForDecay(mass,nbodies,daughter);
	
	/* Make masses for daughters */
	
	

	/*
	int ntry=0;
	do{
		printf("-----------------------------------\n");
		ntry+=1;
		mtot=0.0;
		for(ibody=0;ibody<nbodies;ibody++){
			if(daughter[ibody]->resinfo->decay){
				mass[ibody+1]=daughter[ibody]->resinfo->GenerateMassForDecay(Estarmax);
			}
			else
				mass[ibody+1]=daughter[ibody]->resinfo->mass;
			mtot+=mass[ibody+1];
		}
		printf("mtot=%g, mothermass=%g, minmass=%g, minmass_daugthers=%g\n",
				mtot,mass[0],mother->resinfo->minmass,minmass_daughters);
		if(ntry>100){
			printf("mother info\n");
			mother->Print();
			mother->resinfo->Print();
			//mother->resinfo->PrintSpectralFunction();
			exit(1);
		}
	}while(mtot>mass[0]);
	*/
	
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
		Misc::lorentz(u,p[ibody],pprime);
		dptr->active=true;
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
		dptr->SetMass();
		//dptr->Setp0();
		dptr->SetY();
		//dptr->eta0=mother->eta0;
		dptr->phi0=mother->phi0;
	}
}
