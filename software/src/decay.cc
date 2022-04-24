#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/randy.h"
#include "msu_sampler/resonances.h"
#include "msu_boltzmann/msupart.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"

void CMSU_Boltzmann::Decay(CMSUPart *mother,int &nbodies,array<CMSUPart *,5> &daughter){
	int ibody,alpha;
	double mtot,Estarmax;
	CMSUPart *dptr;
	vector<double> mass(6);
	vector<FourVector> p(5);
	FourVector u,pprime;

	mass[0]=mother->GetMass();
	
	/* Make masses for daughters */
	
	double minmass_daughters=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		minmass_daughters+=daughter[ibody]->resinfo->minmass;
	}
	Estarmax=mass[0]-minmass_daughters;
	printf("check in, mothermass=%g, minmass_daughters=%g, Estarmax=%g\n",
		mass[0],minmass_daughters,Estarmax);

	Estarmax=mass[0]-minmass_daughters;
	
	do{
		mtot=0.0;
		for(ibody=0;ibody<nbodies;ibody++){
			if(daughter[ibody]->resinfo->decay){
				mass[ibody+1]=daughter[ibody]->resinfo->GenerateMass_T0(Estarmax);
			}
			else
				mass[ibody+1]=daughter[ibody]->resinfo->mass;
			mtot+=mass[ibody+1];
		}
		printf("mtot=%g, mothermass=%g, minmass=%g, minmass_daugthers=%g\n",
			mtot,mass[0],mother->resinfo->minmass,minmass_daughters);
	}while(mtot>mass[0]);
	printf("check out\n");
	
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
