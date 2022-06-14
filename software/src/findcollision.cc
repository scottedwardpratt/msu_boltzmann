#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/cell.h"
#include "msu_sampler/resonances.h"
#include "msu_commonutils/constants.h"

bool CMSU_Boltzmann::CheckKinematics(CMSUPart *part1,CMSUPart *part2,
	double &Minv2,double &pibsquared,double &taucoll){
	
	double p1dotp2=0.0,m1squared,m2squared,rsquared=0.0;
	double t1,t2,z1,z2,tau1,tau2,denom,p1dotr=0.0,p2dotr=0.0;
	double dtau1overm1,dtau2overm2;
	FourVector r;

	const int g[4]={1,-1,-1,-1};
	int alpha;
	bool bjtranslate=false;
	bool possible=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}

	m1squared=part1->msquared;
	m2squared=part2->msquared;
	
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=part1->r[alpha]-part2->r[alpha];
		rsquared+=g[alpha]*r[alpha]*r[alpha];
		p1dotp2+=g[alpha]*part1->p[alpha]*part2->p[alpha];
		p1dotr+=g[alpha]*part1->p[alpha]*r[alpha];
		p2dotr+=g[alpha]*part2->p[alpha]*r[alpha];
	}

	Minv2=m1squared+m2squared+2.0*p1dotp2;
		
	denom=p1dotp2*p1dotp2-m1squared*m2squared;

	pibsquared=PI*(-rsquared+(2.0*p1dotp2*p1dotr*p2dotr-p1dotr*p1dotr*m2squared-p2dotr*p2dotr*m1squared)/denom);


	dtau1overm1=(p1dotr*m2squared-p2dotr*p1dotp2)/denom;
	dtau2overm2=-(p2dotr*m1squared-p1dotr*p1dotp2)/denom;
	t1=part1->r[0]+dtau1overm1*part1->p[0];
	t2=part2->r[0]+dtau2overm2*part2->p[0];
	z1=part1->r[3]+dtau1overm1*part1->p[3];
	z2=part2->r[3]+dtau2overm2*part2->p[3];

	if(t1>0 && t2>0  && fabs(z1)<t1 && fabs(z2)<t2){
		tau1=sqrt(t1*t1-z1*z1);
		tau2=sqrt(t2*t2-z2*z2);
		taucoll=0.5*(tau1+tau2);
		if(taucoll>part1->tau0 && taucoll>part2->tau0){
			if(taucoll<part1->tauexit && taucoll<part2->tauexit)
				possible=true;
		}
	}

	if(bjtranslate)
		part1->BjorkenUnTranslate();
	return possible;
}

double CMSU_Boltzmann::GetSigma(CMSUPart *part1,CMSUPart *part2,double Minv2,
		double &sigma_scatter,double &sigma_merge,double &sigma_annihilation,double &sigma_inel,
		vector<double> &dsigma_merge){
	double sigmatot=0,M,Gamma,b,degenR,degen1,degen2,q2,dsigma,sf;
	int ir1,ir2,irflip;
	Cmerge *merge;
	
	if((part1->balanceID>=0 && part2->balanceID<0) || (part2->balanceID>=0 && part1->balanceID<0)){
		sigma_scatter=SIGMABF/double(NSAMPLE);
		sigmatot=sigma_scatter;
		return sigmatot;
	}
	else{
		sigma_scatter=SIGMADEFAULT/double(NSAMPLE);
	}

	if(BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0){
		sigma_annihilation=GetAnnihilationSigma(part1,part2);
	}
	else
		sigma_annihilation=0.0;

	if(SIGMAINELASTIC)
		sigma_inel=SIGMAINELASTIC/double(NSAMPLE);
	else
		sigma_inel=0.0;

	ir1=part1->resinfo->ires;
	ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1;
		ir1=ir2;
		ir2=irflip;
	}
	merge=reslist->MergeArray[ir1][ir2];
	dsigma_merge.clear();
	while(merge!=NULL){
		degen1=part1->resinfo->degen;
		degen2=part2->resinfo->degen;
		M=sqrt(Minv2);
		if(M>merge->resinfo->minmass){
			q2=Misc::triangle2(Minv2,part1->msquared,part2->msquared);
			sf=merge->resinfo->GetSpectralFunction(M);
			b=merge->branching;
			Gamma=merge->resinfo->width;
			degenR=merge->resinfo->degen;
			dsigma=b*sf*Gamma*2.0*PI*PI*(degenR/(degen1*degen2))*HBARC_GEV*HBARC_GEV/q2;
			if(dsigma>40.0)
				dsigma=40.0;
			dsigma=dsigma/double(NSAMPLE);
		}
		else
			dsigma=0.0;
		dsigma_merge.push_back(dsigma);
		sigma_merge+=dsigma;
		merge=merge->next;
	}
	sigmatot=sigma_scatter+sigma_inel+sigma_annihilation+sigma_merge;
	/*
	if(part1->resinfo->pid==-211 && part2->resinfo->pid==211)
		printf("m1=%g, m2=%g, Minv=%g, q=%g, sigmatot=%g\n",
			part1->resinfo->mass,part2->resinfo->mass,sqrt(Minv2),sqrt(q2),10.0*sigmatot);
		*/	
	return sigmatot;
}

bool CMSU_Boltzmann::FindCollision(CMSUPart *part1,CMSUPart *part2,double &taucoll){
	if((part1->balanceID>=0) && (part2->balanceID>=0)){
		return false;
	}
	bool checkpossible;
	double sigma_scatter=0.0,sigma_merge=0.0,sigma_annihilation=0.0,sigma_inel=0.0;
	vector<double> dsigma_merge;
	bool collide=false;
	double pibsquared,Minv2,sigmatot;
	dsigma_merge.clear();

	checkpossible=CheckKinematics(part1,part2,Minv2,pibsquared,taucoll);
	if(checkpossible){
		sigmatot=GetSigma(part1,part2,Minv2,sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel,
			dsigma_merge);
		if(pibsquared<sigmatot){
			AddAction_Collision(part1,part2,taucoll,pibsquared,
						sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel,
						dsigma_merge);
			collide=true;
		}
	}
	return collide;
}		

void CMSU_Boltzmann::FindAllCollisions(){
	double taucoll;
	CMSUPartMap::iterator ppos1,ppos2;
	CMSUPart *part1,*part2;
	CActionMap::iterator epos;
	for(ppos1=PartMap.begin();ppos1!=PartMap.end();++ppos1){
		part1=ppos1->second;
		part1->KillActions();
		part1->active=true;
		part1->ChangeCell(part1->FindCell());
		if(part1->cell!=NULL){
			part1->FindCellExit();
		}
		if(part1->resinfo->decay)
			part1->FindDecay();
	}
	ppos2=PartMap.begin();
	++ppos2;
	while(ppos2!=PartMap.end()){
		part2=ppos2->second;
		ppos1=PartMap.begin();
		while(ppos1!=ppos2){
			part1=ppos1->second;
			if(part1->balanceID<0 || part2->balanceID<0){
				if(part1->cell!=NULL && part2->cell!=NULL){
					FindCollision(part1,part2,taucoll);
				}
			}
			++ppos1;
		}
		++ppos2;
	}
}
