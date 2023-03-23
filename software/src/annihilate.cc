#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/cell.h"
#include "msu_eos/resonances.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"


int CMSU_Boltzmann::Annihilate(CMSUPart *part1,CMSUPart *part2,int &ndaughters,array<CMSUPart*,5> &daughter){
	CMSUPart *dptr;
	int i,alpha,ibody,netq,nets,nK0bar,nK0,nKplus,nKminus,npi0,npiplus,npiminus;
	int nu,nubar,nd,ndbar,ns,nsbar,nbodies=5;
	FourVector u;
	double etabar,rbar[4]={0.0};
	CMSU_BoltzmannCell *newcell;
	FourVector P,pprime;
	//
	ndaughters=nbodies;
	vector<double> mass(6);
	vector<FourVector> p(5);

	CresInfo *resinfo1=part1->resinfo,*resinfo2=part2->resinfo,*resinfo;
	if(resinfo1->baryon<0){
		resinfo=resinfo2;
		resinfo2=resinfo1;
		resinfo1=resinfo;
	}
	netq=resinfo1->charge+resinfo2->charge;
	nets=resinfo1->strange+resinfo2->strange;
	// resinfo 1 = baryon, resinfo2 = antibaryon
	nu=resinfo1->charge+1;
	ns=-resinfo1->strange;
	nd=3-nu-ns;
	nubar=-resinfo2->charge+1;
	nsbar=resinfo2->strange;
	ndbar=3-nubar-nsbar;

	vector<int> quark(nbodies);
	vector<int> antiq(nbodies);

	for(i=0;i<nu;i++)
		quark[i]=1;
	for(i=nu;i<nu+nd;i++)
		quark[i]=-1;
	for(i=nu+nd;i<3;i++){
		quark[i]=0;
	}

	for(i=nbodies-1;i>nbodies-1-nubar;i--)
		antiq[i]=1;
	for(i=nbodies-1-nubar;i>nbodies-1-nubar-ndbar;i--)
		antiq[i]=-1;
	for(i=nbodies-1-nubar-ndbar;i>1;i--)
		antiq[i]=0;

	if(randy->ran()<0.5){
		quark[3]=1;
		antiq[1]=1;
	}
	else{
		quark[3]=-1;
		antiq[1]=-1;
	}
	if(randy->ran()<0.5){
		quark[nbodies-1]=1;
		antiq[0]=1;
	}
	else{
		quark[nbodies-1]=-1;
		antiq[0]=-1;
	}
	nKplus=nKminus=nK0bar=nK0=npiplus=npiminus=npi0=0;
	for(i=0;i<nbodies;i++){
		if(quark[i]==1 && antiq[i]==1){
			npi0+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(111);
		}
		else if(quark[i]==1 && antiq[i]==-1){
			npiplus+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(211);
		}
		else if(quark[i]==-1 && antiq[i]==1){
			npiminus+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(-211);
		}
		else if(quark[i]==-1 && antiq[i]==-1){
			npi0+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(111);
		}
		else if(quark[i]==0 && antiq[i]==0){
			npi0+=1; // only happens for Omega-AntiOmega annihilation
			daughter[i]->resinfo=reslist->GetResInfoPtr(111);
		}
		else if(quark[i]==1 && antiq[i]==0){
			nKplus+=1;
			daughter[i]->resinfo=daughter[i]->resinfo=reslist->GetResInfoPtr(321);
		}
		else if(quark[i]==0 && antiq[i]==1){
			nKminus+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(-321);
		}
		else if(quark[i]==0 && antiq[i]==-1){
			nK0bar+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(311);
		}
		else if(quark[i]==-1 && antiq[i]==0){
			nK0+=1;
			daughter[i]->resinfo=reslist->GetResInfoPtr(-311);
		}
	}

	if(netq!= nKplus+npiplus-nKminus-npiminus){
		snprintf(message,sizeof(message),"charges don't add up, netq=%d\n",netq);
		CLog::Fatal(message);
	}
	if(nets!= nKplus+nK0-nKminus-nK0bar){
		snprintf(message,sizeof(message),"charges don't add up, nets=%d\n",nets);
		CLog::Fatal(message);
	}

	for(i=1;i<=nbodies;i++)
		mass[i]=daughter[i-1]->resinfo->mass;
	mass[0]=0.0;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=part1->p[alpha]+part2->p[alpha];
	}
	mass[0]=sqrt(P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3]);

	msudecay->decay_nbody->SetMasses(nbodies,mass);
	msudecay->decay_nbody->GenerateMomenta(p);
	
	rbar[1]=0.5*(part1->r[1]+part2->r[1]);
	rbar[2]=0.5*(part1->r[2]+part2->r[2]);
	etabar=0.5*(part1->eta+part2->eta);
	rbar[0]=tau*cosh(etabar);
	rbar[3]=tau*sinh(etabar);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=P[alpha]/mass[0];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::Boost(u,p[ibody],pprime);
		dptr->active=true;
		for(alpha=0;alpha<4;alpha++){
			dptr->p[alpha]=pprime[alpha];
			dptr->r[alpha]=rbar[alpha];
		}
		//dptr->msquared=pow(dptr->resinfo->mass,2);
		dptr->msquared=dptr->resinfo->mass*dptr->resinfo->mass;
		dptr->Setp0();
		dptr->SetY();
		dptr->tau0=tau;
		dptr->eta=etabar;
		if(fabs(dptr->eta)>ETAMAX)
			dptr->CyclicReset();
		//dptr->eta0=etabar;
		dptr->phi0=atan2(dptr->r[2],dptr->r[1]);
		dptr->active=true;
		newcell=dptr->FindCell();
		dptr->ChangeCell(newcell);

		dptr->ChangeMap(&PartMap);
		if(fabs(dptr->eta)>ETAMAX){
			snprintf(message,sizeof(message),"eta out of range in Annihilate\n");
			CLog::Fatal(message);
		}
		if(dptr->p[0]<0.0){
			snprintf(message,sizeof(message),"dptr->p[0]=%g\n",dptr->p[0]);
			CLog::Fatal(message);
		}
	}
	mass.clear();
	p.clear();
	return ndaughters;
}

double CMSU_Boltzmann::GetAnnihilationSigma(CMSUPart *part1,CMSUPart *part2){
	const double g[4]={1,-1,-1,-1};
	double Plab,p1dotp2,triangle,sigma_annihilation,rstrange,m1squared,m2squared;
	int alpha;
	m1squared=part1->p[0]*part1->p[0]-part1->p[1]*part1->p[1]-part1->p[2]*part1->p[2]
	-part1->p[3]*part1->p[3];
	m2squared=part2->p[0]*part2->p[0]-part2->p[1]*part2->p[1]-part2->p[2]*part2->p[2]
	-part2->p[3]*part2->p[3];
	p1dotp2=0.0;
	for(alpha=0;alpha<4;alpha++){
		p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
	}
	//Plab=sqrt((p1dotp2*p1dotp2/(part2->msquared))-part1->msquared);
	triangle=p1dotp2*p1dotp2-m1squared*m2squared;
	Plab=0.5*(m1squared+m2squared)*triangle/(m1squared*m2squared);
	Plab=sqrt(Plab);
	sigma_annihilation=6.7*pow(Plab,-0.7)/double(NSAMPLE);
	rstrange=0.5*sqrt(sigma_annihilation);
	rstrange*=pow(ANNIHILATION_SREDUCTION,abs(part1->resinfo->strange))+pow(ANNIHILATION_SREDUCTION,abs(part2->resinfo->strange));
	sigma_annihilation=rstrange*rstrange;
	return sigma_annihilation;
}

bool CMSU_Boltzmann::CancelAnnihilation(CMSUPart *part1,CMSUPart *part2){
	double mupi,muK,muB,muQtot,betaEtot,netK,netpi;
	double betaB,betameson,EB;
	double reduction_factor=1.0;
	FourVector P,UB;
	int btype1,btype2;
	int ix1,iy1,ix2,iy2;
	double uBx1,uBy1,uBx2,uBy2;
	unsigned int iitau;
	bool cancel;
	CMuTInfo *mti1,*mti2;

	//btype1=part1->resinfo->Btype;
	//btype2=part2->resinfo->Btype;

	iitau=lrint(floor(tau/MUTCALC_DELTAU));
	if(iitau>=muTinfo.size()){
		cancel=false;
		return cancel;
	}
	CMuTInfo::GetIxIy(part1->r[1],part1->r[2],ix1,iy1);
	CMuTInfo::GetIxIy(part1->r[1],part1->r[2],ix2,iy2);
	if(iitau>=muTinfo.size() || ix1>=CMuTInfo::NXY || ix2>=CMuTInfo::NXY || iy1>=CMuTInfo::NXY || iy2>=CMuTInfo::NXY || ix1<0 || ix2<0 || iy1<0 || iy2<0){
		cancel=false;
		return cancel;
	}
	
	if(tau<CMuTInfo::taumin[ix1][iy1] || tau<CMuTInfo::taumin[ix2][iy2]){
		//printf("annihilation canceled for being in hydro area\n");
		cancel=true;
		return cancel;
	}

	mti1=muTinfo[iitau][ix1][iy1];
	mti2=muTinfo[iitau][ix2][iy2];
	btype1=part1->resinfo->Btype;
	btype2=part2->resinfo->Btype;
	if(btype1==-1 || btype2==-1 || mti1->Tpi<0.05 || mti2->Tpi<0.05){
		cancel=false;
		//printf("annihilation proceeding, due to small T or funky btype\n");
		//printf("Tpi=%g,%g,  btype=%d,%d\n",mti1->Tpi,mti2->Tpi,btype1,btype2);
		return cancel;
	}
	if(!mti1->sufficientNB[btype1] || !mti2->sufficientNB[btype2] || !mti1->sufficientNpi || !mti2->sufficientNpi || !mti1->sufficientNK || !mti2->sufficientNK){
		cancel=false;
		//printf("annihilation proceeding because of insufficent Npi-NK-NB\n");
		//printf("     Npi=%d,%d, NK=%d,%d, NB=%d,%d\n",mti1->Npi,mti2->Npi,mti1->NK,mti2->NK,
		//	mti1->NB[btype1],mti2->NB[btype2]);
		return cancel;
	}

	netK=fabs(part1->resinfo->strange)+fabs(part2->resinfo->strange);
	if(netK>5)
		netK=4;
	netpi=5.0-netK;
	mupi=0.5*(mti1->mupi+mti2->mupi);
	muK=0.5*(mti1->muK+mti2->muK);
	muB=mti1->muB[btype1]+mti2->muB[btype2];
	muQtot=netpi*mupi+netK*muK-muB;
	for(int alpha=0;alpha<4;alpha++)
		P[alpha]=part1->p[alpha]+part2->p[alpha];

	uBx1=mti1->UxB[btype1];
	if(part1->r[1]<0.0)
		uBx1=-uBx1;
	uBy1=mti1->UyB[btype1];
	if(part1->r[2]<0.0)
		uBy1=-uBy1;
	
	uBx2=mti2->UxB[btype2];
	if(part2->r[1]<0.0)
		uBx2=-uBx2;
	uBy2=mti2->UyB[btype2];
	if(part2->r[2]<0.0)
		uBy2=-uBy2;

	UB[1]=0.5*(uBx1+uBx2);
	UB[2]=0.5*(uBy1+uBy2);
	UB[3]=sinh(0.5*(part1->eta+part2->eta));
	UB[0]=sqrt(1.0+UB[1]*UB[1]+UB[2]*UB[2]+UB[3]*UB[3]);
	EB=DotProduct(UB,P);
	betaB=0.5*(1.0/mti1->TB[btype1]+1.0/mti2->TB[btype2]);
	betameson=(netpi*(1.0/mti1->Tpi+1.0/mti2->Tpi)+netK*(1.0/mti1->TK+1.0/mti2->TK))/10.0;
	betaEtot=EB*(betameson-betaB);
	reduction_factor=1.0-exp(muQtot-betaEtot);
	
	/*	
	if(reduction_factor<0.0){
		printf("________________________________________________________\n");
		printf("-- tau=%g, reduction factor=%g, (%d,%d), netK=%g\n",
			tau,reduction_factor,part1->resinfo->pid,part2->resinfo->pid,netK);
		printf("-- btype1=%d, btype2=%d\n",btype1,btype2);
		printf("-- Tpi=%g, TK=%g\n",0.5*(mti1->Tpi+mti2->Tpi),0.5*(mti1->TK+mti2->TK));
		printf("-- TB1=%g, TB2=%g\n",mti1->TB[btype1],mti2->TB[btype2]);
		printf("-- mupi=%g, muK=%g\n",0.5*(mti1->mupi+mti2->mupi),0.5*(mti1->muK+mti2->muK));
		printf("-- muB1+muB2=%g\n",mti1->muB[btype1]+mti2->muB[btype2]);
		printf("-- muQtot=%g, betaEtot=%g\n",muQtot,betaEtot);
		printf("________________________________________________________\n");
	}
	*/
	
	if(mti1->TB[btype1]!=mti1->TB[btype1] || mti2->TB[btype2]!=mti1->TB[btype2]){
		CLog::Fatal("temperature screwed up\n");
	}
	if(betaB!=betaB || betameson!=betameson){
		CLog::Fatal("betaB or betameson screwed up\n");
	}
	if(reduction_factor!=reduction_factor){
		CLog::Fatal("reduction_factor="+to_string(reduction_factor)+"\n");
	}
		
	if(randy->ran()<reduction_factor){
		//printf("annihilation proceeding, tau=%g, reduction_factor=%g\n",tau,reduction_factor);
		cancel=false;
		return cancel;
	}

	/*
	printf("------------------------------------\n");
	printf("annihilation canceled, tau=%g, reduction_factor=%g, betameson=%g, betaB=%g, EB=%g\n",tau,reduction_factor,betameson,betaB,EB);
	printf("xyz=(%g,%g,%g),    (%g,%g,%g)\n",part1->r[1],part1->r[2],part1->r[3],part2->r[1],part2->r[2],part2->r[3]);
	printf("UB=(%g,%g,%g)\n",UB[1],UB[2],UB[3]);
	printf("P=(%g,%g,%g,%g)\n",P[0],P[1],P[2],P[3]);
	printf("P1=(%g,%g), (%g,%g)\n",part1->p[1],part1->p[2],part2->p[1],part2->p[2]);
	printf("------------------------------------\n");
	*/
	cancel=true;
	return cancel;
}

