#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_hydrobalance/charge.h"
#include "msu_sampler/part.h"
#include "msu_sampler/sampler.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"

using namespace std;

void CMSU_Boltzmann::GenHadronsFromCharges(){
	int maxbid,bid;
	CCharge *chargea,*chargeb;
	pair<CChargeMap::iterator,CChargeMap::iterator> icpair_even,icpair_odd;
	CChargeMap::iterator itc0,itc1;
	CMSUPartMap::iterator itp;
	itc1=chargemap.end(); itc1--;
	maxbid=itc1->first;

	for(bid=0;bid<maxbid;bid+=2){
		icpair_even=chargemap.equal_range(bid);
		itc0=icpair_even.first;
		if(itc0!=chargemap.end()){
			chargea=itc0->second;
			icpair_odd=chargemap.equal_range(bid+1);
			itc1=icpair_odd.first;
			if(itc1!=chargemap.end()){
				chargeb=itc1->second;
				GenHadronsFromCharge(bid,chargea);
				GenHadronsFromCharge(bid+1,chargeb);
			}
		}
	}
	/*
	vector<int> nprotons;
	nprotons.resize(maxbid+1);
	for(bid=0;bid<maxbid;bid+=1)
		nprotons[bid]=0;

	CMSUPartMap::iterator ppos;
	CMSUPart *part;
	CresInfo *resinfo;
	int np=0,npi=0,nhad=0;
	for(ppos=PartMap.begin();ppos!=PartMap.end();++ppos){
		part=ppos->second;
		resinfo=part->resinfo;
		if(part->balanceID>=0){
			nhad+=1;
			if(abs(resinfo->pid)==2212)
				np+=1;
			if(abs(resinfo->pid)==211)
				npi+=1;
		}
		if(abs(resinfo->pid)==2212 && part->balanceID>=0)
			nprotons[part->balanceID]+=1;
	}
	printf("From BF map, nhad=%d, nprotons=%d, npions=%d, maxbid=%d\n",nhad,np,npi,maxbid);
	for(bid=0;bid<maxbid;bid+=2){
		if(nprotons[bid]*nprotons[bid+1]!=0)
			printf("%d %d\n",nprotons[bid],nprotons[bid+1]);
	}*/
}

void CMSU_Boltzmann::GenHadronsFromCharge(int balanceID,CCharge *charge){
	int ires;
	double delN,bweight;
	double rapidity,mass;
	Chyper *hyper=&(charge->hyper);
	CMSUPart *part;
	FourVector p;
	Csampler *sampler;
	CChargeMap::iterator it;
	Eigen::Vector3d Qprime,Q,q;
	CresInfoMap::iterator itr;
	CresInfo *resinfo;
	sampler=charge->hyper.sampler;
	
	Q(0)=charge->q[0];
	Q(1)=charge->q[1];
	Q(2)=charge->q[2];
	Qprime=sampler->chiinv0*Q;

	for(itr=reslist->resmap.begin();itr!=reslist->resmap.end();++itr){
		resinfo=itr->second;
		ires=resinfo->ires;
		if(resinfo->baryon!=0 || resinfo->charge!=0 || resinfo->strange!=0){ 
			q[0]=resinfo->q[0]; q[1]=resinfo->q[1]; q[2]=resinfo->q[2];
			delN=sampler->density0i[ires]*(q.dot(Qprime)); // number of hadrons to create
			bweight=charge->weight*delN/fabs(delN);
			randy->increment_netprob(fabs(delN*NSAMPLE_UDS2BAL));
			while(randy->test_threshold(0.0)){
				sampler->GetP(&(charge->hyper),hyper->T0,resinfo,p);
				mass=resinfo->mass;
				part=GetDeadPart();
				rapidity=charge->eta+asinh(p[3]/sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]));
				part->InitBalance(resinfo->pid,charge->x,charge->y,charge->tau,charge->eta,p[1],p[2],mass,rapidity,bweight,balanceID);
				randy->increase_threshold();
			}
		}
	}
}


void CMSU_Boltzmann::ReadCharges(int ichargefile){
	string dirname="udsdata/"+qualifier;
	char chargefile[10];
	sprintf(chargefile,"%d",ichargefile);
	//string filename=dirname+"/"+parmap.getS("CHARGESINFO_FILENAME","uds.dat");
	string filename=dirname+"/"+"uds"+chargefile+".dat";
	Chyper *hyper;
	int maxbid=0;
	char dummy[120];
	vector<double> etaboost;
	CChargeMap::iterator it;
	CCharge *charge;
	int balanceID,qu,qd,qs,bidcharge;
	double u0,ux,uy,x,y,tau_read,eta,w,dOmega0,dOmegaX,dOmegaY,pitildexx;
	double pitildexy,pitildeyy;
	Csampler *sampler;
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,120,fptr);
	chargemap.clear();
	do{
		fscanf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&balanceID,&qu,&qd,&qs,&w,&tau_read,&eta,&x,&y,&ux,&uy,&dOmega0,&dOmegaX,&dOmegaY,&pitildexx,&pitildeyy,&pitildexy);
		if(!feof(fptr)){
			charge=new CCharge();
			hyper=&(charge->hyper);
			hyper->T0=parmap->getD("FREEZEOUT_TEMP",0.155);
			charge->q[0]=qu;
			charge->q[1]=qd;
			charge->q[2]=qs;
			charge->weight=w;
			charge->tau=tau_read;
			charge->eta=eta;
			charge->x=x;
			charge->y=y;
			hyper->tau=tau_read;
			hyper->r[1]=x;
			hyper->r[2]=y;
			hyper->dOmega[0]=dOmega0;
			hyper->dOmega[1]=dOmegaX;
			hyper->dOmega[2]=dOmegaY;
			u0=sqrt(1.0+ux*ux+uy*uy);
			hyper->u[1]=ux;
			hyper->u[2]=uy;
			hyper->udotdOmega=u0*dOmega0-ux*dOmegaX-uy*dOmegaY;
			hyper->pitilde[1][1]=pitildexx;
			hyper->pitilde[2][2]=pitildeyy;
			hyper->pitilde[1][2]=hyper->pitilde[2][1]=pitildexy;
			sampler=Csampler::mastersampler->ChooseSampler(hyper);
			hyper->sampler=sampler;
			hyper->T0=sampler->Tf;
			hyper->P=sampler->P0;
			hyper->epsilon=sampler->epsilon0;
			chargemap.insert(CChargePair(balanceID,charge));
			if(balanceID>maxbid)
				maxbid=balanceID;
		}
	}while(!feof(fptr));
	fclose(fptr);
	
	etaboost.resize((maxbid+1)/2);
	for(bidcharge=0;bidcharge<(maxbid+1)/2;bidcharge+=1){
		etaboost[bidcharge]=ETAMAX*(1.0-2.0*randy->ran());
	}
	for(it=chargemap.begin();it!=chargemap.end();++it){
		charge=it->second;
		balanceID=it->first;
		bidcharge=floorl(balanceID/2);
		charge->eta+=etaboost[bidcharge];
		while(charge->eta>ETAMAX){
			charge->eta-=2.0*ETAMAX;
		}
		while(charge->eta<-ETAMAX){
			charge->eta+=2.0*ETAMAX;
		}
		hyper=&(charge->hyper);
	}
	etaboost.clear();
	//CalcChiTotFromQ();
}

void CMSU_Boltzmann::DeleteCharges(){
	CChargeMap::iterator it,itnext;
	CCharge *charge;
	it=chargemap.begin();
	while(it!=chargemap.end()){
		itnext=it;
		++itnext;
		charge=it->second;
		if(charge!=NULL)
			delete charge;
		it=itnext;
	}
	chargemap.clear();
}

void CMSU_Boltzmann::IncrementChiTotFromCharges(){
	pair<CChargeMap::iterator,CChargeMap::iterator> icpair_even,icpair_odd;
	CChargeMap::iterator itc;
	int a,b,bid,maxbid;
	Eigen::Vector3d qa;
	Eigen::Vector3d qb;
	itc=chargemap.end(); itc--;
	maxbid=itc->first;
	CCharge *chargea,*chargeb;
	for(bid=0;bid<maxbid;bid+=2){
		icpair_even=chargemap.equal_range(bid);
		itc=icpair_even.first;
		chargea=itc->second;
		icpair_odd=chargemap.equal_range(bid+1);
		itc=icpair_odd.first;
		chargeb=itc->second;
		for(a=0;a<3;a++){
			for(b=0;b<3;b++){
				chitotQ(a,b)-=0.5*(chargea->q[a]*chargeb->q[b]+chargea->q[b]*chargeb->q[a]);
			}
		}
	}
}

void CMSU_Boltzmann::IncrementChiTotFromHadrons(){
	CMSUPartMap bfpartmap;
	int maxbid=-1;
	pair<CMSUPartMap::iterator,CMSUPartMap::iterator> itpair_even,itpair_odd;
	CMSUPartMap::iterator it,ita0,itaf,itb0,itbf,ita,itb;
	CMSUPart *parta,*partb,*part;
	for(it=PartMap.begin();it!=PartMap.end();++it){
		part=it->second;
		if(part->balanceID>maxbid)
			maxbid=part->balanceID;
		if(part->balanceID>=0)
			bfpartmap.insert(CMSUPartPair(part->balanceID,part));
	}
	
	int bid;
	double dchi;
	
	for(bid=0;bid<=maxbid;bid+=2){
		itpair_even=bfpartmap.equal_range(bid);
		ita0=itpair_even.first;
		itaf=itpair_even.second;
		itpair_odd=bfpartmap.equal_range(bid+1);
		itb0=itpair_odd.first;
		itbf=itpair_odd.second;
		if(ita0!=itaf && itb0!=itbf){
			for(ita=ita0;ita!=itaf;++ita){
				parta=ita->second;
				for(itb=itb0;itb!=itbf;++itb){
					partb=itb->second;
					for(int a=0;a<3;a++){
						for(int b=0;b<3;b++){
							dchi=0.5*parta->resinfo->q[a]*partb->resinfo->q[b]*parta->bweight*partb->bweight
								+0.5*parta->resinfo->q[b]*partb->resinfo->q[a]*parta->bweight*partb->bweight;
							chitotH(a,b)-=0.5*dchi;
							chitotH(b,a)-=0.5*dchi;
						}
					}
				}
			}
		}
	}
	bfpartmap.clear();
}
