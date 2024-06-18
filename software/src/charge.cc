#include "msu_boltzmann/msu_boltzmann.h"
#include "bfduke/bfcharge.h"
#include "msu_sampler/part.h"
#include "msu_sampler/sampler.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"

using namespace std;
using namespace NMSUPratt;

void CMSU_Boltzmann::GenHadronsFromCharges(){
	int maxbid,bid;
	CHBCharge *chargea,*chargeb;
	pair<CHBChargeMap::iterator,CHBChargeMap::iterator> icpair_even,icpair_odd;
	CHBChargeMap::iterator itc0,itc1;
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
}

void CMSU_Boltzmann::GenHadronsFromCharge(int balanceID,CHBCharge *charge){
	int ires;
	double delN,bweight;
	double rapidity,mass;
	Chyper *hyper=&(charge->hyper);
	CMSUPart *part;
	FourVector p;
	Csampler *sampler;
	CHBChargeMap::iterator it;
	Eigen::Vector3d Qprime,Q,q;
	CresInfoMap::iterator itr;
	CresInfo *resinfo;
	if(hyper->T0>mastersampler->TFmin){
		sampler=mastersampler->ChooseSampler(hyper);
		sampler->CalcNHadrons(hyper);
		sampler->CalcChiWithFugacity(hyper);

		Q(0)=charge->q[0];
		Q(1)=charge->q[1];
		Q(2)=charge->q[2];
		//Qprime=sampler->chiinv0*Q;
		Qprime=sampler->chiinv*Q;

		for(itr=reslist->resmap.begin();itr!=reslist->resmap.end();++itr){
			resinfo=itr->second;
			ires=resinfo->ires;
			if(resinfo->baryon!=0 || resinfo->charge!=0 || resinfo->strange!=0){
				if(!balancearrays->PPBAR_ONLY 
				|| (balancearrays->PPBAR_ONLY && resinfo->baryon!=0 && abs(resinfo->pid)!=2112)){
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
						if(abs(resinfo->pid)==211)
							Npions_fromcharges+=1;
						if(abs(resinfo->pid)==2212)
							Nprotons_fromcharges+=1;
						randy->increase_threshold();
					}
				}
			}
		}
	}
}

void CMSU_Boltzmann::ReadCharges(int ichargefile){
	string dirname="udsdata/"+qualifier;
	char chargefile[20];
	snprintf(chargefile,20,"%d",ichargefile);
	string filename=dirname+"/"+"uds"+chargefile+".txt";
	Chyper *hyper;
	int maxbid=0;
	char dummy[120];
	vector<double> etaboost;
	CHBChargeMap::iterator it;
	CHBCharge *charge;
	int balanceID,qu,qd,qs,bidcharge;
	double u0,ux,uy,x,y,tau_read,eta,w,dOmega0,dOmegaX,dOmegaY,pitildexx,T0;
	double fugacity_u,fugacity_d,fugacity_s;
	double pitildexy,pitildeyy;
	Csampler *sampler;
	CLog::Info("opening uds file "+filename+"\n");
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,120,fptr);
	chargemap.clear();
	int iread=0;
	do{
		fscanf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&balanceID,&qu,&qd,&qs,&w,
		&tau_read,&eta,&x,&y,
		&T0,&ux,&uy,&dOmega0,&dOmegaX,&dOmegaY,
		&pitildexx,&pitildeyy,&pitildexy,
		&fugacity_u,&fugacity_d,&fugacity_s);
		fgets(dummy,120,fptr);
		//printf("iread=%d, bid=%d, tau=%g, x=%g, y=%g, pitildexy=%g\n",iread,balanceID,tau_read,x,y,pitildexy);
		iread+=1;
		if(!feof(fptr)){
			charge=new CHBCharge();
			hyper=&(charge->hyper);
			hyper->T0=T0;
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
			hyper->fugacity_u=fugacity_u;
			hyper->fugacity_d=fugacity_d;
			hyper->fugacity_s=fugacity_s;
			sampler=NULL;
			double tfmin=mastersampler->TFmin;
			if(T0>tfmin){
				sampler=Csampler::mastersampler->ChooseSampler(hyper);
				hyper->T0=sampler->Tf;
				hyper->P=sampler->P0;
				hyper->epsilon=sampler->epsilon0;
			}
			hyper->sampler=sampler;
			if(T0>tfmin){
				chargemap.insert(CHBChargePair(balanceID,charge));
				if(balanceID>maxbid)
					maxbid=balanceID;
			}
			else
				delete charge;
		}
	}while(!feof(fptr));
	CLog::Info("read in "+to_string(iread)+" uds charges\n");
	fclose(fptr);
	Misc::Pause();
	
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
	CHBChargeMap::iterator it,itnext;
	CHBCharge *charge;
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
	pair<CHBChargeMap::iterator,CHBChargeMap::iterator> icpair_even,icpair_odd;
	CHBChargeMap::iterator itc;
	int a,b,bid,maxbid;
	Eigen::Vector3d qa;
	Eigen::Vector3d qb;
	itc=chargemap.end(); itc--;
	maxbid=itc->first;
	CHBCharge *chargea,*chargeb;
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
