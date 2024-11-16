#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/msupart.h"
#include "msu_eos/resonances.h"
#include "msu_boltzmann/cell.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/constants.h"

using namespace std;
using namespace NMSUPratt;

//#define __SCATTERING_ON__

// If two particles pass one another, Collide will determine whether to scatter and how

int CMSU_Boltzmann::Collide_Scatter(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product){
	bool bjtranslate=false;
	int colltype;

	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}
	CMSUPart *part3=nullptr,*part4=nullptr;
	part3=product[0];
	part4=product[1];
	part3->Copy(part1);
	part4->Copy(part2);
	nproducts=2;

	Scatter(part1,part2,part3,part4);
	part3->balanceID=part1->balanceID;
	part4->balanceID=part2->balanceID;

	if(bjtranslate)
		part3->BjorkenUnTranslate();
	if(part1->balanceID>=0 || part2->balanceID>=0){
		if(bjtranslate)
			part1->BjorkenUnTranslate();
		colltype=-2;
	}
	else{
		colltype=2;
		nscatter+=1;
	}
	return colltype;
}

int CMSU_Boltzmann::Collide_Merge(CMSUPart *part1,CMSUPart *part2,double sigma_merge,vector<double> &dsigma_merge,int &nproducts,array<CMSUPart*,5> &product){
	bool bjtranslate=false,msuccess;
	int ir1,ir2,irflip,colltype,imerge;
	double r,sigmatot;
	Cmerge *merge;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}
	CMSUPart *part3=product[0];
	part3->CopyPositionInfo(part2);
	ir1=part1->resinfo->ires;
	ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}

	merge=reslist->MergeArray[ir1][ir2];
	r=randy->ran();
	imerge=0;
	sigmatot=dsigma_merge[0];
	while(r>sigmatot/sigma_merge && merge!=nullptr){
		if(merge==nullptr || sigmatot>sigma_merge){
			snprintf(message,CLog::CHARLENGTH,"In CMSU_Boltzmann::Collide_Merge, merge is nullptr?? or sigmatot/sigma_merge=%g is >1, sigma_merge=%g, sigmatot=%g\n",sigmatot/sigma_merge,sigma_merge,sigmatot);
			part1->Print();
			part2->Print();
			CLog::Fatal(message);
		}
		imerge+=1;
		merge=merge->next;
		sigmatot+=dsigma_merge[imerge];
	}
	if(merge==nullptr){
		nproducts=0;
		if(bjtranslate){
			part1->BjorkenUnTranslate();
		}
		colltype=0;
	}
	else{
		msuccess=Merge(part1,part2,part3,merge->resinfo);
		if(msuccess){
			nproducts=1;
			colltype=1;
			nmerge+=1;
		}
		else{
			if(bjtranslate){
				part1->BjorkenUnTranslate();
			}
			nproducts=0;
			colltype=0;
		}
	}
	return colltype;
}

int CMSU_Boltzmann::Collide_Annihilate(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product){
	bool bjtranslate=false;
	int colltype;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}

	if(CancelAnnihilation(part1,part2)){
		if(bjtranslate)
			part1->BjorkenUnTranslate();
		nregenerate+=1;
		nproducts=0;
		colltype=0;
		ncancel_annihilate+=1;
	}
	else{
		Annihilate(part1,part2,nproducts,product);
		if(bjtranslate && fabs(product[0]->eta)>ETAMAX){
			for(int iprod=0;iprod<nproducts;iprod++)
				product[iprod]->CyclicReset();
		}
		nannihilate+=1;
		colltype=4;
	}
	return colltype;
}

// warning Inelastic Scattering NEEDS to be tested
int CMSU_Boltzmann::Collide_Inelastic(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product){
	snprintf(message,CLog::CHARLENGTH,"INELASTIC IS UNTESTED!!!!\n");
	CLog::Fatal(message);
	const int NWMAX=5000;
	const double g[4]={1,-1,-1,-1};
	CMSUPart *part3=product[0];
	CMSUPart *part4=product[1];
	nproducts=2;
	int colltype;
	double inel_weight[NWMAX]={0.0};
	//double inel_d=0.0;
	double q_prime,degen1_i,degen2_i,M,P2,p1dotp2,r,wtot;
	list<CInelasticInfo>::iterator inel;
	list<CInelasticInfo> inel_list;
	int iw,ir1,ir2,irflip,alpha,G_Value,netq,netb,nets;
	bool bjtranslate=false,G_Parity=false,success;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}

	ir1=part1->resinfo->ires; ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}

	if(part1->resinfo->G_Parity && part2->resinfo->G_Parity){
		G_Parity = true;
		G_Value = part1->resinfo->G_Parity * part2->resinfo->G_Parity;
	}
	p1dotp2=0.0;
	for(alpha=0;alpha<4;alpha++){
		p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
	}
	P2=part1->msquared+part2->msquared+2.0*p1dotp2;
	M=sqrt(P2);
			// First calculate denominator

	if(inelasticlist->UseInelasticArray){
		inel_list = inelasticlist->InelasticArray[ir1][ir2];
	}else{
		netq = part1->resinfo->charge+part2->resinfo->charge;
		netb = part1->resinfo->baryon+part2->resinfo->baryon;
		nets = part1->resinfo->strange+part2->resinfo->strange;
		inel_list = inelasticlist->ThermalArray[abs(netb)][abs(netq)][abs(nets)][Misc::Sign(netb)][Misc::Sign(netq)][Misc::Sign(nets)];
	}
	//inel_d = Q0*Q0;
	inel = inel_list.begin();
	iw=0;
	r=randy->ran();
	success=false;
	//inel_d=0.0;
	while(inel!=inel_list.end()){
		if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || (!G_Parity)){
			if(inel->resinfo_1->mass+inel->resinfo_2->mass<M){
				degen1_i=inel->resinfo_1->degen;
				degen2_i=inel->resinfo_2->degen;
				q_prime = Misc::triangle(M,inel->resinfo_1->mass,inel->resinfo_2->mass);
				inel_weight[iw] = degen1_i*degen2_i*q_prime;
			}
			else{
				inel_weight[iw]=0.0;
			}
			//inel_d+=inel_weight[iw];
		}
	}
	wtot=0.0;
	iw=0;
	inel=inel_list.begin();
	success=false;
	do{
		wtot+=inel_weight[iw];
		if(wtot>r){
			success=true;
		}
		else{
			iw+=1;
			++inel;
		}
	}while(!success && inel!=inel_list.end());
	if(success){
		netb=part1->resinfo->baryon+part2->resinfo->baryon;
		netq=part1->resinfo->charge+part2->resinfo->charge;
		nets=part1->resinfo->strange+part2->resinfo->strange;
		InelasticScatter(part1,part2,part3,part4,*inel);
		if(bjtranslate)
			part3->BjorkenUnTranslate();
		netb-=part3->resinfo->baryon+part4->resinfo->baryon;
		netq-=part3->resinfo->charge+part4->resinfo->charge;
		nets-=part3->resinfo->strange+part4->resinfo->strange;
		if(netb!=0 || netq!=0 || nets!=0){
			snprintf(message,CLog::CHARLENGTH,"WARNING: charge not conserved in inel collision, netb=%d, netq=%d, nets=%d\n",netb,netq,nets);
			CLog::Info(message);
		}
		nproducts=2;
		ninelastic+=1;
		colltype=3;
	}
	else{
		nproducts=0;
		colltype=0;
	}

	ninelastic+=1;
	return colltype;
}


