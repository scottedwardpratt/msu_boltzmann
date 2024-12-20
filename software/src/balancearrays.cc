#include "msu_boltzmann/balancearrays.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/acceptance.h"
#include "msu_commonutils/randy.h"
#include "msu_eos/resonances.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

CBalanceArrays::CBalanceArrays(CMSU_Boltzmann *boltzmannset){
	boltzmann=boltzmannset;
	CBFNumer::message=message;
	CBFDenom::message=message;
	reslist=boltzmann->reslist;
	parmap=&(boltzmann->parmap);
	InitArrays();
}

void CBalanceArrays::InitArrays(){
	CresInfoMap::iterator itr;
	NSAMPLE_HYDRO2UDS=parmap->getD("HYDRO2UDS_NSAMPLE",1);
	NSAMPLE_UDS2BAL=parmap->getD("MSU_BOLTZMANN_NSAMPLE_UDS2BAL",1);
	FROM_UDS=parmap->getB("BF_FROM_UDS",false);
	BF_YMAX=parmap->getD("BF_YMAX",1.0);
	BF_YMIN=-BF_YMAX;
	BF_PHICUT=parmap->getD("BF_PHICUT",15);
	NPHI=parmap->getI("BF_NPHIBINS",180);
	DELY=parmap->getD("BF_DELY",0.05);
	PPBAR_ONLY=parmap->getB("BF_PPBAR_ONLY",false);
	gammap=gammapnorm=v2=normtest=v2norm,v2prime=v2primenorm=v2perfect=v2perfectnorm=0.0;
	Nchi=parmap->getI("NCHI",200);
	CreateBFArrays();
	acceptance_description=parmap->getS("BF_ACCEPTANCE","PERFECT");
	NoKsNoPhi=parmap->getB("BF_NoKsNoPhi",false);
	if(acceptance_description=="CHEAP"){
		acceptance=new CAcceptance_CHEAP(parmap);
	}
	else if(acceptance_description=="STAR"){
		acceptance=new CAcceptance_STAR(parmap);
	}
	else if(acceptance_description=="ALICE"){
		acceptance=new CAcceptance_ALICE(parmap);
	}
	else if(acceptance_description=="ALICE_PERFECT"){
		acceptance=new CAcceptance_ALICE_Perfect(parmap);
	}
	else{
		snprintf(message,CLog::CHARLENGTH,"Define BF_ACCEPTANCE in mod_parameters.txt\n");
		CLog::Fatal(message);
	}
	CBFNumer::acceptance=acceptance;
	NSAMPLE_HYDRO2UDS=NSAMPLE_HYDRO2UDS/(2.0*boltzmann->ETAMAX);
	NEVENTS=0;
	GAMMAP_OS=GAMMAP_SS=0.0;
	NORM_GAMMAP_OS=NORM_GAMMAP_SS=0;
}

void CBalanceArrays::Reset(){
	numer_pipi->Reset();
	numer_piK->Reset();
	numer_pip->Reset();
	numer_KK->Reset();
	numer_Kp->Reset();
	numer_pp->Reset();
	numer_allcharges->Reset();
	denom_pi->Reset();
	denom_K->Reset();
	denom_p->Reset();
	denom_allcharges->Reset();
	bf_pipi->Reset();
	bf_piK->Reset();
	bf_pip->Reset();
	bf_KK->Reset();
	bf_Kp->Reset();
	bf_pp->Reset();
	bf_allcharges->Reset();
	numer_allcharges_phi0->Reset();
	numer_allcharges_phi45->Reset();
	numer_allcharges_phi90->Reset();
	bf_allcharges_phi0->Reset();
	bf_allcharges_phi45->Reset();
	bf_allcharges_phi90->Reset();
	denom_allcharges_phi0->Reset();
	denom_allcharges_phi45->Reset();
	denom_allcharges_phi90->Reset();
}

void CBalanceArrays::IncrementGammaP(CMSUPart *parta,CMSUPart *partb,double effa,double effb){
	double phia,phib;
	double QaQb;
	QaQb=(parta->resinfo->charge*parta->bweight)*(partb->resinfo->charge*partb->bweight);
	phia=atan2(parta->p[2],parta->p[1]);
	phib=atan2(partb->p[2],partb->p[1]);
	gammap-=cos(phia+phib)*QaQb*effa*effb;
	gammapnorm+=effa*effb;
	normtest-=effa*effb*QaQb;
}

void CBalanceArrays::PrintBFNumer(){
	
}

void CBalanceArrays::PrintBFDenom(){
	snprintf(message,CLog::CHARLENGTH,"pi+,pi-: nplus=%g, minus=%g\n",denom_pi->Nplus,denom_pi->Nminus);
	snprintf(message,CLog::CHARLENGTH,"%sK+,K-  : nplus=%g, minus=%g\n",message,denom_K->Nplus,denom_K->Nminus);
	snprintf(message,CLog::CHARLENGTH,"%sp,pbar : nplus=%g, minus=%g\n",message,denom_p->Nplus,denom_p->Nminus);
	snprintf(message,CLog::CHARLENGTH,"%sN+,N-  : nplus=%g, minus=%g\n",message,denom_allcharges->Nplus,denom_allcharges->Nminus);
	CLog::Info(message);
}

void CBalanceArrays::CreateBFArrays(){
	numer_pipi=new CBFNumer(parmap);
	numer_piK=new CBFNumer(parmap);
	numer_pip=new CBFNumer(parmap);
	numer_KK=new CBFNumer(parmap);
	numer_Kp=new CBFNumer(parmap);
	numer_pp=new CBFNumer(parmap);
	numer_allcharges=new CBFNumer(parmap);
	numer_allcharges_phi0=new CBFNumer(parmap);
	numer_allcharges_phi45=new CBFNumer(parmap);
	numer_allcharges_phi90=new CBFNumer(parmap);
	
	numer_pipi->name="pipi";
	numer_piK->name="piK";
	numer_pip->name="pip";
	numer_KK->name="KK";
	numer_Kp->name="Kp";
	numer_pp->name="pp";
	numer_allcharges->name="allcharges";
	numer_allcharges_phi0->name="allcharges_phi0";
	numer_allcharges_phi45->name="allcharges_phi45";
	numer_allcharges_phi90->name="allcharges_phi90";
	
	bf_pipi=new CBFNumer(parmap);
	bf_piK=new CBFNumer(parmap);
	bf_pip=new CBFNumer(parmap);
	bf_KK=new CBFNumer(parmap);
	bf_Kp=new CBFNumer(parmap);
	bf_pp=new CBFNumer(parmap);
	bf_allcharges=new CBFNumer(parmap);
	bf_allcharges_phi0=new CBFNumer(parmap);
	bf_allcharges_phi45=new CBFNumer(parmap);
	bf_allcharges_phi90=new CBFNumer(parmap);
	
	bf_pipi->name="pipi";
	bf_piK->name="piK";
	bf_pip->name="pip";
	bf_KK->name="KK";
	bf_Kp->name="Kp";
	bf_pp->name="pp";
	bf_allcharges->name="allcharges";
	bf_allcharges_phi0->name="allcharges_phi0";
	bf_allcharges_phi45->name="allcharges_phi45";
	bf_allcharges_phi90->name="allcharges_phi90";
	
	denom_pi=new CBFDenom(parmap);
	denom_K=new CBFDenom(parmap);
	denom_p=new CBFDenom(parmap);
	denom_allcharges=new CBFDenom(parmap);
	denom_allcharges_phi0=new CBFDenom(parmap);
	denom_allcharges_phi45=new CBFDenom(parmap);
	denom_allcharges_phi90=new CBFDenom(parmap);
}

void CBalanceArrays::ConstructBFs(){
	bool NoQ=false;
	if(!PPBAR_ONLY){
		ConstructBF(numer_pipi,denom_pi,bf_pipi,2.0,NoQ);
		ConstructBF(numer_piK,denom_pi,bf_piK,1.0,NoQ);
		ConstructBF(numer_pip,denom_pi,bf_pip,1.0,NoQ);
		ConstructBF(numer_KK,denom_K,bf_KK,2.0,NoQ);
		ConstructBF(numer_Kp,denom_K,bf_Kp,1.0,NoQ);
	}
	ConstructBF(numer_pp,denom_p,bf_pp,2.0,NoQ);
	NoQ=true;
	if(!PPBAR_ONLY){
		ConstructBF(numer_allcharges,denom_allcharges,bf_allcharges,2.0,NoQ);
		ConstructBF(numer_allcharges_phi0,denom_allcharges_phi0,bf_allcharges_phi0,1.0,NoQ);
		ConstructBF(numer_allcharges_phi45,denom_allcharges_phi45,bf_allcharges_phi45,1.0,NoQ);
		ConstructBF(numer_allcharges_phi90,denom_allcharges_phi90,bf_allcharges_phi90,1.0,NoQ);
	}
	snprintf(message,CLog::CHARLENGTH,"v2=%g, v2prime=%g\n",v2/v2norm,v2prime/v2primenorm);
	CLog::Info(message);
}

void CBalanceArrays::ConstructBF(CBFNumer *numer,CBFDenom *denom,CBFNumer *bf,double doublecount,bool NoQ){
	
	bf->Bqinv=numer->Bqinv;
	bf->Bqout=numer->Bqout;
	bf->Bqside=numer->Bqside;
	bf->Bqlong=numer->Bqlong;
	bf->Beta=numer->Beta;
	bf->Beta1=numer->Beta1;
	bf->By=numer->By;
	bf->By1=numer->By1;
	bf->Bphi=numer->Bphi;
	
	bf->npairs=numer->npairs;
	int ibin;
	double norm;
	double N=denom->Nplus+denom->Nminus;
	if(FROM_UDS){
		N*=2.0*NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL*NSAMPLE_UDS2BAL;
	}
	if(!NoQ){
		for(ibin=0;ibin<numer->Nqbins;ibin++){
			bf->Bqinv[ibin]=doublecount*numer->Bqinv[ibin]/(N*numer->Dqinv);
			bf->Bqout[ibin]=doublecount*numer->Bqout[ibin]/(N*numer->Dqinv);
			bf->Bqside[ibin]=doublecount*numer->Bqside[ibin]/(N*numer->Dqinv);
			bf->Bqlong[ibin]=doublecount*numer->Bqlong[ibin]/(N*numer->Dqinv);
			bf->Cqinv[ibin]=doublecount*numer->Cqinv[ibin]/(N*numer->Dqinv);
			bf->Cqout[ibin]=doublecount*numer->Cqout[ibin]/(N*numer->Dqinv);
			bf->Cqside[ibin]=doublecount*numer->Cqside[ibin]/(N*numer->Dqinv);
			bf->Cqlong[ibin]=doublecount*numer->Cqlong[ibin]/(N*numer->Dqinv);
		}
	}
	norm=0.0;
	for(ibin=0;ibin<numer->Netabins;ibin++){
		bf->Beta[ibin]=doublecount*numer->Beta[ibin]/(N*numer->Deta);
		bf->Beta1[ibin]=doublecount*numer->Beta1[ibin]/(N*numer->Deta);
		bf->Betas[ibin]=doublecount*numer->Betas[ibin]/(N*numer->Deta);
		bf->Ceta[ibin]=doublecount*numer->Ceta[ibin]/(N*numer->Deta);
		bf->Ceta1[ibin]=doublecount*numer->Ceta1[ibin]/(N*numer->Deta);
		bf->Cetas[ibin]=doublecount*numer->Cetas[ibin]/(N*numer->Deta);
	}
	for(ibin=0;ibin<numer->Nybins;ibin++){
		bf->By[ibin]=doublecount*numer->By[ibin]/(N*numer->Dy);
		bf->Cy[ibin]=doublecount*numer->Cy[ibin]/(N*numer->Dy);
		bf->By1[ibin]=doublecount*numer->By1[ibin]/(N*numer->Dy);
		bf->Cy1[ibin]=doublecount*numer->Cy1[ibin]/(N*numer->Dy);
	}
	for(ibin=0;ibin<numer->Nphibins;ibin++){
		bf->Bphi[ibin]=doublecount*numer->Bphi[ibin]/(N*numer->Dphi);
		bf->Cphi[ibin]=doublecount*numer->Cphi[ibin]/(N*numer->Dphi);
		norm+=numer->Dphi*bf->Bphi[ibin];
	}
	int iy,iphi;
	for(iy=0;iy<numer->Netabins;iy++){
		for(iphi=0;iphi<numer->Nphibins;iphi++){
			bf->Byphi[iy][iphi]=numer->Byphi[iy][iphi];
			bf->Cyphi[iy][iphi]=numer->Cyphi[iy][iphi];
		}
	}
	snprintf(message,CLog::CHARLENGTH,"%7s: normalization=%g, npairs=%lld\n",bf->name.c_str(),norm,bf->npairs);
	CLog::Info(message);
}

void CBalanceArrays::WriteNumers(){
	string numertype="numer";
	bool NoQ=false;
	if(!PPBAR_ONLY){
		numer_pipi->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_piK->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_pip->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_KK->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_Kp->WriteNumer(bf_results_dirname,numertype,NoQ);
	}
	numer_pp->WriteNumer(bf_results_dirname,numertype,NoQ);
	NoQ=true;
	if(!PPBAR_ONLY){
		numer_allcharges->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_allcharges_phi0->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_allcharges_phi45->WriteNumer(bf_results_dirname,numertype,NoQ);
		numer_allcharges_phi90->WriteNumer(bf_results_dirname,numertype,NoQ);
	}
}

void CBalanceArrays::WriteBFs(){
	string numertype="bf";
	bool NoQ=false;
	if(!PPBAR_ONLY){
		bf_pipi->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_piK->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_pip->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_KK->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_Kp->WriteNumer(bf_results_dirname,numertype,NoQ);
	}
	bf_pp->WriteNumer(bf_results_dirname,numertype,NoQ);
	NoQ=true;
	if(!PPBAR_ONLY){
		bf_allcharges->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_allcharges_phi0->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_allcharges_phi45->WriteNumer(bf_results_dirname,numertype,NoQ);
		bf_allcharges_phi90->WriteNumer(bf_results_dirname,numertype,NoQ);
	}
}

void CBalanceArrays::WriteGammaP(){
	double Ntot,mult,gp;
	string filename;
	filename=bf_results_dirname+"/gammap.txt";
	FILE *fptr=fopen(filename.c_str(),"w");
	if(FROM_UDS){
		Ntot=denom_allcharges->Nplus+denom_allcharges->Nminus;
		mult=Ntot*BF_YMAX/(boltzmann->ETAMAX*NEVENTS);
		gp=gammap/(2.0*Ntot*NSAMPLE_UDS2BAL*NSAMPLE_UDS2BAL*NSAMPLE_HYDRO2UDS);
		gp=gp/(0.5*mult);
		normtest=0.25*normtest/(0.5*Ntot*NSAMPLE_UDS2BAL*NSAMPLE_UDS2BAL*NSAMPLE_HYDRO2UDS);
		//normtest=normtest/(0.5*mult);
	}
	else{
		Ntot=denom_allcharges->Nplus+denom_allcharges->Nminus;
		mult=Ntot/(NEVENTS*boltzmann->ETAMAX/BF_YMAX);
		gp=gammap/(Ntot*0.5*mult);
	}
	fprintf(fptr,"%g %g %g %g\n",gp,2.0*gammap/gammapnorm,normtest,mult);
	snprintf(message,CLog::CHARLENGTH,"gammap=%g =? %g, gammapnorm=%g, normtest=%g, mult=%g, NSAMPLE_HYDRO2UDS=%g, NSAMPLE_UDS2BAL=%d, NEVENTS=%d, Ntot=%g\n",
	gp,2.0*gammap/gammapnorm,gammapnorm,normtest,mult,NSAMPLE_HYDRO2UDS,NSAMPLE_UDS2BAL,NEVENTS,Ntot);
	CLog::Info(message);
	fclose(fptr);
}

void CBalanceArrays::WriteDenoms(){
	FILE *fptr;
	string filename=bf_results_dirname+"/denom.txt";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"denom_pi: %g %g\n",denom_pi->Nplus,denom_pi->Nminus);
	fprintf(fptr,"denom_K: %g %g\n",denom_K->Nplus,denom_K->Nminus);
	fprintf(fptr,"denom_p: %g %g\n",denom_p->Nplus,denom_p->Nminus);
	fprintf(fptr,"denom_allcharges: %g %g\n",denom_allcharges->Nplus,denom_allcharges->Nminus);
	fprintf(fptr,"denom_allcharges_phi0: %g %g\n",denom_allcharges_phi0->Nplus,denom_allcharges_phi0->Nminus);
	fprintf(fptr,"denom_allcharges_phi45: %g %g\n",denom_allcharges_phi45->Nplus,denom_allcharges_phi45->Nminus);
	fprintf(fptr,"denom_allcharges_phi90: %g %g\n",denom_allcharges_phi90->Nplus,denom_allcharges_phi90->Nminus);
	fprintf(fptr,"dNdy_pi: %g\n",denom_pi->dNdy);
	fprintf(fptr,"dNdy_K: %g\n",denom_K->dNdy);
	fprintf(fptr,"dNdy_p: %g\n",denom_p->dNdy);
	fclose(fptr);
}

void CBalanceArrays::ProcessBFPartMap(){
	int balanceID,pida,pidb,maxbid=-1;
	double delymax=1.6,dely,ya,yb;
	pair<CMSUPartMap::iterator,CMSUPartMap::iterator> itpair_even,itpair_odd;
	CMSUPartMap::iterator it,ita0,itaf,itb0,itbf,ita,itb;
	CMSUPart *parta,*partb,*part;
	for(it=boltzmann->PartMap.begin();it!=boltzmann->PartMap.end();++it){
		part=it->second;
		if(part->balanceID>maxbid)
			maxbid=part->balanceID;
		if(part->balanceID>=0)
			bfpartmap.insert(CMSUPartPair(part->balanceID,part));
	}
	snprintf(message,CLog::CHARLENGTH,"maxbid=%d, bfpartmap.size=%d\n",maxbid,int(bfpartmap.size()));
	CLog::Info(message);
	
	for(balanceID=0;balanceID<maxbid;balanceID+=2){
		itpair_even=bfpartmap.equal_range(balanceID);
		itpair_odd=bfpartmap.equal_range(balanceID+1);
		ita0=itpair_even.first;
		itaf=itpair_even.second;
		itb0=itpair_odd.first;
		itbf=itpair_odd.second;
		if(ita0!=itaf && itb0!=itbf){
			for(ita=ita0;ita!=itaf;++ita){
				parta=ita->second;
				pida=abs(parta->resinfo->pid);
				if(pida==2212 || pida==211 || pida==321){
					ya=parta->y;
					for(itb=itb0;itb!=itbf;++itb){
						partb=itb->second;
						pidb=abs(partb->resinfo->pid);
						if(pidb==2212 || pidb==211 || pidb==321){
							yb=partb->y;
							dely=fabs(ya-yb);
							delymax=acceptance->GetDelYMax(pida,pidb);
							if(dely<delymax){
								IncrementNumer(parta,partb);
							}
						}
					}
				}
			}
		}
	}
	bfpartmap.clear();
}

void CBalanceArrays::ProcessPartMap(){
  // makes denom + correlations from cascade
	multimap<double,CMSUPart *> ppartmap;
	double ya,yb,dely,delymax,MSU_BOLTZMANN_ETAMAX=boltzmann->ETAMAX;
	CMSUPartMap::iterator it;
	multimap<double,CMSUPart *>::iterator ita,itb;
	int pida,pidb;
	CMSUPart *parta,*partb;
	pair<CMSUPartMap::iterator,CMSUPartMap::iterator> itpair;
	snprintf(message,CLog::CHARLENGTH,"processing %d parts in PartMap\n",int(boltzmann->PartMap.size()));
	CLog::Info(message);
	NEVENTS+=1;
	if(FROM_UDS){
		for(it=boltzmann->PartMap.begin();it!=boltzmann->PartMap.end();++it){
			parta=it->second;
			if(abs(parta->resinfo->charge)==1 && parta->balanceID<0){
				if(!PPBAR_ONLY || abs(parta->resinfo->pid)==2212){
					IncrementDenom(parta);
				}
			}
		}
	}
	else{
		for(it=boltzmann->PartMap.begin();it!=boltzmann->PartMap.end();++it){
			parta=it->second;
			if(parta->balanceID<0){
				if(!PPBAR_ONLY || abs(parta->resinfo->pid)==2212){
					ya=atanh(parta->p[3]/parta->p[0]);
					while(ya<-MSU_BOLTZMANN_ETAMAX)
						ya+=2.0*MSU_BOLTZMANN_ETAMAX;
					while(ya>MSU_BOLTZMANN_ETAMAX)
						ya-=2.0*MSU_BOLTZMANN_ETAMAX;
					ppartmap.insert(pair<double,CMSUPart* >(ya,parta));
				}
			}
		}
		ita=ppartmap.begin();
		do{
			parta=ita->second;
			pida=abs(parta->resinfo->pid);
			if(pida==211 || pida==321 || pida==2212){
				pida=parta->resinfo->pid;
				IncrementDenom(parta);
				ya=ita->first;
				itb=ita;
				++itb;
				if(itb==ppartmap.end())
					itb=ppartmap.begin();
				dely=0.0;
				do{
					partb=itb->second;
					pidb=abs(partb->resinfo->pid);
					if(pidb==211 || pidb==321 || pidb==2212){
						yb=itb->first;
						dely=yb-ya;
						if(dely<0.0)
							dely+=2.0*MSU_BOLTZMANN_ETAMAX;
				
						// See Jinjin's thesis, page 31
						pidb=partb->resinfo->pid;
						delymax=acceptance->GetDelYMax(pida,pidb);
						if(dely<delymax){
							IncrementNumer(parta,partb);
						}
					}
					++itb;
					if(itb==ppartmap.end())
						itb=ppartmap.begin();
				}while(dely<2.0*BF_YMAX && dely<MSU_BOLTZMANN_ETAMAX);
			}
			++ita;
		}while(ita!=ppartmap.end());
		ppartmap.clear();
	}
}

void CBalanceArrays::ProcessV2Perfect(){
	double phi,eta,pmag;
	CMSUPart boostedpart;
	CMSUPartMap::iterator it;
	CMSUPart *part;
	for(it=boltzmann->PartMap.begin();it!=boltzmann->PartMap.end();++it){
		part=it->second;
		if(abs(part->resinfo->charge)>0){
			if(part->balanceID<0){
				phi=atan2(part->p[2],part->p[1]);
				v2perfect+=cos(2.0*phi);
				v2perfectnorm+=1.0;
				boostedpart.Copy(part);
				while(boostedpart.y<-acceptance->ETAMAX){
					boostedpart.BoostRap(2.0*acceptance->ETAMAX);
				}
				while(boostedpart.y>acceptance->ETAMAX){
					boostedpart.BoostRap(-2.0*acceptance->ETAMAX);
				}
				pmag=sqrt(boostedpart.p[1]*boostedpart.p[1]+boostedpart.p[2]*boostedpart.p[2]+boostedpart.p[3]*boostedpart.p[3]);
				eta=atanh(boostedpart.p[3]/pmag);
				if(fabs(eta)<acceptance->ETAMAX){
					dNdeta+=0.5/boltzmann->ETAMAX;
				}
			}
		}
	}
}

void CBalanceArrays::IncrementDenom(CMSUPart *part){
	int pid;
	bool accepta;
	double effa,dely,ya,phia;
	Crandy *randy=boltzmann->randy;
	CMSUPart parta;
	pid=part->resinfo->pid;
	int iy=5+floorl(part->y);
	if(abs(pid)==211)
		denom_pi->dNdy+=1.0/(4.0*boltzmann->NSAMPLE*boltzmann->ETAMAX);
	if(abs(pid)==321)
		denom_K->dNdy+=1.0/(4.0*boltzmann->NSAMPLE*boltzmann->ETAMAX);
	if(abs(pid)==2212)
		denom_p->dNdy+=1.0/(4.0*boltzmann->NSAMPLE*boltzmann->ETAMAX);
	if(abs(pid)==211 || abs(pid)==321 || abs(pid)==2212){
		parta.Copy(part);
		ya=BF_YMIN+randy->ran()*(BF_YMAX-BF_YMIN);
		dely=ya-parta.y;
		parta.BoostRap(dely);

		acceptance->CalcAcceptance(accepta,effa,&parta);
		if(accepta){
			if(abs(pid)==211){
				denom_pi->Increment(&parta,effa);
			}
			if(abs(pid)==321)
				denom_K->Increment(&parta,effa);
			if(abs(pid)==2212){
				denom_p->Increment(&parta,effa);
			}
		}
		
		acceptance->CalcAcceptanceNoID(accepta,effa,&parta);	
		if(accepta){	
			if(iy>=0 && iy<10)
				rapdist[iy]+=effa;
			denom_allcharges->Increment(&parta,effa);
			phia=atan2(parta.p[2],parta.p[1]);
			v2+=cos(2.0*phia)*effa;
			v2norm+=effa;
			phia=phia*180.0/PI;
			if(phia<0.0)
				phia=-phia;
			if(phia>90.0)
				phia=180-phia;
			if(phia<BF_PHICUT)
				denom_allcharges_phi0->Increment(&parta,effa);
			if(phia>45.0-BF_PHICUT && phia<45.0+BF_PHICUT)
				denom_allcharges_phi45->Increment(&parta,effa);
			if(phia>90.0-BF_PHICUT)
				denom_allcharges_phi90->Increment(&parta,effa);
		}
	}
}

void CBalanceArrays::IncrementNumer(CMSUPart *parta,CMSUPart *partb){
	double effa,effb,effaNoID,effbNoID,ya,dely,phia,phib,delyb=0.0,Minv;
	bool accepta,acceptb,acceptaNoID,acceptbNoID;
	double MSU_BOLTZMANN_ETAMAX=boltzmann->ETAMAX;
	int pida,pidb;
	Crandy *randy=boltzmann->randy;
	CMSUPart partaa,partbb;
	pida=parta->resinfo->pid;
	pidb=partb->resinfo->pid;
	if(abs(pida)!=211 && abs(pida)!=321 && abs(pida)!=2212){
		return;
	}
	if(abs(pidb)!=211 && abs(pidb)!=321 && abs(pidb)!=2212){
		return;
	}
	if(PPBAR_ONLY && (abs(pida)!=2212 || abs(pidb)!=2212)){
		return;
	}
	
	if(abs(parta->resinfo->charge)==1 && abs(partb->resinfo->charge)==1 ){
		partaa.Copy(parta);
		ya=BF_YMIN+randy->ran()*(BF_YMAX-BF_YMIN);
		dely=ya-partaa.y;
		partaa.BoostRap(dely);
		acceptance->CalcAcceptance(accepta,effa,&partaa);
		acceptance->CalcAcceptanceNoID(acceptaNoID,effaNoID,&partaa);
		if(accepta || acceptaNoID){
			delyb=0.0;
			if(partb->y+dely>MSU_BOLTZMANN_ETAMAX){
				delyb=-2.0*MSU_BOLTZMANN_ETAMAX;
			}
			if(partb->y+dely<-MSU_BOLTZMANN_ETAMAX)
				delyb=2.0*MSU_BOLTZMANN_ETAMAX;
			partbb.Copy(partb);
			partbb.BoostRap(dely+delyb);
			acceptance->CalcAcceptance(acceptb,effb,&partbb);
			acceptance->CalcAcceptanceNoID(acceptbNoID,effbNoID,&partbb);

			if(NoKsNoPhi){
				if(abs(pida)==321 && abs(pidb)==321){
					Minv=GetMinv(&partaa,&partbb);
					if(Minv>1010.0 && Minv<1030)
						accepta=acceptb=false;
				}
				if(abs(pida)==211 && abs(pidb)==211){
					Minv=GetMinv(&partaa,&partbb);
					if(Minv>497.146 && Minv<498.146)
						accepta=acceptb=false;
				}
			}
			if(accepta && acceptb){
				if(abs(pida)==211 && abs(pidb)==211){
					numer_pipi->Increment(&partaa,&partbb,effa,effb);
				}
				if( (abs(pida)==211 && abs(pidb)==321) || (abs(pidb)==211 && abs(pida)==321) ){
					numer_piK->Increment(&partaa,&partbb,effa,effb);
				}
				if( (abs(pida)==211 && abs(pidb)==2212) || (abs(pidb)==211 && abs(pida)==2212) ){
					numer_pip->Increment(&partaa,&partbb,effa,effb);
				}
				if(abs(pida)==321 && abs(pidb)==321){
					numer_KK->Increment(&partaa,&partbb,effa,effb);
				}
				if( (abs(pida)==321 && abs(pidb)==2212) || (abs(pidb)==321 && abs(pida)==2212) ){
					numer_Kp->Increment(&partaa,&partbb,effa,effb);
				}
				if(abs(pida)==2212 && abs(pidb)==2212){
					numer_pp->Increment(&partaa,&partbb,effa,effb);
				}
			}
			
			if(acceptaNoID && acceptbNoID){
				numer_allcharges->Increment(&partaa,&partbb,effaNoID,effbNoID);
				phia=atan2(partaa.p[2],partaa.p[1]);
				phib=atan2(partbb.p[2],partbb.p[1]);
				v2prime+=effaNoID*cos(2.0*phia)+effbNoID*cos(2.0*phib);
				v2primenorm+=effaNoID+effbNoID;
				phia=phia*180.0/PI;
				if(phia<0.0){
					phia=-phia;
				}
				if(phia>90.0){
					phia=180-phia;
				}
				if(phia<BF_PHICUT){
					numer_allcharges_phi0->Increment(&partaa,&partbb,effaNoID,effbNoID);
				}
				if(phia>45.0-BF_PHICUT && phia<45.0+BF_PHICUT){
					numer_allcharges_phi45->Increment(&partaa,&partbb,effaNoID,effbNoID);
				}
				if(phia>90.0-BF_PHICUT){
					numer_allcharges_phi90->Increment(&partaa,&partbb,effaNoID,effbNoID);
				}
				phib=phib*180.0/PI;
				if(phib<0.0){
					phib=-phib;
				}
				if(phib>90.0){
					phib=180-phib;
				}
				if(phib<BF_PHICUT){
					numer_allcharges_phi0->Increment(&partbb,&partaa,effbNoID,effaNoID);
				}
				if(phib>45.0-BF_PHICUT && phib<45.0+BF_PHICUT){
					numer_allcharges_phi45->Increment(&partbb,&partaa,effbNoID,effaNoID);
				}
				if(phia>90.0-BF_PHICUT){
					numer_allcharges_phi90->Increment(&partaa,&partbb,effbNoID,effaNoID);
				}

				IncrementGammaP(&partaa,&partbb,effaNoID,effbNoID);
			}			
		}
	}
}

void CBalanceArrays::SetQualifier(string qualifier_set){
	qualifier=qualifier_set;
	string command="mkdir -p modelruns/"+boltzmann->run_name+"/"+qualifier;
	system(command.c_str());
	if(FROM_UDS){
		bf_results_dirname="modelruns/"+boltzmann->run_name+"/"+qualifier+"/results_type1";
	}
	else{
		bf_results_dirname="modelruns/"+boltzmann->run_name+"/"+qualifier+"/results_type2";
	}
	bf_results_dirname=bf_results_dirname+"/subruns/subrun"+to_string(boltzmann->subrun_number);
	command="mkdir -p "+bf_results_dirname;
	system(command.c_str());

}

double CBalanceArrays::GetMinv(CMSUPart *parta,CMSUPart *partb){
	FourVector *pa=&(parta->p);
	FourVector *pb=&(partb->p);
	FourVector P;
	for(int alpha=0;alpha<4;alpha++)
		P[alpha]=(*pa)[alpha]+(*pb)[alpha];
	double minv=P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3];
	return sqrt(minv);
}