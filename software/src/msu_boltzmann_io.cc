#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/misc.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/msupart.h"
#include "msu_eos/resonances.h"
#include "msu_boltzmann/cell.h"
#include "msu_sampler/classdefs.h"
#include "msu_sampler/part.h"

using namespace std;

void CMSU_Boltzmann::InputPartList(CpartList *partlist){
	int ipart,balanceID;
	double oldeta,oldrapidity,Et;
	double tau0,eta,mass,rapidity,weight=1.0;
	//int nbaryons=0,nhyperons=0;
	CMSUPart *newpart;
	Cpart *part;
	for(ipart=0;ipart<partlist->nparts;ipart++){
		part=&(partlist->partvec[ipart]);
		tau0=sqrt(part->r[0]*part->r[0]-part->r[3]*part->r[3]);
		oldeta=asinh(part->r[3]/tau0);
		mass=sqrt(part->msquared);
		oldrapidity=atanh(part->p[3]/part->p[0]);
		weight=1.0;
		balanceID=-1;
		newpart=GetDeadPart();
		eta=ETAMAX*(1.0-2.0*randy->ran());
		rapidity=oldrapidity+(eta-oldeta);
		Et=(part->msquared+part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		part->p[3]=Et*sinh(eta);
		part->p[0]=Et*cosh(eta);
		newpart->InitBalance(part->pid,part->r[1],part->r[2],tau0,eta,part->p[1],part->p[2],mass,rapidity,weight,balanceID);
		/*if(newpart->resinfo->baryon!=0)
			nbaryons+=1;
		if(newpart->resinfo->baryon!=0 && newpart->resinfo->strange!=0)
			nhyperons+=1;
			*/
	}
	//printf("nbaryons=%d, nhyperons=%d\n",nbaryons,nhyperons);
}

double CMSU_Boltzmann::WriteOSCAR(int ievent){
	CMSU_BoltzmannBinaryPartInfo bpart;
	double dnchdy=0;
	int ipart;
	CMSUPart *part;
	CMSUPartMap::iterator ppos;
	
	int nparts=PartMap.size();
	sprintf(message,"writing %d particles to %s\n",nparts,oscarfilename.c_str());
	CLog::Info(message);
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"boltzmann output\n");
		}
	}
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parmap->getD("GLAUBER_B",0.0),
	parmap->getD("GLAUBER_B",0.0));
	ppos=PartMap.begin();
	for(ipart=0;ipart<nparts;ipart++){
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->pid;
			bpart.tau=part->tau0;
			bpart.x=part->r[1];
			bpart.y=part->r[2];
			bpart.eta=part->eta;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			bpart.weight=part->weight;
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else
			fprintf(oscarfile,"%5d %5d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %g\n",
		ipart,part->resinfo->pid,part->p[1],part->p[2],part->p[3],part->p[0],sqrt(part->msquared),part->r[1],part->r[2],part->r[3],part->r[0],part->weight);
		if(ppos==PartMap.end()){
			sprintf(message,"ppos shouldn't be here\n");
			CLog::Fatal(message);
		}
		++ppos;
	}
	return dnchdy/(2.0*ETAMAX);
}

void CMSU_Boltzmann::ReadOSCARHeader(){
	int ndead=3,idead;
	char dummy[200];
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.txt";
	if(BINARY_RW)
		oscarfile=fopen(oscarfilename.c_str(),"rb");
	else{
		oscarfile=fopen(oscarfilename.c_str(),"r");
		for(idead=0;idead<ndead;idead++)
			fgets(dummy,200,oscarfile);
	}
}
/*
int CMSU_Boltzmann::ReadOSCAR(int ievent){
	CMSU_BoltzmannBinaryPartInfo bpart;
	CresInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,ID;
	int nparts_read,ipart=0;
	int ievent_read;
	double bmin,bmax; // impact parameter
	CMSUPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			sprintf(message,"trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			CLog::Fatal(message);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			tau0=bpart.tau;
			r[1]=bpart.x;
			r[2]=bpart.y;
			eta=bpart.eta;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			resinfo=reslist->GetResInfoPtr(ID);
			mass=resinfo->mass;
			weight=bpart.weight;
		}
		else{
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
			&ipart,&ID,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0],&weight);
			tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
			eta=asinh(r[3]/tau0);
			rapidity=asinh(p[3]/p[0]);
		}
		mother->Init(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight);
	}
	return nparts_read;
}*/

double CMSU_Boltzmann::WriteBalanceParts(int ievent){
	CMSU_BoltzmannBinaryBalancePartInfo bpart;
	double sigma=0;
	int nsigma=0;
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"boltzmann output\n");
		}
	}
	int nparts=PartMap.size();
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parmap->getD("GLAUBER_B",0.0),
	parmap->getD("GLAUBER_B",0.0));
	double dnchdy=0,rapidity;
	int ipart;
	CMSUPart *part;
	CMSUPartMap::iterator ppos;
	ipart=0;
	ppos=PartMap.begin();
	do{
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->pid;
			bpart.balanceID=part->balanceID;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			sigma+=part->y*part->y;
			nsigma+=1;
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else{
			rapidity=0.5*log((part->p[0]+part->p[3])/(part->p[0]-part->p[3]));
			if(bpart.balanceID!=-1){
				fprintf(oscarfile,"%5d %5d %7d %12.6e %12.6e %12.6e\n",
				ipart,part->resinfo->pid,part->balanceID,part->p[1],part->p[2],rapidity);
			}
		}
		++ppos;
		ipart+=1;
	} while(ipart<nparts);
	sprintf(message,"WriteBalance -- sigma=%g\n",sqrt(sigma/double(nsigma)));
	CLog::Info(message);
	return dnchdy/(2.0*ETAMAX);
}

int CMSU_Boltzmann::ReadBalanceParts(int ievent){
	CMSU_BoltzmannBinaryBalancePartInfo bpart;
	CresInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,ID,balanceID;
	int nparts_read,ipart=0;
	int ievent_read;
	double bmin,bmax; // impact parameter
	CMSUPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			sprintf(message,"trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			CLog::Fatal(message);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			balanceID=bpart.balanceID;
		}
		else{
			sprintf(message,"should only work for binary rw\n");
			CLog::Fatal(message);
			//fscanf(oscarfile,"%d %d %lf %lf %lf",&ipart,&ID,&p[1],&p[2],&rapidity);
		}
		r[1]=r[2]=0.0;
		eta=rapidity;
		while(fabs(eta)>ETAMAX)
			eta-=2.0*ETAMAX*fabs(eta)/eta;
		tau0=10.0;
		weight=1.0;
		resinfo=reslist->GetResInfoPtr(ID);
		mass=resinfo->mass;
		if(resinfo->charge!=0 || resinfo->decay){
			mother->InitBalance(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight,balanceID);
		}
	}
	return nparts_read;
}

void CMSU_Boltzmann::WriteDens(){
	string densfilename="model_output/"+run_name+"/"+qualifier+"/dens.txt";
	FILE *densfile = fopen(densfilename.c_str(),"w");
	fprintf(densfile,"#ix iy  dens[itau=0] dens[itau=1]...\n");
	double dxy;
	int ix,iy,ieta,iitau;
	for(ix=0;ix<CMuTInfo::NXY;ix++){
		for(iy=0;iy<CMuTInfo::NXY;iy++){
			fprintf(densfile,"%3d %3d",ix,iy);
			for(iitau=0;iitau<DENSWRITE_NTAU;iitau++){
				dxy=0.0;
				for(ieta=0;ieta<2*NETA;ieta++){
					dxy+=cell[ix][iy][ieta]->dens[iitau];
				}
				fprintf(densfile," %6.0f",dxy);
			}
			fprintf(densfile,"\n");
		}
	}
	fclose(densfile);
}

void CMSU_Boltzmann::WriteMuTInfo(){
	char dummy[500];
	string dirname;
	int ix,iy,iitau,ntau,btype;
	double tau_print;
	char filename[50];
	bool sufficientN;
	FILE *fptr;
	CMuTInfo *mti;
	CalcMuTU();
	dirname=parmap->getS("MSU_BOLTZMANN_MUTCALC_OUTPUT_DIR","mutcalc_results");
	ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
	for(iitau=0;iitau<ntau;iitau++){
		tau_print=(iitau+1)*MUTCALC_DELTAU;
		for(ix=0;ix<CMuTInfo::NXY;ix++){
			for(iy=0;iy<CMuTInfo::NXY;iy++){
				mti=muTinfo[iitau][ix][iy];
			}
		}
		sprintf(filename,"%s/mutinfo_pi_tau%g.txt",dirname.c_str(),tau_print);
		fptr=fopen(filename,"w");
		fgets(dummy,500,fptr);
		fprintf(fptr,"#  ix    iy     Npi       Tpi    Uxpi    Uypi     mupi     rho    epsilon   Uxpi_alt    Uyp_alt\n");
		for(ix=0;ix<CMuTInfo::NXY;ix++){
			for(iy=0;iy<CMuTInfo::NXY;iy++){
				mti=muTinfo[iitau][ix][iy];
				sufficientN=true;
				if(mti->Npi<CMuTInfo::NMINCALC)
					sufficientN=false;
				if(sufficientN){
					fprintf(fptr,"%2d %2d %8d %g %g %g %g %g %g %g %g\n",
						ix,iy,mti->Npi,mti->Tpi,mti->Uxpi,mti->Uypi,
						mti->mupi,mti->rhopi,mti->epsilonpi,mti->Uxpi_alt,mti->Uypi_alt);
				}
				else{
					fprintf(fptr,"%d %d %d %g %g %g %g %g %g %g %g\n",ix,iy,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
				}
			}
		}
		fclose(fptr);

		sprintf(filename,"%s/mutinfo_K_tau%g.txt",dirname.c_str(),tau_print);
		fptr=fopen(filename,"w");
		fgets(dummy,500,fptr);
		fprintf(fptr,"#  ix    iy     NK      TK    UxK    UyK      muK      rho      epsilon  UxK_alt   UyK_alt\n");
		for(ix=0;ix<CMuTInfo::NXY;ix++){
			for(iy=0;iy<CMuTInfo::NXY;iy++){
				mti=muTinfo[iitau][ix][iy];
				sufficientN=true;
				if(mti->NK<CMuTInfo::NMINCALC)
					sufficientN=false;
				if(sufficientN){
					fprintf(fptr,"%2d %2d %d %g %g %g %g %g %g %g %g\n",
						ix,iy,mti->NK,mti->TK,mti->UxK,mti->UyK,mti->muK,mti->rhoK,mti->epsilonK,mti->UxK_alt,mti->UyK_alt);
				}
				else{
					fprintf(fptr,"%d %d %8d %g %g %g %g %g %g %g %g\n",ix,iy,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
				}
			}
		}
		fclose(fptr);

		for(btype=0;btype<8;btype++){
			sprintf(filename,"%s/mutinfo_B%d_tau%g.txt",dirname.c_str(),btype,tau_print);
			fptr=fopen(filename,"w");
			fgets(dummy,500,fptr);
			fprintf(fptr,"#  ix    iy     NB      TB    UxB       UyB       muB       rhoB     epsilonB\n");
			for(ix=0;ix<CMuTInfo::NXY;ix++){
				for(iy=0;iy<CMuTInfo::NXY;iy++){
					mti=muTinfo[iitau][ix][iy];
					sufficientN=true;
					if(mti->NB[btype]<CMuTInfo::NMINCALC)
						sufficientN=false;
					if(sufficientN){
						fprintf(fptr,"%d %d %8d %g %g %g %g %g %g %g %g\n",
							ix,iy,mti->NB[btype],mti->TB[btype],mti->UxB[btype],mti->UyB[btype],mti->muB[btype],
							mti->rhoB[btype],mti->epsilonB[btype],mti->UxB_alt[btype],mti->UyB_alt[btype]);
					}
					else{
						fprintf(fptr,"%2d %2d %d %g %g %g %g %g %g %g %g\n",ix,iy,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
					}
				}
			}
			fclose(fptr);
		}

	}
}

void CMSU_Boltzmann::ReadMuTInfo(){
	char dummy[500];
	string dirname;
	int ix,iy,iitau,N,NB,ntau,btype;
	double T,Ux,Uy,Ux_alt,Uy_alt,mu,rho,epsilon;
	double tau_read;
	char filename[60];
	bool READN=false;
	FILE *fptr;
	CMuTInfo *mti;
	ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
	dirname=parmap->getS("MSU_BOLTZMANN_MUTCALC_INPUT_DIR","mutcalc_input");
	for(iitau=0;iitau<ntau;iitau++){
		tau_read=(iitau+1)*MUTCALC_DELTAU;
		sprintf(filename,"%s/mutinfo_pi_tau%g.txt",dirname.c_str(),tau_read);
		fptr=fopen(filename,"r");
		if(fptr){
			fgets(dummy,500,fptr);
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
			while(!feof(fptr)){
				mti=muTinfo[iitau][ix][iy];
				if(N>CMuTInfo::NMINCALC)
					mti->sufficientNpi=true;
				else
					mti->sufficientNpi=false;
				mti->Tpi=T;
				mti->mupi=mu;
				mti->Uxpi=Ux;
				mti->Uypi=Uy;
				mti->Npi=0;
				if(READN){
					mti->Npi=N;
					mti->epsilonpi=epsilon;
				}
				fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
			}
			fclose(fptr);
		}

		sprintf(filename,"%s/mutinfo_K_tau%g.txt",dirname.c_str(),tau_read);
		fptr=fopen(filename,"r");
		if(fptr){
			fgets(dummy,500,fptr);
			fscanf(fptr,"%d %d %d  %lf %lf %lf %lf %lf %lf %lf %lf\n",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
			while(!feof(fptr)){
				mti=muTinfo[iitau][ix][iy];
				if(N>CMuTInfo::NMINCALC)
					mti->sufficientNK=true;
				else
					mti->sufficientNK=false;
				mti->TK=T;
				mti->muK=mu;
				mti->UxK=Ux;
				mti->UyK=Uy;
				mti->NK=0;
				if(READN){
					mti->NK=N;
					mti->epsilonK=epsilon;
				}
				fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
			}
			fclose(fptr);
		}

		for(btype=0;btype<8;btype++){
			sprintf(filename,"%s/mutinfo_B%d_tau%g.txt",dirname.c_str(),btype,tau_read);
			fptr=fopen(filename,"r");
			if(fptr){
				fgets(dummy,500,fptr);
				fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&ix,&iy,&NB,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
				while(!feof(fptr)){
					mti=muTinfo[iitau][ix][iy];
					if(NB>CMuTInfo::NMINCALC)
						mti->sufficientNB[btype]=true;
					else
						mti->sufficientNB[btype]=false;
					mti->TB[btype]=T;
					mti->muB[btype]=mu;
					if(fabs(mu)>0.0001)
					printf("btype=%d, muB=%g\n",btype,mu);
					mti->UxB[btype]=Ux;
					mti->UyB[btype]=Uy;
					mti->NB[btype]=0;
					if(READN){
						mti->NB[btype]=NB;
						mti->epsilonB[btype]=epsilon;
					}
					fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&ix,&iy,&NB,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
				}
				fclose(fptr);
			}
		}

	}
}
