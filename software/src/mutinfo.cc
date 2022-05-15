#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/constants.h"
#include "msu_boltzmann/mutinfo.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/cell.h"
#include "msu_sampler/resonances.h"
//using namespace MSU_EOS;

CMSU_Boltzmann *CMuTInfo::boltzmann=NULL;
int CMuTInfo::NETEVENTS=0;
int CMuTInfo::NMINCALC=5;
int CMuTInfo::NXY=0;
double CMuTInfo::DXY=0.0;
vector<vector<double>> CMuTInfo::taumin{};
vector<CresInfo *> CMuTInfo::Bresinfo{};
vector<double> CMuTInfo::massB{0.938,1.18937,1.31483,1.11568,1.232,1.385,1.530,1.67243};
vector<int> CMuTInfo::degenB{8,12,8,4,32,24,16,8};

CMuTInfo::CMuTInfo(double tau_set){
	tau=tau_set;
	Epi=EK=0.0;
	Pxpi=Pypi=PxK=PyK=0.0;
	Txxpi=Tyypi=Txypi=0.0;
	TxxK=TyyK=TxyK=0.0;
	Npi=NK=0;
	Tpi=TK=145.0;
	sufficientNpi=sufficientNK=false;

	NB.resize(8);
	TB.resize(8);
	muB.resize(8);
	TxxB.resize(8);
	TyyB.resize(8);
	TxyB.resize(8);
	PxB.resize(8);
	PyB.resize(8);
	EB.resize(8);
	epsilonB.resize(8);
	rhoB.resize(8);
	UxB.resize(8);
	UyB.resize(8);
	UxB_alt.resize(8);
	UyB_alt.resize(8);
	sufficientNB.resize(8);

	for(int btype=0;btype<8;btype++){
		NB[btype]=0;
		EB[btype]=PxB[btype]=PyB[btype]=0.0;
		TxxB[btype]=TyyB[btype]=TxyB[btype]=0.0;
		muB[btype]=0.0;
		TB[btype]=100.0;
		epsilonB[btype]=0.0;
		sufficientNB[btype]=false;
	}
}

void CMuTInfo::Print(){
	char message[500];
	sprintf(message,"-------- MuT Info, tau=%g ----------\n",tau);

	sprintf(message,"%sNpi=%d, Epi/N=%g, Pxpi/Npi=%g, Pypi/Npi=%g\n",
		message,Npi,Epi/Npi,Pxpi/Npi,Pypi/Npi);
	sprintf(message,"%sTpi=%g, mupi=%g\n",message,Tpi,mupi);

	sprintf(message,"%sNK=%d, EK/NK=%g, PxK/NK=%g, PyK/NK=%g\n",
		message,NK,EK/NK,PxK/NK,PyK/NK);
	sprintf(message,"%sTK=%g, muK=%g\n",message,TK,muK);

	CLog::Info(message);

	for(int btype=0;btype<8;btype++){
		sprintf(message,"btype=%d\n",btype);
		sprintf(message,"%sNB=%d, EK/NK=%g, PxK/NK=%g, PyK/NK=%g\n",
			message,NB[btype],EB[btype]/NB[btype],PxB[btype]/NB[btype],PyB[btype]/NB[btype]);
		sprintf(message,"%sTK=%g, muK=%g\n",message,TB[btype],muB[btype]);
		CLog::Info(message);
	}
}

void CMuTInfo::CalcAllMuTU(){
	double Txx,Tyy,Txy,T00,T0x,T0y,gamma,degen;
	int btype;
	char message[200];

	double volume=4.0*tau*2.0*boltzmann->ETAMAX*DXY*DXY*double(NETEVENTS);   // factor or 4 due to combining quadrants

	if(Npi>=NMINCALC){
		Txx=Txxpi/volume;
		Tyy=Tyypi/volume;
		Txy=Txypi/volume;
		T00=Epi/volume;
		T0x=Pxpi/volume;
		T0y=Pypi/volume;
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,Uxpi,Uypi,epsilonpi);
		Uxpi_alt=Pxpi/(PionMassGeV*Npi);
		Uypi_alt=Pypi/(PionMassGeV*Npi);
		//TestEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,Uxpi,Uypi,epsilonpi);
		gamma=sqrt(1.0+Uxpi*Uxpi+Uypi*Uypi);
		rhopi=double(Npi)/(gamma*volume);
		degen=3;
		GetMuT(PionMassGeV,degen,rhopi,epsilonpi,Tpi,mupi);
	}
	else{
		Tpi=-1.0;
		mupi=0.0;
	}

	if(NK>=NMINCALC){
		Txx=TxxK/volume;
		Tyy=TyyK/volume;
		Txy=TxyK/volume;
		T00=EK/volume;
		T0x=PxK/volume;
		T0y=PyK/volume;
		GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,UxK,UyK,epsilonK);
		UxK_alt=PxK/(KaonMassGeV*NK);
		UyK_alt=PyK/(KaonMassGeV*NK);
		gamma=sqrt(1.0+UxK*UxK+UyK*UyK);
		rhoK=double(NK)/(gamma*volume);
		degen=4;
		GetMuT(KaonMassGeV,degen,rhoK,epsilonK,TK,muK);
	}
	else{
		TK=-1.0;
		muK=0.0;
	}
	for(btype=0;btype<8;btype++){
		if(NB[btype]>=NMINCALC){
			Txx=TxxB[btype]/volume;
			Tyy=TyyB[btype]/volume;
			Txy=TxyB[btype]/volume;
			T00=EB[btype]/volume;
			T0x=PxB[btype]/volume;
			T0y=PyB[btype]/volume;
			GetEpsilonU(T00,T0x,T0y,Txx,Tyy,Txy,UxB[btype],UyB[btype],epsilonB[btype]);
			UxB_alt[btype]=PxB[btype]/(massB[btype]*NB[btype]);
			UyB_alt[btype]=PyB[btype]/(massB[btype]*NB[btype]);
			gamma=sqrt(1.0+UxB[btype]*UxB[btype]+UyB[btype]*UyB[btype]);
			rhoB[btype]=double(NB[btype])/(gamma*volume);
			GetMuT(massB[btype],degenB[btype],rhoB[btype],epsilonB[btype],TB[btype],muB[btype]);
			if(TB[btype]!=TB[btype] || muB[btype]!=muB[btype]){
				sprintf(message,"btype=%d: Disaster, m=%g, rho=%g, degen=%d, epsilon=%g\n",
					btype,massB[btype],rhoB[btype],degenB[btype],epsilonB[btype]);
				CLog::Info(message);
				sprintf(message,"NB=%d, T00/rho=%g, epsilon/rho=%g, EB/NB=%g, TB=%g, muB=%g\n",NB[btype],T00/rhoB[btype],epsilonB[btype]/rhoB[btype],EB[btype]/NB[btype],TB[btype],muB[btype]);
				CLog::Fatal(message);
			}
		}
		else{
			TB[btype]=-1.0;
			muB[btype]=0.0;
		}
	}
}

void CMuTInfo::GetEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
double &Ux,double &Uy,double &epsilon){
	double Qx,Qy,dQxdUx,dQxdUy,dQydUx,dQydUy;
	double gamma,Det,dUx,dUy;
	double A,B,C,dAdUx,dAdUy,dBdUx,dBdUy,dCdUx,dCdUy,dUmag;
	bool success=false;
	int ntry=0;

	Ux=T0x/T00;
	Uy=T0y/T00;
	do{
		gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
		A=-gamma*T00;
		B=(1.0+gamma/(1.0+gamma))*(T0x*Ux+T0y*Uy);
		C=-(1.0/(1.0+gamma))*(Ux*Txx*Ux+2.0*Ux*Txy*Uy+Uy*Tyy*Uy);
		Qx=gamma*T0x-Ux*Txx-Uy*Txy+(A+B+C)*Ux;
		Qy=gamma*T0y-Uy*Tyy-Ux*Txy+(A+B+C)*Uy;
		if(Qx*Qx+Qy*Qy>1.0E-14){
			dAdUx=A*Ux/(gamma*gamma);
			dAdUy=A*Uy/(gamma*gamma);
			dBdUx=(1.0+(gamma/(1.0+gamma)))*T0x+(Ux/(gamma*(1.0+gamma)*(1.0+gamma)))*(T0x*Ux+T0y*Uy);
			dBdUy=(1.0+(gamma/(1.0+gamma)))*T0y+(Uy/(gamma*(1.0+gamma)*(1.0+gamma)))*(T0x*Ux+T0y*Uy);
			dCdUx=-C*Ux/(gamma*(1.0+gamma))-(2.0/(1.0+gamma))*(Txx*Ux+Txy*Uy);
			dCdUy=-C*Uy/(gamma*(1.0+gamma))-(2.0/(1.0+gamma))*(Tyy*Uy+Txy*Ux);
			dQxdUx=-Txx+(Ux/gamma)*T0x+Ux*(dAdUx+dBdUx+dCdUx);
			dQxdUy=-Txy+(Uy/gamma)*T0x+Ux*(dAdUy+dBdUy+dCdUy);
			dQydUx=-Txy+(Ux/gamma)*T0y+Uy*(dAdUx+dBdUx+dCdUx);
			dQydUy=-Tyy+(Uy/gamma)*T0y+Uy*(dAdUy+dBdUy+dCdUy);
			dQxdUx+=(A+B+C);
			dQydUy+=(A+B+C);
			Det=dQxdUx*dQydUy-dQxdUy*dQydUx;
			dUx=-(dQydUy*Qx-dQxdUy*Qy)/Det;
			dUy=-(-dQydUx*Qx+dQxdUx*Qy)/Det;
			if(fabs(dUx*dUx+dUy*dUy)>0.5){
				dUmag=sqrt(dUx*dUx+dUy*dUy);
				dUx*=0.5/dUmag;
				dUy*=0.5/dUmag;
			}
			Ux+=dUx;
			Uy+=dUy;
			ntry+=1;
		}
		else
			success=true;
	}while(!success && ntry<200);
	if(!success){
		char message[300];
		sprintf(message,"Yikes, No Convergence in CMuTInfo::GetEpsilonU\n");
		sprintf(message,"%sQx=%g, Qy=%g, Ux=%g, Uy=%g\n",message,Qx,Qy,Ux,Uy);
		sprintf(message,"%sT00=%g, T0x=%g, T0y=%g\n",message,T00,T0x,T0y);
		sprintf(message,"%sTxx=%g, Tyy=%g, Txy=%g\n",message,Txx,Tyy,Txy);
		CLog::Fatal(message);
	}
	gamma=sqrt(1.0+Ux*Ux+Uy*Uy);
	epsilon=gamma*T00*gamma-2.0*gamma*T0x*Ux-2.0*gamma*T0y*Uy
	+Ux*Txx*Ux+Uy*Tyy*Uy+2.0*Ux*Txy*Uy;
}

void CMuTInfo::TestEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
	double Ux,double Uy,double epsilon){
	FourVector u,n;
	double udotn,epsilontest=0.0;
	int alpha,beta,gamma,delta;
	udotn=sqrt(1.0+Ux*Ux+Uy*Uy);
	u[0]=udotn;
	u[1]=Ux;
	u[2]=Uy;
	u[3]=0.0;
	n[0]=1.0;
	n[1]=n[2]=n[3]=0.0;
	double g[4],L[4][4],Linv[4][4],T[4][4],Ttilde[4][4],Gtest[4][4];
	g[0]=1.0; g[1]=g[2]=g[3]=-1.0;
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			L[alpha][beta]=Linv[alpha][beta]=T[alpha][beta]=Ttilde[alpha][beta]=Gtest[alpha][beta]=0.0;
		}
	}

	T[0][0]=T00;
	T[0][1]=T[1][0]=T0x;
	T[0][2]=T[2][0]=T0y;
	T[1][1]=Txx;
	T[2][2]=Tyy;
	T[1][2]=T[2][1]=Txy;
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			epsilontest+=u[alpha]*g[alpha]*T[alpha][beta]*g[beta]*u[beta];
			if(alpha==beta){
				L[alpha][beta]+=g[alpha];
				Linv[alpha][beta]+=g[alpha];
			}
			L[alpha][beta]+=2.0*u[alpha]*n[beta];
			Linv[alpha][beta]+=2.0*n[alpha]*u[beta];
			L[alpha][beta]-=(u[alpha]+n[alpha])*(u[beta]+n[beta])/(1.0+udotn);
			Linv[alpha][beta]-=(u[alpha]+n[alpha])*(u[beta]+n[beta])/(1.0+udotn);
		}
	}
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			for(gamma=0;gamma<4;gamma++){
				Gtest[alpha][beta]+=Linv[alpha][gamma]*g[gamma]*L[gamma][beta];
				for(delta=0;delta<4;delta++){
					if(alpha==beta)
					Ttilde[alpha][beta]+=Linv[alpha][gamma]*g[gamma]*T[gamma][delta]*g[delta]*L[delta][beta];
				}
			}
		}
	}
	char message[400];
	sprintf(message,"epsilon=%g =? %g =? %g\n",epsilon,Ttilde[0][0],epsilontest);
	CLog::Info(message);
	sprintf(message,"------ Ttilde --------\n");
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			sprintf(message,"%s%10.5f ",message,1000*Ttilde[alpha][beta]);
		}
		sprintf(message,"%s\n",message);
	}
	CLog::Info(message);
	sprintf(message,"------ Gtest --------\n");
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			sprintf(message,"%s%10.5f ",message,Gtest[alpha][beta]);
		}
		sprintf(message,"%s\n",message);
	}
	CLog::Info(message);
}

void CMuTInfo::GetMuT(double mass,int degen,double rho_target,double epsilon_target,double &T,double &mu){
	double E,dEdT,ETarget,epsilon0,dedT,P,rho0,dT;
	int ntry=0;
	char message[100];
	ETarget=epsilon_target/rho_target;
	if(ETarget<mass+0.01){
		T=0.005;
		mu=ETarget/T;
	}
	else{
		if(ETarget/mass<1.1){
			T=(2.0/3.0)*(ETarget-mass);
			rho0=pow((mass*T)/(2.0*PI*HBARC_GEV*HBARC_GEV),1.5);
		}
		else{
			T=0.09;
			do{
				ntry+=1;
				MSU_EOS::freegascalc_onespecies(T,mass,epsilon0,P,rho0,dedT);
				E=epsilon0/rho0;
				dEdT=dedT/rho0-epsilon0*epsilon0/(rho0*rho0*T*T);
				dT=(ETarget-E)/dEdT;
				if(fabs(dT)>0.4*T)
					dT=0.4*T*dT/fabs(dT);
				T+=dT;
			}while(fabs(dT)>1.0E-5 && ntry<30);
			if(ntry==50 || T!=T){
				sprintf(message,"CMuTInfo::GetMuT did not converge!!!, T=%g, dT=%g\n",
					T,dT);
				CLog::Info("T="+to_string(T)+", ETarget="+to_string(ETarget)+", mass="+to_string(mass)+"\n");
				CLog::Info("rho_target="+to_string(rho_target)+", epsilon_target="+to_string(epsilon_target)+"\n");
				CLog::Fatal(message);
			}
		}
		MSU_EOS::freegascalc_onespecies(T,mass,epsilon0,P,rho0,dedT);
		mu=log(rho_target/(rho0*degen));
	}
	if(T!=T || mu!=mu){
		printf("disaster, T=%g, mu=%g, rho_target=%g, epsilon_target=%g, e/rho=%g\n",T,mu,rho_target,epsilon_target,epsilon_target/rho_target);
		exit(1);
	}
}

void CMuTInfo::GetIxIy(double x,double y,int &ix,int &iy){
	ix=lrint(floor(fabs(x)/DXY));
	iy=lrint(floor(fabs(y)/DXY));
}

int CMuTInfo::GetBtypeOctetDecuplet(int pid){
	int btype=-1;
	pid=abs(pid);
	if(pid==2112 || pid==2212)
		btype=0;
	else if(pid==3222 || pid==3212 || pid==3112)
		btype=1;
	else if(pid==3322 || pid==3312)
		btype=2;
	else if(pid==3122)
		btype=3;
	else if(pid==2224 || pid==2214 || pid==2114 || pid==1114)
		btype=4;
	else if(pid==3224 || pid==3214 || pid==3114)
		btype=5;
	else if(pid==3324 || pid==3314)
		btype=6;
	else if(pid==3334)
		btype=7;
	return btype;
}

