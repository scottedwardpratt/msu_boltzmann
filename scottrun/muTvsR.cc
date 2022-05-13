#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

int main(){
	int ix,iy,NR=24,ir,N;
	double delR=1.0,DXY=1.0,x,y,tau,r;
	double Ux,Uy,mu,T,rho,epsilon;
	double Ux_alt,Uy_alt;

	FILE *fptr,*output;
	vector<double> muB,muK,mupi,Tpi,TK,TB;
	vector<double> Upi,UK,UB,Upi_alt,UK_alt,UB_alt;
	vector<double> rhopi,rhoK,rhoB,epi,eK,eB;

	muB.resize(NR,0.0); muK.resize(NR,0.0); mupi.resize(NR,0.0);
	TB.resize(NR,0.0); TK.resize(NR,0.0); Tpi.resize(NR,0.0);
	Upi.resize(NR,0.0); UK.resize(NR,0.0); UB.resize(NR,0.0);
	Upi_alt.resize(NR,0.0); UK_alt.resize(NR,0.0),UB_alt.resize(NR,0.0);
	rhopi.resize(NR,0.0); rhoK.resize(NR,0.0); rhoB.resize(NR,0.0);
	epi.resize(NR,0.0); eK.resize(NR,0.0); eB.resize(NR,0.0);
	vector<int> npts,Nparts;
	npts.resize(NR,0);
	Nparts.resize(NR,0);
	char filename[200],dummy[200];

	for(tau=1;tau<30.001;tau+=1){
		

		//Pions
		for(ir=0;ir<NR;ir++){
			npts[ir]=0;
			Nparts[ir]=0;
			mupi[ir]=Tpi[ir]=Upi[ir]=rhopi[ir]=epi[ir]=0.0;
		}

		sprintf(filename,"mutcalc_results/mutinfo_pi_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
			if(!feof(fptr)){
				x=(ix+0.5)*DXY;
				y=(iy+0.5)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(floor(r/delR));
				if(ir<NR && N>0){
					npts[ir]+=1;
					mupi[ir]+=mu;
					Tpi[ir]+=T;
					Upi[ir]+=sqrt(Ux*Ux+Uy*Uy);
					Upi_alt[ir]+=sqrt(Ux_alt*Ux_alt+Uy_alt*Uy_alt);
					rhopi[ir]+=rho;
					epi[ir]+=epsilon;
					Nparts[ir]+=N;
				}
			}
		}while(!feof(fptr));
		fclose(fptr);

		sprintf(filename,"mutcalc_results/muTvsR_pi_tau%g.txt",tau);
		output=fopen(filename,"w");
		for(ir=0;ir<NR;ir++){
			if(npts[ir]>0){
				mupi[ir]=mupi[ir]/double(npts[ir]);
				Upi[ir]=Upi[ir]/double(npts[ir]);
				Upi_alt[ir]=Upi_alt[ir]/double(npts[ir]);
				Tpi[ir]=Tpi[ir]/double(npts[ir]);
				rhopi[ir]=rhopi[ir]/double(npts[ir]);
				epi[ir]=epi[ir]/double(npts[ir]);
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.5f %7.4f %7.5f %7.5f %7.5f %7.5f\n",r,Nparts[ir],Tpi[ir],Upi[ir],mupi[ir],rhopi[ir],epi[ir],Upi_alt[ir]);
			}
			else{
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.5f %7.4f %7.5f %7.4f %7.5f %7.5f\n",r,0,0.0,0.0,0.0,0.0,0.0,0.0);
			}
			npts[ir]=0;
		}
		fclose(output);

		// Kaons
		for(ir=0;ir<NR;ir++){
			npts[ir]=0;
			Nparts[ir]=0;
			muK[ir]=TK[ir]=UK[ir]=rhoK[ir]=eK[ir]=0.0;
		}
		
		sprintf(filename,"mutcalc_results/mutinfo_K_tau%g.txt",tau);
		fptr=fopen(filename,"r");
		fgets(dummy,200,fptr);
		do{
			fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);

			if(!feof(fptr)){
				x=(ix+0.5)*DXY;
				y=(iy+0.5)*DXY;
				r=sqrt(x*x+y*y);
				ir=lrint(floor(r/delR));
				if(ir<NR && N>0){
					//printf("ir=%d\n",ir);
					npts[ir]+=1;
					muK[ir]+=mu;
					TK[ir]+=T;
					UK[ir]+=sqrt(Ux*Ux+Uy*Uy);
					UK_alt[ir]+=sqrt(Ux_alt*Ux_alt+Uy_alt*Uy_alt);
					rhoK[ir]+=rho;
					eK[ir]+=epsilon;
					Nparts[ir]+=N;
				}
			}
		}while(!feof(fptr));
		fclose(fptr);

		sprintf(filename,"mutcalc_results/muTvsR_K_tau%g.txt",tau);
		output=fopen(filename,"w");
		for(ir=0;ir<NR;ir++){
			if(npts[ir]>0){
				muK[ir]=muK[ir]/double(npts[ir]);
				UK[ir]=UK[ir]/double(npts[ir]);
				UK_alt[ir]=UK_alt[ir]/double(npts[ir]);
				TK[ir]=TK[ir]/double(npts[ir]);
				rhoK[ir]=rhoK[ir]/double(npts[ir]);
				eK[ir]=eB[ir]/double(npts[ir]);
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.5f %7.4f %7.5f %7.5f %7.5f %7.5f\n",r,Nparts[ir],TK[ir],UK[ir],muK[ir],rhoK[ir],eK[ir],UK_alt[ir]);
			}
			else{
				r=(ir+0.5)*delR;
				fprintf(output,"%6.2f %7d %7.5f %7.4f %7.5f %7.5f %7.5f %7.5f\n",r,0,0.0,0.0,0.0,0.0,0.0,0.0);
			}
			npts[ir]=0;
		}
		fclose(output);

		// Baryons and Hyperons

		for(int btype=0;btype<8;btype++){
			for(ir=0;ir<NR;ir++){
				npts[ir]=0;
				Nparts[ir]=0;
				muB[ir]=TB[ir]=UB[ir]=rhoB[ir]=eB[ir]=0.0;
			}
			sprintf(filename,"mutcalc_results/mutinfo_B%d_tau%g.txt",btype,tau);
			fptr=fopen(filename,"r");
			fgets(dummy,200,fptr);
			do{
				fscanf(fptr,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&ix,&iy,&N,&T,&Ux,&Uy,&mu,&rho,&epsilon,&Ux_alt,&Uy_alt);
				if(!feof(fptr)){
					x=(ix+0.5)*DXY;
					y=(iy+0.5)*DXY;
					r=sqrt(x*x+y*y);
					ir=lrint(floor(r/delR));
					if(ir<NR && N>0){
						npts[ir]+=1;
						muB[ir]+=mu;
						TB[ir]+=T;
						UB[ir]+=sqrt(Ux*Ux+Uy*Uy);
						UB_alt[ir]+=sqrt(Ux_alt*Ux_alt+Uy_alt*Uy_alt);
						rhoB[ir]+=rho;
						eB[ir]+=epsilon;
						Nparts[ir]+=N;
					}
				}
			}while(!feof(fptr));
			fclose(fptr);


			sprintf(filename,"mutcalc_results/muTvsR_B%d_tau%g.txt",btype,tau);
			output=fopen(filename,"w");
			for(ir=0;ir<NR;ir++){
				if(npts[ir]>0){
					muB[ir]=muB[ir]/double(npts[ir]);
					TB[ir]=TB[ir]/double(npts[ir]);
					UB[ir]=UB[ir]/double(npts[ir]);
					UB_alt[ir]=UB_alt[ir]/double(npts[ir]);
					rhoB[ir]=rhoB[ir]/double(npts[ir]);
					eB[ir]=eB[ir]/double(npts[ir]);
					r=(ir+0.5)*delR;
					fprintf(output,"%6.2f %7d %7.5f %7.4f %7.5f %7.5f %7.5f %7.5f\n",r,Nparts[ir],TB[ir],UB[ir],muB[ir],rhoB[ir],eB[ir],UB_alt[ir]);
				}
				else{
					r=(ir+0.5)*delR;
					fprintf(output,"%6.2f %7d %7.5f %7.4f %7.5f %7.5f %7.5f %7.5f\n",r,0,0.0,0.0,0.0,0.0,0.0,0.0);
				}
				npts[ir]=0;
			}
			fclose(output);
		}

	}
	return 0;
}
