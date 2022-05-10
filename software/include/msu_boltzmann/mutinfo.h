#ifndef __MUTINFO_H__
#define __MUTINFO_H__
#include "msu_commonutils/commondefs.h"
#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_boltzmann/cell.h"
using namespace std;


class CMuTInfo{
public:
	CMuTInfo(double tau_set);
	double tau;
	double Pxpi,Pypi,Epi,PxK,PyK,EK;
	double epsilonpi,epsilonK,rhopi,rhoK;
	double Tpi,TK,mupi,muK;
	double Txxpi,Tyypi,Txypi;
	double TxxK,TyyK,TxyK;
	double Uxpi,Uypi,UxK,UyK;
	double Uxpi_alt,Uypi_alt,UxK_alt,UyK_alt;
	int Npi,NK;

	vector<double> epsilonB,rhoB,UxB,UyB,UxB_alt,UyB_alt;
	vector<double> EB,PxB,PyB,TxxB,TyyB,TxyB;
	vector<double> TB,muB;
	vector<int> NB;
	bool sufficientNpi,sufficientNK;
	vector<bool> sufficientNB;
	
	void UpdateNPE(CMSU_BoltzmannCell *cell);
	void CalcAllMuTU();
	static void GetMuT(double mass,int degen,double rho_target,double epsilon_target,
		double &T,double &mu);
	//static bool GetMuT_Baryon(double rhoB_target,double rhoBS_target,double epsilon_target,
		//double &T,double &muB,double &muBS);
	static void GetEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
		double &Ux,double &Uy,double &epsilon);
	static void TestEpsilonU(double T00,double T0x,double T0y,double Txx,double Tyy,double Txy,
double Ux,double Uy,double epsilon);
	void Print();
	static void GetIxIy(double x,double y,int &ix,int &iy);
	//void MuTCalc_PionsWithBose();
	static int NETEVENTS;
	static int NMINCALC;
	static vector<CresInfo*> Bresinfo;
	static vector<vector<double>> taumin;
	static vector<double> massB;
	static vector<int> degenB;
	static int NXY;
	static double DXY;
	static CMSU_Boltzmann *boltzmann;
	static int GetBtypeOctetDecuplet(int pid);
};

#endif