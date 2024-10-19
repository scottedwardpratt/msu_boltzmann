#ifndef __MSUPART_H__
#define __MSUPART_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <map>
#include <unordered_map>
#include <sys/stat.h>
#include "msu_sampler/part.h"
#include "msu_commonutils/commondefs.h"
#include "msu_boltzmann/boltzmanndefs.h"

using namespace std;

namespace NMSUPratt{
	
	class CMSUPart{
	public:
		//These are identical to members of Cpart
		CMSUPart();
		~CMSUPart();
		int pid;
		double msquared;
		FourVector p,r;
		CMSU_BoltzmannCell *cell,*nextcell;
		double tau0,tau_lastint,tauexit,taudecay;
		double y,eta;
		double weight,bweight;
		double eta0,phi0;
		int listid,nscatt,balanceID,key;
		int actionmother; //refers to action from which particle was created
		CresInfo *resinfo;
		bool active;
		CMSUPartMap *currentmap; // PartList for a Cell, or boltzmann->DeadPartList


		void Print();
		double GetMass();
		void AddPart(int pid,FourVector &p,FourVector &r);
		void Setp0();
		void SetMsquared();
		void Boost(FourVector &u);
		void BoostP(FourVector &u);
		void BoostR(FourVector &u);
	
		// Additional below are unique for Boltzmann

		CMSUPart(int keyset);
	
		void InitDead(int keyset);
		void InitBalance(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity,double bweight,int balanceid);
		void Propagate(double tau);
		void FindDecay();
		void FindCellExit();

		void Init(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity,double bweight);
		void SetMass();
		void Copy(CMSUPart *part);
		void CopyPositionInfo(CMSUPart *part);
		void CopyMomentumInfo(CMSUPart *part);
		void CyclicReset();
		void SetInitialKey();
		void CheckRapidity();

		void Kill();
		void FindActions();
		void SubtractAction(CAction *actionptr);
		void AddAction(CAction *actionptr);
		void KillActions();
		void CheckMap(CMSUPartMap *expectedpartmap);
		void CheckCell();
		void ChangeMap(CMSUPartMap *newmap);
		void BjorkenTranslate();
		void BjorkenUnTranslate();
		void CalcDCA(double *dca);
		void BoostRap(double dely);
		//~CMSUPart();

		// These are the actions involving these particles
		CActionMap actionmap;

		CMSU_BoltzmannCell *FindCell();

		static CMSU_Boltzmann *boltzmann;
		static CBalance *cb;
		static char *message;
		double GetEta(double tau);
		double GetPseudoRapidity();
		double GetRapidity();
		double GetMT();
		void SetY();
		void SetEta(double neweta);
		void GetHBTPars(double &t,double &rout,double &rside,double &rlong);
		void FindCollisions();

		CMSUPartMap::iterator GetPos(CMSUPartMap *pmap);
		CMSUPartMap::iterator DeleteFromMap(CMSUPartMap *partmap);
		void ChangeCell(CMSU_BoltzmannCell *newcell);
		void RemoveFromCell();

		void ChangePartMap(CMSUPartMap *newmap);
		CMSUPartMap::iterator DeleteFromCurrentMap();
		void AddToMap(CMSUPartMap *newmap);
		void AddToMap(CMSUPartMap::iterator guess,CMSUPartMap *newmap);
	};

	class CMSU_BoltzmannBinaryPartInfo{
	public:
		int ID;
		double tau,x,y,eta;
		double px,py,rapidity,bweight,weight;
		void Print(){
			printf("ID=%d, tau=%g, x=%g, y=%g, eta=%g\n",ID,tau,x,y,eta);
			printf("px=%g,py=%g, rapidity=%g, bweight=%g\n",px,py,rapidity,weight);
		}
	};

	class CMSU_BoltzmannBinaryBalancePartInfo{
	public:
		int ID,balanceID;
		double px,py,rapidity,bweight;
		void Print(){
			printf("ID=%d, balanceID=%d\n",ID,balanceID);
			printf("px=%g,py=%g, rapidity=%g,bweight=%g\n",px,py,rapidity,bweight);
		}
	};

}

#endif
