#ifndef __MSU_BOLTZMANN_ACTION_H__
#define __MSU_BOLTZMANN_ACTION_H__

#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/cell.h"
#include "msu_eos/resonances.h"

using namespace std;
using namespace NMSUPratt;

//!An action in the CMSU_Boltzmann model.
/*!
\version 1.0
\author Scott Pratt
\date March 2011

This class handles any actions that the model takes during execution. Examples of "actions" that the model takes are a resonance decaying, a particle crossing a cell boundary, a collision, new particles being generated, etc. In this way, a complex system of interacting particles is reduced to a scheduled list of actions. Scheduling is handled using a C++ map container of CAction objects, keyed by the boost-invariant time tau (\f$\tau\f$) at which they are scheduled to occur. Note that this map is revised consistently, as future actions often change dramatically as a result of the current action.

Actions are allocated and are moved from the map of future actions (CMSU_Boltzmann::ActionMap) to the list of completed actions (CMSU_Boltzmann::DeadActionMap) once they have been performed.
*/

namespace NMSUPratt{

	class CAction{
	public:
		double tau;
		double pibsquared,sigma_scatter,sigma_merge,sigma_annihilation,sigma_inel;
		vector<double> dsigma_merge;
		int listid;
		double key;
		int type; // =0 for activation, 1 for decay, 2 for collision, ....  6 for ExitCell
		// These are the particles in the action
		CMSUPartMap partmap;

		bool Kill();
		void AddPart(CMSUPart *partptr);
		void Print();

		void Perform();
		void PerformDensCalc();
		void PerformMuTCalcUpdateNPE();
		void PerformSECalc();
		void PerformActivate();
		void PerformExitCell();
		void PerformDecay();
		void PerformCollide();
		void PerformCollide_BALANCE();
		void PerformResetCollisions();
		array <CresInfo *,5> daughterresinfo;
		//void PerformSwallowParticles();
		CAction();
		void InitDead(int keyset);
		CAction(int keyset);
		~CAction();

		static CMSU_Boltzmann *boltzmann;
		static char *message;

		CActionMap::iterator GetPos(CActionMap *actionmap);
		void MoveToActionMap();
		void RemoveFromActionMap();
		void AddToMap(CActionMap *newmap);
		void AddToMap(CActionMap::iterator guess,CActionMap *newmap);
		void CheckPartList();
		CActionMap *currentmap;
		array<CMSUPart *,5> product;
	};

}

#endif