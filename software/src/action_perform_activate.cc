#include "action.h"
#include "part.h"

void CAction::PerformActivate(){
	CPart *part;
	CPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	part->active=true;
	part->ChangeMap(&(boltzmann->PartMap));
	part->ChangeCell(part->FindCell());
	part->CyclicReset();
	part->tau_lastint=tau;
	part->actionmother=boltzmann->nactions;
	part->FindActions();
	boltzmann->nactivate+=1;
	part->SetMass();
}
