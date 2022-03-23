#include "action.h"
#include "part.h"

void CAction::PerformActivate(){
	CPart *part;
	CPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	part->active=true;
	part->ChangeMap(&(b3d->PartMap));
	part->ChangeCell(part->FindCell());
	part->CyclicReset();
	part->tau_lastint=tau;
	part->actionmother=b3d->nactions;
	part->FindActions();
	b3d->nactivate+=1;
	part->SetMass();
}
