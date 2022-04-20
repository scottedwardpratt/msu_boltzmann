#include "msu_boltzmann/action.h"
#include "msu_boltzmann/msupart.h"

void CAction::PerformActivate(){
	CMSUPart *part;
	CMSUPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	part->active=true;
	printf("precheck\n");
	part->ChangeMap(&(boltzmann->PartMap));
	part->ChangeCell(part->FindCell());
	part->CyclicReset();
	part->tau_lastint=tau;
	part->actionmother=boltzmann->nactions;
	printf("check a\n");
	part->FindActions();
	printf("check b\n");
	boltzmann->nactivate+=1;
	part->SetMass();
}
