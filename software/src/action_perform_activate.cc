#include "msu_boltzmann/action.h"
#include "msu_boltzmann/msupart.h"

using namespace std;
using namespace NMSUPratt;

void CAction::PerformActivate(){
	CMSUPart *part;
	CMSUPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	part->active=true;
	part->ChangeMap(&(boltzmann->PartMap));
	part->ChangeCell(part->FindCell());
	part->CyclicReset();
	part->tau_lastint=tau;
	part->actionmother=boltzmann->nactionstot;
	part->Setp0();
	part->FindActions();
	boltzmann->nactivate+=1;
}
