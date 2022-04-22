#include "msu_boltzmann/action.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/cell.h"
#include "msu_sampler/resonances.h"

void CAction::PerformDecay(){
	CMSUPart *mother,*dptr;
	CMSUPartMap::iterator ppos;
	int ibody,nbodies;
	double mt,etamax=boltzmann->ETAMAX,mothermass;
	double deleta;
	printf("howdy\n");
	ppos=partmap.begin();
	mother=ppos->second;
	boltzmann->GetDeadParts(product);

	mothermass=mother->GetMass();
	if(mother->cell!=NULL && mother->cell!=mother->FindCell() && tau<boltzmann->TAUCOLLMAX){
		mother->CheckRapidity();
		mother->cell->Print();
		mother->FindCell()->Print();
		mother->Print();
		sprintf(message,"Cells don't match for decaying mother\n");
		CLog::Fatal(message);
	}
	if(tau>boltzmann->TAUCOLLMAX || mother->cell==NULL){
		while(boltzmann->BJORKEN && (mother->eta<-etamax || mother->eta>etamax)){
			if(mother->eta<-etamax) deleta=2.0*etamax*ceil((-etamax-mother->eta)/(2.0*etamax));
			if(mother->eta>etamax) deleta=-2.0*etamax*ceil((mother->eta-etamax)/(2.0*etamax));
			mother->eta+=deleta;
			mother->y+=deleta;
			mt=mother->GetMT();
			mother->p[0]=mt*cosh(mother->y);
			mother->p[3]=mt*sinh(mother->y);
			mother->r[0]=tau*cosh(mother->eta);
			mother->r[3]=tau*sinh(mother->eta);
		}
	}
	printf("howdy c\n");
	mother->resinfo->DecayGetResInfoPtr(mothermass,nbodies,daughterresinfo);
	printf("howdy d\n");
	for(ibody=0;ibody<nbodies;ibody++){
		product[ibody]->resinfo=daughterresinfo[ibody];
	}
	printf("howdy e\n");
	boltzmann->Decay(mother,nbodies,product);
	printf("howdy f\n");
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for(ibody=0;ibody<nbodies;ibody++){
		dptr=product[ibody];
		dptr->active=true;
		dptr->bweight=mother->bweight;
		dptr->balanceID=mother->balanceID;
		dptr->nscatt=0;
		dptr->tau_lastint=tau;
		dptr->actionmother=boltzmann->nactions;
		dptr->ChangeCell(dptr->FindCell());
		if(dptr->currentmap!=&(boltzmann->PartMap))
			dptr->ChangeMap(&(boltzmann->PartMap));
		dptr->FindActions();
	}
	mother->Kill();
	boltzmann->ndecay+=1;
	printf("bye bye\n");
}
