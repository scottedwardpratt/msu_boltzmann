#include "msu_boltzmann/action.h"
#include "msu_boltzmann/msupart.h"

void CAction::Perform(){
	CMSUPartMap::iterator ppos;
	CMSUPart *part;
	CMSU_BoltzmannCell *cell;

	//sprintf(message,"Performing Action of type %d\n",type);
	//CLog::Info(message);

	boltzmann->nactions+=1;
	Kill();
	boltzmann->tau=tau;

	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Propagate(tau);
		if(part->currentmap != &(boltzmann->PartMap)){
			sprintf(message,"wrong map for part in Action List\n");
			CLog::Info(message);
			part->Print();
			part->active=true;
			part->ChangeMap(&(boltzmann->PartMap));
			cell=part->FindCell();
			if(part->cell!=cell){
				part->ChangeCell(cell);
			}
		}
	}
	if(tau+1.0E-4<boltzmann->tau){
		sprintf(message,"FATAL:: action earlier than tau!!!!, boltzmann->tau=%15.10e, action tau=%15.10e\n",boltzmann->tau,tau);
		CLog::Fatal(message);
	}
	if(type==6){
		PerformExitCell();
	}
	else if(type==0){
		PerformActivate();
	}
	else if(type==1){
		PerformDecay();
	}
	else if(type==2){
		PerformCollide();
	}
	else if(type==4){
		PerformDensCalc();
	}
	else if(type==5){
		PerformMuTCalcUpdateNPE();
	}
	else{
		sprintf(message,"FATAL: action type = %d is unknown, exiting\n",type);
		CLog::Fatal(message);
	}

	
}
