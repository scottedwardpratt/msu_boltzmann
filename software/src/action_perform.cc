#include "msu_boltzmann/action.h"
#include "msu_boltzmann/msupart.h"

using namespace std;
using namespace NMSUPratt;

void CAction::Perform(){
	CMSUPartMap::iterator ppos;
	CMSUPart *part;
	CMSU_BoltzmannCell *cell;
	//if(tau>20 && type!=1){
		//CLog::Info("Performing Action of type "+to_string(type)+", tau="+to_string(tau)+"\n");
		//CLog::Info("nactionstot="+to_string(boltzmann->nactionstot)+"\n");
	//}


	boltzmann->nactionstot+=1;
	boltzmann->tau=tau;
	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Propagate(tau);
		if(part->currentmap != &(boltzmann->PartMap)){
			CLog::Info("wrong map for part in Action List\n");
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
		snprintf(message,CLog::CHARLENGTH,"FATAL:: action earlier than tau!!!!, boltzmann->tau=%15.10e, action tau=%15.10e\n",boltzmann->tau,tau);
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
		snprintf(message,CLog::CHARLENGTH,"FATAL: action type = %d is unknown, exiting\n",type);
		CLog::Fatal(message);
	}
	Kill();
	//CLog::Info("action performed\n");
}
