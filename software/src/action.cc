#include  "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/action.h"

using namespace std;
using namespace NMSUPratt;

CMSU_Boltzmann *CAction::boltzmann=NULL;
char *CAction::message=new char[500];

CAction::CAction(){
}

CAction::CAction(int keyset){
	InitDead(keyset);
}

CAction::~CAction(){
}

void CAction::InitDead(int keyset){
	key=keyset;
	listid=key;
	tau=0.0;
	type=-1;
	currentmap=&boltzmann->DeadActionMap;
	boltzmann->DeadActionMap.insert(CActionPair(key,this));
	partmap.clear();
	boltzmann->nactionstot+=1;
}

// type=0(creation) 1(decay) 2(collision) 3(VizWrite) 4(DensCalc)

CActionMap::iterator CAction::GetPos(CActionMap *emap){
	pair<CActionMap::iterator,CActionMap::iterator> epospair;
	CActionMap::iterator epos;
	epospair=emap->equal_range(key);
	epos=epospair.first;
	while(epos->second!=this && epos!=epospair.second){
		++epos;
	}
	if(epos->second!=this){
		snprintf(message,CLog::CHARLENGTH,"CAction::GetPos cannot find this action, key=%g\n",key);
		CLog::Info(message);
		return emap->end();
	}
	else
		return epos;
}

void CAction::MoveToActionMap(){
	CActionMap::iterator epos,eepos;
	if(currentmap==&(boltzmann->ActionMap)){
		snprintf(message,CLog::CHARLENGTH,"trying to move action to ActionMap even though action is already in ActionMap\n");
		snprintf(message,CLog::CHARLENGTH,"%sWrong current map\n",message);
		CLog::Fatal(message);
	}
	if(currentmap==&(boltzmann->DeadActionMap)){
		epos=GetPos(currentmap);
		if(epos!=currentmap->end()){
			boltzmann->DeadActionMap.erase(epos);
			//partmap.clear();
		}
		else{
			snprintf(message,CLog::CHARLENGTH,"cannot find epos for action in DeadActionMap!!!\n");
			snprintf(message,CLog::CHARLENGTH,"%sDeadActionMap.size=%d\n",message,int(boltzmann->DeadActionMap.size()));
			CLog::Fatal(message);
		}
		key=tau;
		AddToMap(&(boltzmann->ActionMap));
	}
}

bool CAction::Kill(){
	CMSUPart *part;
	CActionMap::iterator eepos,epos;
	CMSUPartMap::iterator ppos;
	epos=GetPos(currentmap);
	if(epos==boltzmann->ActionMap.end()){
		key=0;
		ppos=partmap.begin();
		while(ppos!=partmap.end()){
			part=ppos->second;
			epos=part->actionmap.begin();
			while(epos!=part->actionmap.end()){
				eepos=epos;
				++eepos;
				if(epos->second==this)
					part->actionmap.erase(epos);
				epos=eepos;
			}
			++ppos;
		}
		partmap.clear();
		dsigma_merge.clear();
		currentmap=&(boltzmann->DeadActionMap);
		return false;
	}
	else{
		currentmap->erase(epos);
		boltzmann->nactionkills+=1;
		key=0;
		AddToMap(boltzmann->DeadActionMap.end(),&boltzmann->DeadActionMap);
		//AddToMap(&boltzmann->DeadActionMap);
		ppos=partmap.begin();
		while(ppos!=partmap.end()){
			part=ppos->second;
			epos=part->actionmap.begin();
			while(epos!=part->actionmap.end()){
				eepos=epos;
				++eepos;
				if(epos->second==this)
					part->actionmap.erase(epos);
				epos=eepos;
			}
			++ppos;
		}
		dsigma_merge.clear();
		partmap.clear();
		return true;
	}
}

void CAction::AddToMap(CActionMap *newmap){
	currentmap=newmap;
	newmap->insert(CActionPair(key,this));
}

void CAction::AddToMap(CActionMap::iterator guess,CActionMap *newmap){
	currentmap=newmap;
	newmap->insert(guess,CActionPair(key,this));
}

void CAction::AddPart(CMSUPart *part){
	partmap.insert(CMSUPartPair(part->key,part));
}

void CAction::Print(){
	snprintf(message,CLog::CHARLENGTH,"___________ type=%d, tau=%g, nparts=%d ___________\n",type,tau,int(partmap.size()));
	CLog::Info(message);
	CMSUPartMap::iterator ppos;
	CMSUPart *part;
	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Print();
	}
	snprintf(message,CLog::CHARLENGTH,"_________________________________________________\n");
	CLog::Info(message);
}

void CAction::CheckPartList(){
	CMSUPart *part;
	CMSUPartMap::iterator ppos,ppos2;
	ppos=partmap.begin();
	while(ppos!=partmap.end()){
		part=ppos->second;
		ppos2=part->GetPos(&(boltzmann->PartMap));
		if(ppos2==boltzmann->PartMap.end()){
			part->Print();
			snprintf(message,CLog::CHARLENGTH,"%s ____________ CAction::CheckPartList FATAL, action type=%d ________________\n",message,type);
			snprintf(message,CLog::CHARLENGTH,"%s iterator not in expected pmap\n",message);
			CLog::Fatal(message);
		}
		++ppos;
	}
}

void CAction::PerformDensCalc(){
	int itau;
	itau=lrint(floor((tau-0.001)/boltzmann->DENSWRITE_DELTAU));
	if(itau>=boltzmann->DENSWRITE_NTAU){
		snprintf(message,CLog::CHARLENGTH,"trying to perform CAction::DensCalc() for itau>=DENSWRITE_NTAU, =%d",itau);
		CLog::Fatal(message);
	}
	CMSU_BoltzmannCell *cell;
	int ix,iy,ieta;
	for(ix=0;ix<2*boltzmann->NXY;ix++){
		for(iy=0;iy<2*boltzmann->NXY;iy++){
			for(ieta=0;ieta<2*boltzmann->NETA;ieta++){
				cell=boltzmann->cell[ix][iy][ieta];
				cell->dens[itau]+=cell->partmap.size();
			}
		}
	}
}

