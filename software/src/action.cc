#include  "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/action.h"

using namespace std;
using namespace NMSUPratt;

CMSU_Boltzmann *CAction::boltzmann=NULL;
char *CAction::message=new char[500];

CAction::CAction(){
	//int keyset=boltzmann->DeadActionMap.size();
	//printf("initializing with key %d\n",keyset);
	//InitDead(keyset);
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
	if(currentmap!=emap){
		CLog::Info("In CAction::GetPos(), emap!=currentmap\n");
		return emap->end();
	}
	pair<CActionMap::iterator,CActionMap::iterator> epospair;
	CActionMap::iterator epos;
	epospair=emap->equal_range(key);
	epos=epospair.first;
	if(epos==emap->end()){
		return emap->end();
	}
	while(epos->second!=this && epos!=epospair.second){
		++epos;
	}
	if(epos->second!=this){
		snprintf(message,CLog::CHARLENGTH,"CAction::GetPos cannot find this action, key=%g\n",key);
		CLog::Info(message);
		Print();
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
			key=tau;
			AddToMap(&(boltzmann->ActionMap));
			//partmap.clear();
		}
		else{
			snprintf(message,CLog::CHARLENGTH,"In MoveToActionMap(), cannot find epos for action in DeadActionMap!!!\n");
			snprintf(message,CLog::CHARLENGTH,"%sDeadActionMap.size=%d\n",message,int(boltzmann->DeadActionMap.size()));
			CLog::Fatal(message);
		}
	}
}

bool CAction::Kill(){
	CActionMap::iterator eeepos,eepos,epos;
	CMSUPartMap::iterator ppos;
	CMSUPart *part;
	if(currentmap==&(boltzmann->DeadActionMap)){
		CLog::Info("trying to kill dead action aAAAAaaa\n");
		epos=GetPos(&(boltzmann->DeadActionMap));
		if(epos==boltzmann->DeadActionMap.end()){
			printf("currentmap=DeadActionMap, but action not there\n");
			epos=GetPos(&(boltzmann->ActionMap));
			if(epos==boltzmann->ActionMap.end()){
				CLog::Fatal("particle not in ActionMap either\n");
			}
		}
	}
	else if(currentmap==&(boltzmann->ActionMap)){
		epos=GetPos(&(boltzmann->ActionMap));
		if(epos==boltzmann->ActionMap.end()){
			printf("currentmap=ActionMap, but action not there\n");
			eepos=GetPos(&(boltzmann->DeadActionMap));
			if(eepos==boltzmann->DeadActionMap.end()){
				Print();
				printf("ActionMap size=%lu DeadActionMapsize=%lu\n",
				boltzmann->ActionMap.size(),boltzmann->DeadActionMap.size());
				CLog::Info("action not in DeadActionMap either\n");
			}
		}
		else{
			currentmap->erase(epos);
			boltzmann->nactionkills+=1;
			key=listid;
			//AddToMap(boltzmann->DeadActionMap.end(),&boltzmann->DeadActionMap);
			AddToMap(&boltzmann->DeadActionMap);
			ppos=partmap.begin();
			while(ppos!=partmap.end()){
				part=ppos->second;
				eepos=part->actionmap.begin();
				while(eepos!=part->actionmap.end()){
					eeepos=eepos;
					++eeepos;
					if(eepos->second==this)
						part->actionmap.erase(eepos);
					eepos=eeepos;
				}
				++ppos;
			}
			dsigma_merge.clear();
			partmap.clear();
			return true;
		}
	}
	type=-1;
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
	if(currentmap==&(boltzmann->DeadActionMap)){
		CLog::Info("currentmap=DeadActionMap\n");
	}
	else if(currentmap==&(boltzmann->ActionMap)){
		CLog::Info("currentmap=ActionMap\n");
	}
	else{
		CLog::Info("currentmap is neither ActionMap or DeadActionMap\n");
	}
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

