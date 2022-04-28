#ifndef __MSU_BOLTZMANN_DEFS_H__
#define __MSU_BOLTZMANN_DEFS_H__

#include <map>
#include <unordered_map>
#include <vector>

using namespace std;

class CCharge;
class CMSUPart;
class CAction;
class CBranchInfo;
class CMSU_Boltzmann;
class CMSU_BoltzmannCell;
class Cpart;
class CpartList;
class CresInfo;

//typedef unordered_map<long int,CresInfo *> CresInfoMap;
//typedef pair<long int, CresInfo*> CresInfoPair;
typedef vector<CBranchInfo *> CBranchList; //gives branchlist name
typedef multimap<int,CCharge* > CChargeMap;
typedef pair<int,CCharge* > CChargePair;
typedef multimap<int,CMSUPart* > CMSUPartMap;
typedef pair<int,CMSUPart* > CMSUPartPair;
typedef multimap<int,CCharge* > mapic;
typedef multimap<int,CMSUPart* > mapip;
typedef multimap<double,CAction *> CActionMap;
typedef pair<double,CAction*> CActionPair;
typedef pair<int,CCharge* > pairic;

#endif