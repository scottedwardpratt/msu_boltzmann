#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/log.h"
#include <cstring>
using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		CLog::Info("Usage: msuboltz run_name\n");
		exit(-1);
  }
	CparameterMap parmap;
	string run_name=argv[1];
	char message[200];
	int nmerge,nscatter,nannihilate,ncancel_annihilate,nevents,nparts,ievent,ndecay;
	//char logfilename[100];
	//sprintf(logfilename,"msuboltz_log.txt");
	//CLog::Init(logfilename);
	CLog::INTERACTIVE=true;

	string filename="model_output/fixed_parameters.txt";
	parmap.ReadParsFromFile(filename);
	filename="model_output/"+run_name+"/parameters.txt";
	parmap.ReadParsFromFile(filename);
	CmasterSampler ms(&parmap);
	CpartList *pl=new CpartList(&parmap,ms.reslist);

	ms.partlist=pl;
	ms.randy->reset(time(NULL));
	ms.ReadHyper_OSU_2D();

	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_name,&parmap,ms.reslist);
	msuboltz->InitCascade();

	CBalanceArrays *barray=msuboltz->balancearrays;
	
	nparts=0;
	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	nevents=parmap.getI("MSU_BOLTZMANN_NEVENTSMAX",10);

	nmerge=nscatter=nannihilate=ncancel_annihilate=ndecay=0;
	msuboltz->ReadMuTInfo();
	msuboltz->nevents=0;

	CQualifiers qualifiers;
	int iqual=0;
	qualifiers.Read("qualifiers.txt");
	msuboltz->SetQualifier(qualifiers.qualifier[iqual]->qualname);
	qualifiers.SetPars(msuboltz->parmap,iqual);

	for(ievent=0;ievent<nevents;ievent++){
		msuboltz->Reset();
		nparts+=ms.MakeEvent();
		msuboltz->InputPartList(pl);
		pl->Clear();
		msuboltz->PerformAllActions();
		msuboltz->IncrementHadronCount();
		
		nmerge+=msuboltz->nmerge;
		nscatter+=msuboltz->nscatter;
		nannihilate+=msuboltz->nannihilate;
		ncancel_annihilate+=msuboltz->ncancel_annihilate;
		ndecay+=msuboltz->ndecay;
		sprintf(message,"ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
		CLog::Info(message);
		barray->ProcessPartMap();
	}
	sprintf(message,"ndecay/event=%g, nmerge/event=%g, nscatter/event=%g\n",
		double(ndecay)/double(nevents),double(nmerge)/double(nevents),double(nscatter)/double(nevents));
	CLog::Info(message);
	sprintf(message,"nannihilate=%g, ncancel_annihilate=%g\n",
		double(nannihilate)/double(nevents),double(ncancel_annihilate)/double(nevents));
	CLog::Info(message);
	//msuboltz->WriteMuTInfo();
	msuboltz->WriteHadronCount();
	barray->ConstructBFs();
	barray->WriteBFs();
	barray->WriteDenoms();
	barray->WriteGammaP();

	CLog::Info("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
