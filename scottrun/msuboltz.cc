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
	int NN=0,Npi=0,NK=0;
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
	
	nparts=0;
	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	nevents=parmap.getI("MSU_BOLTZMANN_NEVENTSMAX",10);

	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	nmerge=nscatter=nannihilate=ncancel_annihilate=ndecay=0;
	msuboltz->ReadMuTInfo();
	msuboltz->nevents=0;
	Npi=NK=NN=0;
	for(ievent=0;ievent<nevents;ievent++){
		msuboltz->Reset();
		nparts+=ms.MakeEvent();
		msuboltz->InputPartList(pl);
		Npi+=pl->CountResonances(211)+pl->CountResonances(-221)+pl->CountResonances(111);
		NK+=pl->CountResonances(321)+pl->CountResonances(-321)+pl->CountResonances(311)+pl->CountResonances(-311);
		NN+=pl->CountResonances(2212)+pl->CountResonances(-2212)+pl->CountResonances(2112)+pl->CountResonances(-2112);
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
	}
	sprintf(message,"Npi=%d, NK=%d, NN=%d, NN/Npi=%g\n",Npi,NK,NN,double(NN)/double(Npi));
	CLog::Info(message);
	sprintf(message,"ndecay/event=%g, nmerge/event=%g, nscatter/event=%g\n",
		double(ndecay)/double(nevents),double(nmerge)/double(nevents),double(nscatter)/double(nevents));
	CLog::Info(message);
	sprintf(message,"nannihilate=%g, ncancel_annihilate=%g\n",
		double(nannihilate)/double(nevents),double(ncancel_annihilate)/double(nevents));
	CLog::Info(message);
	msuboltz->WriteMuTInfo();
	msuboltz->WriteHadronCount();

	CLog::Info("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
