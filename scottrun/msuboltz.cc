#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/log.h"
#include <cstring>
using namespace std;

int main(){
	CparameterMap parmap;
	char message[200];
	string run_name="default_0";
	int nmerge,nscatter,nannihilate,nevents,nparts,ievent,iqual;
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
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		nmerge=nscatter=nannihilate=0;
		msuboltz->ReadMuTInfo();
		for(ievent=0;ievent<nevents;ievent++){
			msuboltz->Reset();
			nparts+=ms.MakeEvent();
			msuboltz->InputPartList(pl);
			pl->Clear();
			msuboltz->PerformAllActions();
			
			nmerge+=msuboltz->nmerge;
			nscatter+=msuboltz->nscatter;
			nannihilate+=msuboltz->nannihilate;
			sprintf(message,"ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
			CLog::Info(message);
		}
		sprintf(message,"nmerge/event=%g, nscatter/event=%g, nannihilate=%g\n",
			double(nmerge)/double(nevents),double(nscatter)/double(nevents),double(nannihilate)/double(nevents));
		CLog::Info(message);
		msuboltz->WriteMuTInfo();
	}

	CLog::Info("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
