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
	int nmerge,nscatter,nevents,nparts,ievent,iqual;
	//char logfilename[100];
	//sprintf(logfilename,"msuboltz_log.txt");
	printf("check aaa\n");
	//CLog::Init(logfilename);
	CLog::INTERACTIVE=true;
	printf("check aa\n");

	string filename="model_output/fixed_parameters.txt";
	parmap.ReadParsFromFile(filename);
	printf("check a\n");
	filename="model_output/"+run_name+"/parameters.txt";
	parmap.ReadParsFromFile(filename);
	printf("check bbb\n");
	CmasterSampler ms(&parmap);
	printf("check bb\n");
	CpartList *pl=new CpartList(&parmap,ms.reslist);

	ms.partlist=pl;
	ms.randy->reset(time(NULL));
	printf("check b\n");
	ms.ReadHyper_OSU_2D();
	printf("check c\n");

	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_name,&parmap,ms.reslist);
	msuboltz->InitCascade();
	printf("check d\n");
	
	nparts=0;
	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	nevents=parmap.getI("MSU_BOLTZMANN_NEVENTSMAX",10);

	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		nmerge=nscatter=0;
		printf("check e, nevents=%d\n",nevents);
		msuboltz->ReadMuTInfo();
		for(ievent=0;ievent<nevents;ievent++){
			printf("ievent=%d\n",ievent);
			msuboltz->Reset();
			nparts+=ms.MakeEvent();
			printf("check f\n");
			msuboltz->InputPartList(pl);
			pl->Clear();
			printf("check g\n");
			msuboltz->PerformAllActions();
			printf("check h\n");
			nmerge+=msuboltz->nmerge;
			nscatter+=msuboltz->nscatter;
			sprintf(message,"ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
			CLog::Info(message);
		}
		sprintf(message,"nmerge/event=%g, nscatter/event=%g\n",
			double(nmerge)/double(nevents),double(nscatter)/double(nevents));
		CLog::Info(message);
		msuboltz->WriteMuTInfo();
	}

	CLog::Info(string("YIPPEE!!!!! We made it all the way through!\n"));
	return 0;
}
