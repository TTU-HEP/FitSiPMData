#include "../include/NTupleReader.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

class GetScaleWeights
{
private:
    NTupleReader trSupp_;

    void getScaleWeights(NTupleReader& tr)
    {
        if(trSupp_.getNextEvent() &&
           tr.getVar<unsigned long long>("EvtNum") == trSupp_.getVar<unsigned long long>("EvtNum") &&
           tr.getVar<float>("GenMET")          == trSupp_.getVar<float>("GenMET") )
        {
            tr.registerDerivedVec("LHEScaleWeight", new std::vector(trSupp_.getVec<float>("ScaleWeights")));
        }
        else
        {
            std::cout << tr.getVar<unsigned long long>("EvtNum") << "\t" << trSupp_.getVar<unsigned long long>("EvtNum") << "\t" << tr.getVar<float>("GenMET") << "\t" << trSupp_.getVar<float>("MET") << std::endl;
            THROW_NTREXCEPTION("ERROR: Event mismatch between master and supplamental file!!!!");
        }
    }

public:
    GetScaleWeights(TChain *ch) : trSupp_(ch, {"EvtNum"}) {}
    GetScaleWeights(GetScaleWeights&) = delete;
    GetScaleWeights(GetScaleWeights&& gsw) : trSupp_(std::move(gsw.trSupp_)) {}

    void operator()(NTupleReader& tr) { getScaleWeights(tr); }
};

int main()
{
    char baseFile[]         = "testFile.root";
    char supplamantalFile[] = "testFile.root";
    char treeName[]         = "TreeMaker2/PreSelection";
    std::string exampleVar  = "NJets";

    TChain *chBase = new TChain(treeName);
    chBase->Add(baseFile);

    TChain *chSupp = new TChain(treeName);
    chSupp->Add(supplamantalFile);

    try
    {
        NTupleReader tr(chBase, {exampleVar});

        //For NTupleReader users you simply need to register the class with NTupleReader 
        tr.emplaceModule<GetScaleWeights>(chSupp);

        while(tr.getNextEvent())
        {
            const auto& LHEScaleWeight = tr.getVec<float>("LHEScaleWeight");

            std::cout << LHEScaleWeight.size() << std::endl;
        }
    }
    catch(const NTRException& e)
    {
        e.print();
    }

    return 0;
}

