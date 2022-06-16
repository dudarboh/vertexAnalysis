#ifndef VertexAnalysis_h
#define VertexAnalysis_h 1

#include <string>
#include <vector>
#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "EVENT/LCObject.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Vertex.h"
#include "UTIL/Operators.h"
#include "UTIL/LCTOOLS.h"
#include "TH1F.h"


using namespace lcio;
using ROOT::Math::XYZVector;
using ROOT::Math::XYZVectorF;


class VertexAnalysis : public marlin::Processor {
    public:

        VertexAnalysis(const VertexAnalysis&) = delete;
        VertexAnalysis& operator=(const VertexAnalysis&) = delete;
        marlin::Processor* newProcessor() { return new VertexAnalysis; }

        VertexAnalysis();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();

        void prepareRootTree();


        EVENT::MCParticle* getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav);
        void fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& decayChain);
        ReconstructedParticle* getDefaultPfo(EVENT::LCEvent* event, EVENT::ReconstructedParticle* pfo);


    private:
        int _nEvent;
        bool _refitOpt;
        std::string _vtxColName;
        std::string _vtxV0ColName;
        std::string _pfoColName;

        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;

        XYZVector _ip;
        XYZVector _vtxTrue;
        int _nTrueTracks;
        int _nTrueKaons;
        int _nTrueProtons;
        XYZVector _pTrueKaon;
        XYZVector _pTrueProton;
        XYZVectorF _vtxReco;
        double _sigmaX;
        double _sigmaY;
        double _sigmaZ;
        int _nRecoTracks;
        int _nRecoKaons;
        int _nRecoProtons;
        XYZVector _pRecoKaon;
        XYZVector _pRecoProton;
        int _nMissedTracks;
        int _nMissedKaons;
        int _nMissedProtons;
        XYZVector _pMissedKaon;
        XYZVector _pMissedProton;
        int _nConfusedTracks;
        int _nConfusedKaons;
        int _nConfusedProtons;
        XYZVector _pConfusedKaon;
        XYZVector _pConfusedProton;
        bool _isMatched;
        bool _isMissed;
        bool _isConfused;


};

#endif
