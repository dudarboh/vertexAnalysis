#include "VertexAnalysis.hpp"
#include "marlin/VerbosityLevels.h"
using namespace std;

VertexAnalysis aVertexAnalysis;

struct CompareVectors{
    bool operator()(const XYZVector& a, const XYZVector& b) const { return a.r() < b.r(); }
};

struct CompareVertices{
    bool operator()(const Vertex* a, const Vertex* b) const { return XYZVector(a->getPosition()[0], a->getPosition()[1], a->getPosition()[2]).r() < XYZVector(b->getPosition()[0], b->getPosition()[1], b->getPosition()[2]).r(); }
};


VertexAnalysis::VertexAnalysis():marlin::Processor("VertexAnalysis"){
    _description = "Write to the ROOT file comprehensive information about\
                    reconstucted vs true vertex to analyse impact of refitting the tracks";

    registerProcessorParameter("refit_option",
                              "If true works with verticies collections produced after track refitting processor",
                              _refitOpt,
                              bool(false));

}

void VertexAnalysis::init(){
    _nEvent = 0;
    prepareRootTree();
    if (_refitOpt){
        _vtxColName = "BuildUpVertex_2";
        _vtxV0ColName = "BuildUpVertex_V0_2";
        _pfoColName = "updatedPandoraPFOs";
    }
    else{
        _vtxColName = "BuildUpVertex";
        _vtxV0ColName = "BuildUpVertex_V0";
        _pfoColName = "PandoraPFOs";
    }
}

void VertexAnalysis::processEvent( EVENT::LCEvent* event ){
    ++_nEvent;
    streamlog_out(DEBUG8)<<"===============Event "<<_nEvent<<"==============="<<std::endl;
    // LCCollection* mcCol = event->getCollection("MCParticle");
    const std::vector<std::string>* colNames = event->getCollectionNames();
    std::string mcColName;
    if ( std::find(colNames->begin(), colNames->end(), "MCParticle") != colNames->end() ){
        mcColName = "MCParticle";
    }
    else if ( std::find(colNames->begin(), colNames->end(), "MCParticlesSkimmed") != colNames->end() ){
        mcColName = "MCParticlesSkimmed";
    }
    else{
        return;
    }
    XYZVector trueIP;
    trueIP.SetCoordinates( static_cast<MCParticle*>(event->getCollection(mcColName)->getElementAt(0))->getVertex() );
    // LCCollection* pfoCol = event->getCollection(_pfoColName);


    // Write algorithm that finds all true reconstructable vertices:
    // 1) Create list of all MCParticles that have a track
    std::vector<MCParticle*> mcs;
    LCCollection* pfoCol = event->getCollection("PandoraPFOs");
    LCRelationNavigator navRecoToMc( event->getCollection("RecoMCTruthLink") );
    streamlog_out(DEBUG1)<<"List of MC particles with associated tracks:"<<endl;
    for(int i=0; i < pfoCol->getNumberOfElements(); ++i){
        streamlog_out(DEBUG1)<<"*******i="<<i<<"********"<<endl;
        // streamlog_out(DEBUG)<<pfo;
        ReconstructedParticle* pfo = static_cast<ReconstructedParticle*> ( pfoCol->getElementAt(i) );
        MCParticle* mc = getMcMaxTrackWeight(pfo, navRecoToMc);
        if (mc != nullptr){
            // Sometimes single MCParticle creates two PFOs! Don't write dublicate MCParticle in the list!
            if ( std::find(mcs.begin(), mcs.end(), mc) != mcs.end() ) continue;
            mcs.push_back(mc);
            XYZVector pos;
            pos.SetCoordinates(mc->getVertex());
            streamlog_out(DEBUG1)<<"MC:  "<<mc<<"    ("<<(pos - trueIP).r()<<")"<<"    PDG:"<<mc->getPDG()<<endl;
        }
    }

    // 2) Create a map of these MCParticles to all "True vertices"
    std::map< XYZVector, std::vector<MCParticle*>, CompareVectors > vtx2prongs;
    for(size_t i=0; i < mcs.size(); ++i){    
        MCParticle* mc = mcs[i];
        XYZVector pos( mc->getVertex()[0], mc->getVertex()[1], mc->getVertex()[2] );
        vtx2prongs[pos].push_back(mc);
    }

    streamlog_out(DEBUG2)<<"List of raw true vertices:"<<endl;
    for(auto iter=vtx2prongs.begin(); iter!=vtx2prongs.end(); ++iter){
        streamlog_out(DEBUG2)<<" Vertex #"<<std::distance(vtx2prongs.begin(), iter)<<"  ("<<( (iter->first) - trueIP ).r()<<")  MCs:  ";
        for(auto mc : iter->second) streamlog_out(DEBUG2)<<mc<<"    ";
        streamlog_out(DEBUG2)<<endl;
    }

    // 3) Merge all vertices within 3 um. And recalculate vertex position as a weighted average with N tracks per vertex
    for(auto iter1=vtx2prongs.begin(); iter1!=vtx2prongs.end();){
        bool isMerged = false;
        for(auto iter2=std::next(iter1); iter2!=vtx2prongs.end();){
            streamlog_out(DEBUG1)<<"Beginning merge loop Iter1: "<<std::distance(vtx2prongs.begin(), iter1)<<"/"<<vtx2prongs.size()<<endl;
            streamlog_out(DEBUG1)<<"Beginning merge loop Iter2: "<<std::distance(vtx2prongs.begin(), iter2)<<"/"<<vtx2prongs.size()<<endl;
            XYZVector vtx1 = iter1->first;
            XYZVector vtx2 = iter2->first;
            double distance = (vtx1-vtx2).r();
            //Assuming vertex detector resolution 3um, merge vertices within 3um and recalculate vertex position as a weighted average.
            //NOTE: I merge them in arbitrary order. Practially it doesn't matter, but this may cause a difference in some cases.
            if (distance < 0.003){
                isMerged = true;
                streamlog_out(DEBUG1)<<" Vertices are being merged! Resetting the loop iterators!"<<endl;
                std::vector<MCParticle*> prongs1 = iter1->second;
                int nProngs1 = prongs1.size();
                std::vector<MCParticle*> prongs2 = iter2->second;
                int nProngs2 = prongs2.size();
                XYZVector avgPos = (vtx1*nProngs1 + vtx2*nProngs2)/(nProngs1+nProngs2);
                std::vector<MCParticle*> mergedProngs = prongs1;
                mergedProngs.insert( mergedProngs.end(), prongs2.begin(), prongs2.end() );
                vtx2prongs.erase(vtx1);
                vtx2prongs.erase(vtx2);
                vtx2prongs[avgPos] = mergedProngs;
                break;
            }
            else{
                ++iter2;
            }
        }
        if (isMerged) iter1 = vtx2prongs.begin();
        else ++iter1;
    }
    streamlog_out(DEBUG2)<<"List of true vertices after merging:"<<endl;
    for(auto iter=vtx2prongs.begin(); iter!=vtx2prongs.end(); ++iter){
        streamlog_out(DEBUG2)<<" Vertex #"<<std::distance(vtx2prongs.begin(), iter)<<"  ("<<( (iter->first) - trueIP ).r()<<")  MCs:  ";
        for(auto mc : iter->second) streamlog_out(DEBUG2)<<mc<<"    ";
        streamlog_out(DEBUG2)<<endl;
    }


    // 4) Remove any vertices with only single track as they are unfindable
    // auto hasOneTrack = [](std::pair<XYZVector, std::vector<MCParticle*> > item) {return (item.second).size() < 2;};
    // vtx2prongs.erase(std::remove_if(vtx2prongs.begin(), vtx2prongs.end(), hasOneTrack), vtx2prongs.end());
    for (auto iter = vtx2prongs.begin() ; iter != vtx2prongs.end(); ){
        if ( (iter->second).size() < 2 ) iter = vtx2prongs.erase(iter);
        else ++iter;
    }
    streamlog_out(DEBUG2)<<"List of true vertices after merging and removing 1-track vertices:"<<endl;
    for(auto iter=vtx2prongs.begin(); iter!=vtx2prongs.end(); ++iter){
        streamlog_out(DEBUG2)<<" Vertex #"<<std::distance(vtx2prongs.begin(), iter)<<"  ("<<( (iter->first) - trueIP ).r()<<")  MCs:  ";
        for(auto mc : iter->second) streamlog_out(DEBUG2)<<mc<<"    ";
        streamlog_out(DEBUG2)<<endl;
    }


    // 5) Write prim. vertex position and remove it from the list (closest to the 0,0,0). std::map is already sorted with our comparator.
    if ( vtx2prongs.size() != 0 ){
        auto closestToIp = [&trueIP](const std::pair<XYZVector, std::vector<MCParticle*> >& a, const std::pair<XYZVector, std::vector<MCParticle*> >& b){return (a.first - trueIP).r() < (b.first - trueIP).r();};
        auto ipIter = std::min_element(vtx2prongs.begin(), vtx2prongs.end(), closestToIp);
        _ip = ipIter->first;
        vtx2prongs.erase( ipIter );
    }
    else{
        // primary vertex was single-track and removed or non-existant at all, so don't have to remove it.
        _ip = trueIP;
    }
    streamlog_out(DEBUG3)<<"List of true vertices after merging and removing 1-track vertices and IP:"<<endl;
    for(auto iter=vtx2prongs.begin(); iter!=vtx2prongs.end(); ++iter){
        streamlog_out(DEBUG3)<<" Vertex #"<<std::distance(vtx2prongs.begin(), iter)<<"  ("<<( (iter->first) - trueIP ).r()<<")  MCs:  ";
        for(auto mc : iter->second) streamlog_out(DEBUG3)<<mc<<"    ";
        streamlog_out(DEBUG3)<<endl;
    }

    // 6) Create a map of secondary vertices found by LCFIPlus
    std::map< Vertex*, std::vector<MCParticle*>, CompareVertices > recovtx2prongs;
    LCCollection* vtxCol = event->getCollection(_vtxColName);
    streamlog_out(DEBUG3)<<"Number of reco vertices : "<<vtxCol->getNumberOfElements()<<endl;
    for(int i=0; i<vtxCol->getNumberOfElements(); ++i){
        Vertex* vertex = static_cast<Vertex*> (vtxCol->getElementAt(i));
        std::vector<ReconstructedParticle*> prongs = vertex->getAssociatedParticle()->getParticles();
        for(auto prong : prongs){
            MCParticle* mc = nullptr;
            if (_refitOpt){
                ReconstructedParticle* defaultPfo = getDefaultPfo(event, prong);
                mc = getMcMaxTrackWeight(defaultPfo, navRecoToMc);
            }
            else{
                mc = getMcMaxTrackWeight(prong, navRecoToMc);
            }
            if (mc != nullptr) recovtx2prongs[vertex].push_back(mc);
        }
    }


    LCCollection* vtxV0Col = event->getCollection(_vtxV0ColName);
    streamlog_out(DEBUG3)<<"Number of reco V0 vertices : "<<vtxV0Col->getNumberOfElements()<<endl;
    for(int i=0; i<vtxV0Col->getNumberOfElements(); ++i){
        Vertex* vertex = static_cast<Vertex*> (vtxV0Col->getElementAt(i));
        std::vector<ReconstructedParticle*> prongs = vertex->getAssociatedParticle()->getParticles();
        for(auto prong : prongs){
            MCParticle* mc = nullptr;
            if (_refitOpt){
                ReconstructedParticle* defaultPfo = getDefaultPfo(event, prong);
                mc = getMcMaxTrackWeight(defaultPfo, navRecoToMc);
            }
            else{
                mc = getMcMaxTrackWeight(prong, navRecoToMc);
            }
            if (mc != nullptr) recovtx2prongs[vertex].push_back(mc);
        }
    }
    streamlog_out(DEBUG3)<<"List of reco vertices:"<<endl;
    for(auto iter=recovtx2prongs.begin(); iter!=recovtx2prongs.end(); ++iter){
        XYZVectorF pos;
        pos.SetCoordinates( (iter->first)->getPosition() );
        int idx = std::distance(recovtx2prongs.begin(), iter);
        bool isV0 = ( (idx - vtxCol->getNumberOfElements() ) >= 0);
        streamlog_out(DEBUG3)<<" Vertex #"<<idx<<" isV0: "<<isV0<<"  ("<<( pos - trueIP ).r()<<")  MCs:  ";
        for(auto mc : iter->second) streamlog_out(DEBUG3)<<mc<<"    ";
        streamlog_out(DEBUG3)<<endl;
    }
    if (vtx2prongs.size() < recovtx2prongs.size()) streamlog_out(DEBUG8)<<"WARNING EVENT: "<< _nEvent<<std::endl;
    

    auto momComp = [](MCParticle* mc1, MCParticle* mc2){
        XYZVector mom1 (mc1->getMomentum()[0], mc1->getMomentum()[1], mc1->getMomentum()[2]);
        XYZVector mom2 (mc2->getMomentum()[0], mc2->getMomentum()[1], mc2->getMomentum()[2]);
        return mom1.r() < mom2.r();
    };

    // MAIN LOOP over all TRUE vertices to match to reco vertices, extract and write data
    for(auto trueIter=vtx2prongs.begin(); trueIter != vtx2prongs.end(); ++trueIter){
        streamlog_out(DEBUG4)<<"****** Analyzing vertex # "<<std::distance(vtx2prongs.begin(), trueIter)<<" ***********"<<endl;

        std::vector<MCParticle*> trueMcs = trueIter->second;
        // WRITE TRUE INFORMATION
        _vtxTrue = trueIter->first;
        streamlog_out(DEBUG4)<<"True position: "<<_vtxTrue.r()<<" mm"<<endl;
        streamlog_out(DEBUG4)<<"D to IP: "<<(_vtxTrue - _ip).r()<<" mm"<<endl;
        _nTrueTracks = trueMcs.size();
        streamlog_out(DEBUG4)<<"N true tracks: "<<_nTrueTracks<<endl;

        std::vector<MCParticle*> trueKaons;
        std::copy_if(trueMcs.begin(), trueMcs.end(), std::back_inserter(trueKaons), [](MCParticle* mc){return std::abs(mc->getPDG()) == 321;});
        _nTrueKaons = trueKaons.size();
        streamlog_out(DEBUG4)<<"N true kaons: "<<_nTrueKaons<<endl;
        _pTrueKaon = XYZVector();
        if (_nTrueKaons != 0){
            MCParticle* mc = *std::min_element( trueKaons.begin(), trueKaons.end(), momComp);
            _pTrueKaon = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom true kaon: "<<_pTrueKaon.r()<<" GeV"<<endl;

        std::vector<MCParticle*> trueProtons;
        std::copy_if(trueMcs.begin(), trueMcs.end(), std::back_inserter(trueProtons), [](MCParticle* mc){return std::abs(mc->getPDG()) == 2212;});
        _nTrueProtons = trueProtons.size();
        streamlog_out(DEBUG4)<<"N true protons: "<<_nTrueProtons<<endl;
        _pTrueProton = XYZVector();
        if (_nTrueProtons != 0){
            MCParticle* mc = *std::min_element( trueProtons.begin(), trueProtons.end(), momComp);
            _pTrueProton = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom true proton: "<<_pTrueProton.r()<<" GeV"<<endl;


        _isMatched = false;
        _isMissed = false;
        _isConfused = false;
        // there are no secondary reconstructed vertices in this event. Mark all as missed and write dummy numbers
        if ( recovtx2prongs.empty() ){
            streamlog_out(DEBUG4)<<"NO RECONSTRUCTED VERTICES: isMissed = True "<<endl;
            _isMissed = true;
            _vtxReco = XYZVector();
            _sigmaX = 0.;
            _sigmaY = 0.;
            _sigmaZ = 0.;
            _nRecoTracks = 0;
            _nRecoKaons = 0;
            _nRecoProtons = 0;
            _pRecoKaon = XYZVector();
            _pRecoProton = XYZVector();
            _nMissedTracks = _nTrueTracks;
            _nMissedKaons = _nTrueKaons;
            _nMissedProtons = _nTrueProtons;
            _pMissedKaon = _pTrueKaon;
            _pMissedProton = _pTrueProton;
            _nConfusedTracks = 0;
            _nConfusedKaons = 0;
            _nConfusedProtons = 0;
            _pConfusedKaon = XYZVector();
            _pConfusedProton = XYZVector();
            _tree->Fill();
            continue;
        }


        // 7) For each TRUE vertex calculate number of matched MCParticles to each RECONSTRUCTED vertex.
        std::map< Vertex*, int, CompareVertices > recovtx2nMatched;
        for(auto recoIter=recovtx2prongs.begin(); recoIter != recovtx2prongs.end(); ++recoIter){
            int nMatched = 0;
            std::vector<MCParticle*> recoMcs = recoIter->second;
            for (auto mc : recoMcs){
                if ( std::find(trueMcs.begin(), trueMcs.end(), mc) != trueMcs.end() ) nMatched++;
            }
            recovtx2nMatched[recoIter->first] = nMatched;
        }
        // 8) Reconstructed vertex corresponding to this TRUE vertex will be the one with a max number of matched MCParticles 
        auto matchedRecoIter = std::max_element( recovtx2nMatched.begin(), recovtx2nMatched.end(), [](const std::pair<Vertex*, int>& a, const std::pair<Vertex*, int>& b) {return a.second < b.second;} );

        int nMatchedTracks = matchedRecoIter->second;
        streamlog_out(DEBUG4)<<"N matched tracks: "<<nMatchedTracks<<endl;

        if ( nMatchedTracks == 0 ){
            //if among all reco vertices max number of matched tracks is 0 it means this true vertex is missed completely. Do the same in case of no reco vertices
            streamlog_out(DEBUG4)<<"NO MATCHED RECONSTRUCTED VERTICES: isMissed = True "<<endl;
            _isMissed = true;
            _vtxReco = XYZVector();
            _sigmaX = 0.;
            _sigmaY = 0.;
            _sigmaZ = 0.;
            _nRecoTracks = 0;
            _nRecoKaons = 0;
            _nRecoProtons = 0;
            _pRecoKaon = XYZVector();
            _pRecoProton = XYZVector();
            _nMissedTracks = _nTrueTracks;
            _nMissedKaons = _nTrueKaons;
            _nMissedProtons = _nTrueProtons;
            _pMissedKaon = _pTrueKaon;
            _pMissedProton = _pTrueProton;
            _nConfusedTracks = 0;
            _nConfusedKaons = 0;
            _nConfusedProtons = 0;
            _pConfusedKaon = XYZVector();
            _pConfusedProton = XYZVector();
            _tree->Fill();
            continue;
        }
        
        // If we are here, we have a reconstructed vertex with more than one matched track
        _nMissedTracks = _nTrueTracks - nMatchedTracks;
        streamlog_out(DEBUG4)<<"N missed tracks: "<<_nMissedTracks<<endl;
        Vertex* recoVertex = matchedRecoIter->first;
        _vtxReco.SetCoordinates( recoVertex->getPosition() );
        streamlog_out(DEBUG4)<<"Reconstructed position: "<<_vtxReco.r()<<" mm"<<endl;
        streamlog_out(DEBUG4)<<"Reco position to IP: "<<(_vtxReco - _ip).r()<<" mm"<<endl;
        _sigmaX = std::sqrt( recoVertex->getCovMatrix()[0] );
        _sigmaY = std::sqrt( recoVertex->getCovMatrix()[2] );
        _sigmaZ = std::sqrt( recoVertex->getCovMatrix()[5] );
        streamlog_out(DEBUG4)<<"SigmaX: "<<_sigmaX<<" mm"<<endl;
        streamlog_out(DEBUG4)<<"SigmaY: "<<_sigmaY<<" mm"<<endl;
        streamlog_out(DEBUG4)<<"SigmaZ: "<<_sigmaZ<<" mm"<<endl;
        std::vector<MCParticle*> recoMcs = recovtx2prongs[recoVertex];
        _nRecoTracks = recoMcs.size();
        streamlog_out(DEBUG4)<<"N reco tracks: "<<_nRecoTracks<<endl;
        _nConfusedTracks = _nRecoTracks - nMatchedTracks;
        streamlog_out(DEBUG4)<<"N confused tracks: "<<_nConfusedTracks<<endl;
        if (_nMissedTracks == 0 && _nConfusedTracks == 0) _isMatched = true; // perfectly matched w/o missess or confused tracks
        else _isConfused = true;
        streamlog_out(DEBUG4)<<"isMatched: "<<_isMatched<<endl;
        streamlog_out(DEBUG4)<<"isMissed: "<<_isMissed<<endl;
        streamlog_out(DEBUG4)<<"isConfused: "<<_isConfused<<endl;


        std::vector<MCParticle*> recoKaons;
        std::copy_if(recoMcs.begin(), recoMcs.end(), std::back_inserter(recoKaons), [](MCParticle* mc){return std::abs(mc->getPDG()) == 321;});
        _nRecoKaons = recoKaons.size();
        streamlog_out(DEBUG4)<<"N reco kaons: "<<_nRecoKaons<<endl;
        _pRecoKaon = XYZVector();
        if (_nRecoKaons != 0){
            MCParticle* mc = *std::min_element( recoKaons.begin(), recoKaons.end(), momComp);
            _pRecoKaon = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom reco kaon: "<<_pRecoKaon.r()<<" GeV"<<endl;

        std::vector<MCParticle*> recoProtons;
        std::copy_if(recoMcs.begin(), recoMcs.end(), std::back_inserter(recoProtons), [](MCParticle* mc){return std::abs(mc->getPDG()) == 2212;});
        _nRecoProtons = recoProtons.size();
        streamlog_out(DEBUG4)<<"N reco protons: "<<_nRecoProtons<<endl;
        _pRecoProton = XYZVector();
        if (_nRecoProtons != 0){
            MCParticle* mc = *std::min_element( recoProtons.begin(), recoProtons.end(), momComp);
            _pRecoProton = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom reco proton: "<<_pRecoProton.r()<<" GeV"<<endl;

        // FILL INFORMATION ABOUT MISSED PARTICLES
        std::vector<MCParticle*> missedKaons;
        for(auto trueKaon : trueKaons){
            if ( std::find(recoKaons.begin(), recoKaons.end(), trueKaon) == recoKaons.end() ) missedKaons.push_back(trueKaon);
        }
        _nMissedKaons = missedKaons.size();
        streamlog_out(DEBUG4)<<"N missed kaons: "<<_nMissedKaons<<endl;
        _pMissedKaon = XYZVector();
        if (_nMissedKaons != 0){
            MCParticle* mc = *std::min_element( missedKaons.begin(), missedKaons.end(), momComp);
            _pMissedKaon = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom missed kaon: "<<_pMissedKaon.r()<<" GeV"<<endl;

        std::vector<MCParticle*> missedProtons;
        for(auto trueProton : trueProtons){
            if ( std::find(recoProtons.begin(), recoProtons.end(), trueProton) == recoProtons.end() ) missedProtons.push_back(trueProton);
        }
        _nMissedProtons = missedProtons.size();
        streamlog_out(DEBUG4)<<"N missed protons: "<<_nMissedProtons<<endl;
        _pMissedProton = XYZVector();
        if (_nMissedProtons != 0){
            MCParticle* mc = *std::min_element( missedProtons.begin(), missedProtons.end(), momComp);
            _pMissedProton = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom missed proton: "<<_pMissedProton.r()<<" GeV"<<endl;

        // FILL INFORMATION ABOUT CONFUSED PARTICLES

        std::vector<MCParticle*> confusedKaons;
        for(auto recoKaon : recoKaons){
            if ( std::find(trueKaons.begin(), trueKaons.end(), recoKaon) == trueKaons.end() ) confusedKaons.push_back(recoKaon);
        }
        _nConfusedKaons = confusedKaons.size();
        streamlog_out(DEBUG4)<<"N confused kaons: "<<_nConfusedKaons<<endl;
        _pConfusedKaon = XYZVector();
        if (_nConfusedKaons != 0){
            MCParticle* mc = *std::min_element( confusedKaons.begin(), confusedKaons.end(), momComp);
            _pConfusedKaon = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom confused kaon: "<<_pConfusedKaon.r()<<" GeV"<<endl;

        std::vector<MCParticle*> confusedProtons;
        for(auto recoProton : recoProtons){
            if ( std::find(trueProtons.begin(), trueProtons.end(), recoProton) == trueProtons.end() ) confusedProtons.push_back(recoProton);
        }
        _nConfusedProtons = confusedProtons.size();
        streamlog_out(DEBUG4)<<"N confused protons: "<<_nConfusedProtons<<endl;
        _pConfusedProton = XYZVector();
        if (_nConfusedProtons != 0){
            MCParticle* mc = *std::min_element( confusedProtons.begin(), confusedProtons.end(), momComp);
            _pConfusedProton = XYZVector( mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2] );
        }
        streamlog_out(DEBUG4)<<"Mom confused proton: "<<_pConfusedProton.r()<<" GeV"<<endl;

        _tree->Fill();
    }

}

void VertexAnalysis::end(){
    _file->Write();
}

void VertexAnalysis::prepareRootTree(){
    std::string filename;
    if (_refitOpt) filename = "after_refit.root";
    else filename = "before_refit.root";
    _file.reset( new TFile(filename.c_str(), "RECREATE") );
    _tree.reset( new TTree("VertexAnalysis", "VertexAnalysis") );
    //true information
    _tree->Branch("ip_pos", &_ip);
    _tree->Branch("pos_true", &_vtxTrue);
    _tree->Branch("n_true_tracks", &_nTrueTracks);
    _tree->Branch("n_true_kaons", &_nTrueKaons);
    _tree->Branch("n_true_protons", &_nTrueProtons);
    _tree->Branch("p_true_kaon", &_pTrueKaon);
    _tree->Branch("p_true_proton", &_pTrueProton);
    _tree->Branch("pos_reco", &_vtxReco);
    _tree->Branch("sigma_x", &_sigmaX);
    _tree->Branch("sigma_y", &_sigmaY);
    _tree->Branch("sigma_z", &_sigmaZ);
    _tree->Branch("n_reco_tracks", &_nRecoTracks);
    _tree->Branch("n_reco_kaons", &_nRecoKaons);
    _tree->Branch("n_reco_protons", &_nRecoProtons);
    _tree->Branch("p_reco_kaon", &_pRecoKaon);
    _tree->Branch("p_reco_proton", &_pRecoProton);
    _tree->Branch("n_missed_tracks", &_nMissedTracks);
    _tree->Branch("n_missed_kaons", &_nMissedKaons);
    _tree->Branch("n_missed_protons", &_nMissedProtons);
    _tree->Branch("p_missed_kaon", &_pMissedKaon);
    _tree->Branch("p_missed_proton", &_pMissedProton);
    _tree->Branch("n_confused_tracks", &_nConfusedTracks);
    _tree->Branch("n_confused_kaons", &_nConfusedKaons);
    _tree->Branch("n_confused_protons", &_nConfusedProtons);
    _tree->Branch("p_confused_kaon", &_pConfusedKaon);
    _tree->Branch("p_confused_proton", &_pConfusedProton);
    _tree->Branch("is_matched", &_isMatched);
    _tree->Branch("is_missed", &_isMissed);
    _tree->Branch("is_confused", &_isConfused);

}


EVENT::MCParticle* VertexAnalysis::getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav){
    const vector<LCObject*>& mcs = nav.getRelatedToObjects(pfo);
    const vector<float>& weights = nav.getRelatedToWeights(pfo);
    
    streamlog_out(DEBUG1)<<"For pfo "<<pfo<<" weights:"<<endl;
    for(size_t i=0; i< mcs.size(); ++i){
        streamlog_out(DEBUG1)<<"MC: "<<mcs[i]<<"  track weight: "<< (int(weights[i])%10000)/1000.<<"  calo weight: "<<(int(weights[i])/10000)/1000.<<endl;
    }

    if (mcs.size() == 0) return nullptr;
    //get index of highest TRACK weight MC particle
    int i = std::max_element(weights.begin(), weights.end(), [](float a, float b){return (int(a)%10000)/1000. < (int(b)%10000)/1000.;}) - weights.begin();
    // there is no track if max TRACK weight is 0!
    if ( (int(weights[i])%10000)/1000. == 0 ) return nullptr;
    return static_cast<MCParticle*> ( mcs[i] );
}


ReconstructedParticle* VertexAnalysis::getDefaultPfo(EVENT::LCEvent* event, EVENT::ReconstructedParticle* selectedPFO){
    LCCollection* defaultPfoCol = event->getCollection("PandoraPFOs");
    LCCollection* newPfoCol = event->getCollection("updatedPandoraPFOs");
    for(int i=0; i < newPfoCol->getNumberOfElements(); ++i){
        ReconstructedParticle* newPfo = static_cast<ReconstructedParticle*> ( newPfoCol->getElementAt(i) );
        if (selectedPFO == newPfo) return static_cast<ReconstructedParticle*> ( defaultPfoCol->getElementAt(i) );
    }
    return nullptr;
}
