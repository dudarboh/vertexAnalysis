std::map< XYZVector, vector<MCParticle*>, CompareVectors> candidates;
for( auto& pfo : pfos ){
    if (_refitOpt) pfo = getDefaultPfo(event, pfo);
    //track weights in relations are based on old tracks/before refit! Does it have an impact?
    MCParticle* mc = getMcMaxTrackWeight(pfo, navRecoToMc);
    XYZVector pos;
    pos.SetCoordinates( mc->getVertex() );
    //check if MCParticle is already there. There are PFOs that may link to the same MCParticle entering dublicate without checking!
    bool alreadyInList = std::find(candidates[pos].begin(), candidates[pos].end(), mc) != candidates[pos].end();
    if (not alreadyInList) candidates[pos].push_back(mc);
}

//choose true vertex position as the one with highest number of MCParticles (prongs)
auto compByMcSize = [](const std::pair<XYZVector, vector<MCParticle*> >& a, const std::pair<XYZVector, vector<MCParticle*> >& b){return a.second.size() < b.second.size();};
std::pair<XYZVector, vector<MCParticle*> > bestCandidate = ( *std::max_element(candidates.begin(), candidates.end(), compByMcSize) );
int nMaxProngs = bestCandidate.second.size();

auto hasMaxProngs = [&nMaxProngs](const std::pair<XYZVector, vector<MCParticle*> >& a){return a.second.size() == nMaxProngs;};
if (std::count_if(candidates.begin(), candidates.end(), hasMaxProngs) == 1){
    //only single max value -- easy case just take vertex at this index
    _vtxMcPos = bestCandidate.first;
}
else{
    //multiple true vertex position with the same amount of MCParticles/prongs/tracks. Take furthest from the IP
    double dToIp = -1.;
    for(auto const& candidate : candidates){
        XYZVector pos = candidate.first;
        int nProngs = candidate.second.size();
        if (nProngs == nMaxProngs && (pos - _ip).r() > dToIp ){
            dToIp = (pos - _ip).r();
            bestCandidate = candidate;
            _vtxMcPos = pos;           
        }
    }
}
int nTrueRecoProngs = nMaxProngs;
_nTracks = recoVertex->getAssociatedParticle()->getParticles().size();
_nConfusedTracks = _nTracks - nTrueRecoProngs;

//get meson type K/D/B or just other
MCParticle* trueProng = bestCandidate.second[0];
MCParticle* parent = nullptr;
int pdg;
if (trueProng->getParents().size() == 0 ){
    //weird no parents at all case
    _decayParentType = 0;
    _nMissedTracks = 0;
    _tree->Fill();
    continue;
}
//it really must have only 1 parent in the most of the cases..
parent = trueProng->getParents()[0];
XYZVector parentVertex;
parentVertex.SetCoordinates( parent->getVertex() );
pdg = parent->getPDG();
// do not count parents that come from the same vertex (instant decay). We are interested in parent who created 2ndary vertex thus lived some time
// If PDG == 92, it is already hadronization, vertex recosntructed from IP or most likely overlay
while ( parentVertex == _vtxMcPos && pdg != 92 && parent->getParents().size() != 0 ){
    parent = parent->getParents()[0];
    parentVertex.SetCoordinates( parent->getVertex() );
    pdg = parent->getPDG();
}

vector<int> strangeMesons = {130, 310, 311, 321, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 9000319, 9000329};
vector<int> charmMesons = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435};
vector<int> bottomMesons = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 10541, 543, 10543, 20543, 545};
if( std::find(strangeMesons.begin(), strangeMesons.end(), std::abs(pdg)) != strangeMesons.end() ){
    _decayParentType = 3;
}
else if ( std::find(charmMesons.begin(), charmMesons.end(), std::abs(pdg)) != charmMesons.end() ){
    _decayParentType = 4;
}
else if ( std::find(bottomMesons.begin(), bottomMesons.end(), std::abs(pdg)) != bottomMesons.end() ){
    _decayParentType = 5;
}
else{
    _decayParentType = 0;
}

//get number of missed tracks
int nTrueProngs = 0;
_hasKaon = false;
_hasProton = false;
vector<MCParticle*> decayChain;
fillDecayChainDown(parent, decayChain);
for(auto& daughter : decayChain){
    XYZVector daughterVertex;
    daughterVertex.SetCoordinates( daughter->getVertex() );
    if (daughterVertex == _vtxMcPos && daughter->getCharge() != 0 && daughter->getGeneratorStatus() <= 1){
        nTrueProngs++;
        if ( std::abs(daughter->getPDG()) == 321 ) _hasKaon = true;
        else if (std::abs(daughter->getPDG()) == 2212) _hasProton = true;
    }
}

_nMissedTracks = nTrueProngs - nTrueRecoProngs;

cout<<"****** Analyzing vertex # "<<i<<" ***********"<<endl;
cout<<"Reconstructed position: "<<_vtxRecoPos.r()<<" mm"<<endl;
cout<<"True position: "<<_vtxMcPos.r()<<" mm"<<endl;
cout<<"D to IP: "<<(_vtxMcPos - _ip).r()<<" mm"<<endl;
cout<<"total reconstructed prongs: "<<_nTracks<<endl;
cout<<"total TRUE reconstructed prongs: "<<nTrueRecoProngs<<endl;
cout<<"Confusing tracks: "<<_nConfusedTracks<<endl;
cout<<"Decay parent: "<<_decayParentType<<endl;
cout<<"total TRUE prongs: "<<nTrueProngs<<endl;
cout<<"Missed tracks: "<<_nMissedTracks<<endl;
cout<<"Has Kaon: "<<_hasKaon<<endl;
cout<<"Has Proton: "<<_hasProton<<endl;














void VertexAnalysis::fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& decayChain){
    // stop iterating down at particles not created by generator
    // if(mc->getGeneratorStatus() == 0) return;
    decayChain.push_back(mc);
    const vector<MCParticle*> daughters = mc->getDaughters();
    for(auto daughter : daughters){
        bool foundInTheList = std::find(decayChain.begin(), decayChain.end(), daughter) != decayChain.end();
        if ( !foundInTheList ) fillDecayChainDown(daughter, decayChain);
    }
}



















































// get all secondary vertices:
std::map<Vertex*, bool> vertices;
for(int i=0; i<vtxCol->getNumberOfElements(); ++i){
    Vertex* vertex = static_cast<Vertex*> (vtxCol->getElementAt(i));
    vertices[vertex] = false;
}
for(int i=0; i<vtxV0Col->getNumberOfElements(); ++i){
    Vertex* vertex = static_cast<Vertex*> (vtxV0Col->getElementAt(i));
    vertices[vertex] = true;
}

int nBuildUpVertices = std::count_if(vertices.begin(), vertices.end(), [](std::pair<Vertex*, bool> vtx) {return vtx.second == 0;});
int nBuildUpV0Vertices = std::count_if(vertices.begin(), vertices.end(), [](std::pair<Vertex*, bool> vtx) {return vtx.second == 1;});
std::cout<<"Analyzing "<<nBuildUpVertices<<" buildUp vertices and "<<nBuildUpV0Vertices<<" V0 vertices"<<std::endl;

//Loop over build up vertices and extract all the information
for(auto vertex : vertices){
    Vertex* recoVertex = vertex.first;
    _isV0 = vertex.second;
    cout<<"IP: "<<_ip.r()<<" mm"<<endl;
    cout<<"isV0: "<<_isV0<<" mm"<<endl;
    std::cout<<std::endl;
    _vtxReco.SetCoordinates( recoVertex->getPosition() );
    _sigmaX = std::sqrt( recoVertex->getCovMatrix()[0] );
    _sigmaY = std::sqrt( recoVertex->getCovMatrix()[2] );
    _sigmaZ = std::sqrt( recoVertex->getCovMatrix()[5] );
    _nRecoTracks = recoVertex->getAssociatedParticle()->getParticles().size();
    cout<<"Reco position to IP: "<<(_vtxReco - _ip).r()<<" mm"<<endl;
    cout<<"SigmaX: "<<_sigmaX<<" mm"<<endl;
    cout<<"SigmaY: "<<_sigmaY<<" mm"<<endl;
    cout<<"SigmaZ: "<<_sigmaZ<<" mm"<<endl;
    cout<<"N reco tracks: "<<_nRecoTracks<<endl;

    vector<ReconstructedParticle*> recoTracks = recoVertex->getAssociatedParticle()->getParticles();
    vector<MCParticle*> recoMCs;
    vector<MCParticle*> recoKaons;
    vector<MCParticle*> recoProtons;
    vector<XYZVector> momRecoKaons;
    vector<XYZVector> momRecoProtons;
    vector<XYZVector> posRecoMCs;
    for( size_t j=0; j<recoTracks.size(); ++j ){
        ReconstructedParticle* pfo = recoTracks[j];
        MCParticle* mc = nullptr;
        //get MC particle. track weights in relations are based on old tracks/before refit! Does it have an impact?
        if (_refitOpt){
            ReconstructedParticle* defaultPfo = getDefaultPfo(event, pfo);
            mc = getMcMaxTrackWeight(defaultPfo, navRecoToMc);
        }
        else{
            mc = getMcMaxTrackWeight(pfo, navRecoToMc);
        }
        if (mc == nullptr) continue;
        recoMCs.push_back(mc);
        posRecoMCs.push_back( XYZVector(mc->getVertex()[0], mc->getVertex()[1], mc->getVertex()[2]) );
        if ( std::abs(mc->getPDG()) == 321 ){
            recoKaons.push_back(mc);
            // momRecoKaons.push_back( XYZVector(pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2]) );
            momRecoKaons.push_back( XYZVector(mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2]) );
        }
        if ( std::abs(mc->getPDG()) == 2212 ){
            recoProtons.push_back(mc);
            // momRecoProtons.push_back( XYZVector(pfo->getMomentum()[0], pfo->getMomentum()[1], pfo->getMomentum()[2]) );
            momRecoProtons.push_back( XYZVector(mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2]) );
        }
    }

    _nRecoKaons = recoKaons.size();
    _nRecoProtons = recoProtons.size();
    cout<<"N reco kaons: "<<_nRecoKaons<<endl;
    cout<<"N reco protons: "<<_nRecoProtons<<endl;
    auto compVec = [](const XYZVector& a, const XYZVector& b) { return a.r() < b.r(); };
    //save only lowest momentum particle
    _pRecoKaon = XYZVector();
    _pRecoProton = XYZVector();
    if (_nRecoKaons != 0) _pRecoKaon = *std::min_element(momRecoKaons.begin(), momRecoKaons.end(), compVec);
    if (_nRecoProtons != 0) _pRecoProton = *std::min_element(momRecoProtons.begin(), momRecoProtons.end(), compVec);
    cout<<"Mom reco kaon: "<<_pRecoKaon.r()<<" GeV"<<endl;
    cout<<"Mom reco proton: "<<_pRecoProton.r()<<" GeV"<<endl;
    std::cout<<std::endl;


    //Begin search for a true vertex
    // Start with selecting unique position candidates. And count number of occurances of each

    // Group position candidates based on number of MC particles they have
    std::map<XYZVector, int, CompareVectors> truePosCandidates;
    for(auto& pos : posRecoMCs) truePosCandidates[pos]++;

    auto compMapFirst = [](std::pair<XYZVector, int> a, std::pair<XYZVector, int> b){return (a.first).r() < (b.first).r();};
    auto compMapSecond = [](std::pair<XYZVector, int> a, std::pair<XYZVector, int> b){return a.second < b.second;};
    //Find the highest number of particles coming from the vertex
    int nMaxProngs = (*std::max_element(truePosCandidates.begin(), truePosCandidates.end(), compMapSecond )).second;

    // retain only candidates which have the highest number of MCParticles prongs
    for(auto it = truePosCandidates.begin(); it != truePosCandidates.end(); ) {
        if(it->second != nMaxProngs) it = truePosCandidates.erase(it);
        else ++it;
    }

    // In most cases it is only one left. But if not, pick the furthest from the IP
    _vtxTrue = (*std::max_element( truePosCandidates.begin(), truePosCandidates.end(), compMapFirst)).first;
    cout<<"True position to IP: "<<(_vtxTrue - _ip).r()<<" mm"<<endl;

    vector<MCParticle*> trueMCs;
    vector<MCParticle*> trueKaons;
    vector<MCParticle*> trueProtons;
    vector<XYZVector> momTrueKaons;
    vector<XYZVector> momTrueProtons;

    //Get true number of tracks (even those not reconstructed)
    for(int j=0; j<mcCol->getNumberOfElements(); ++j){
        MCParticle* mc = static_cast<MCParticle*> ( mcCol->getElementAt(j) );
        XYZVector pos(mc->getVertex()[0], mc->getVertex()[1], mc->getVertex()[2]);
        if (mc->getGeneratorStatus() == 1 && pos == _vtxTrue && mc->getCharge() != 0){
            trueMCs.push_back(mc);
            if (std::abs( mc->getPDG() ) == 321 ){
                trueKaons.push_back(mc);
                momTrueKaons.push_back( XYZVector(mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2]) );
            }
            if (std::abs( mc->getPDG() ) == 2212 ){
                trueProtons.push_back(mc);
                momTrueProtons.push_back( XYZVector(mc->getMomentum()[0], mc->getMomentum()[1], mc->getMomentum()[2]) );
            }
        }
    }

    _nTrueTracks = trueMCs.size();
    _nTrueKaons = trueKaons.size();
    _nTrueProtons = trueProtons.size();
    cout<<"N true tracks: "<<_nTrueTracks<<endl;
    cout<<"N true kaons: "<<_nTrueKaons<<endl;
    cout<<"N true protons: "<<_nTrueProtons<<endl;

    _pTrueKaon = XYZVector();
    _pTrueProton = XYZVector();
    if (_nTrueKaons != 0) _pTrueKaon = *std::min_element(momTrueKaons.begin(), momTrueKaons.end(), compVec);
    if (_nTrueProtons != 0) _pTrueProton = *std::min_element(momTrueProtons.begin(), momTrueProtons.end(), compVec);
    cout<<"Mom true kaon: "<<_pTrueKaon.r()<<" GeV"<<endl;
    cout<<"Mom true proton: "<<_pTrueProton.r()<<" GeV"<<endl;
    std::cout<<std::endl;

    vector<MCParticle*> missedMCs = trueMCs;
    vector<MCParticle*> missedKaons = trueKaons;
    vector<MCParticle*> missedProtons = trueProtons;
    vector<XYZVector> momMissedKaons;
    vector<XYZVector> momMissedProtons;

    for ( size_t j=0; j<recoMCs.size(); ++j){
        MCParticle* recoMC = recoMCs[j];
        if ( std::find(missedMCs.begin(), missedMCs.end(), recoMC) != missedMCs.end() ){
            auto end = std::remove(missedMCs.begin(), missedMCs.end(), recoMC);
            missedMCs.erase(end, missedMCs.end());
        }
        if ( std::find(missedKaons.begin(), missedKaons.end(), recoMC) != missedKaons.end() ){
            auto end = std::remove(missedKaons.begin(), missedKaons.end(), recoMC);
            missedKaons.erase(end, missedKaons.end());
        }
        if ( std::find(missedProtons.begin(), missedProtons.end(), recoMC) != missedProtons.end() ){
            auto end = std::remove(missedProtons.begin(), missedProtons.end(), recoMC);
            missedProtons.erase(end, missedProtons.end());
        }
    }
    for(auto kaon : missedKaons){
        momMissedKaons.push_back( XYZVector(kaon->getMomentum()[0], kaon->getMomentum()[1], kaon->getMomentum()[2]) );
    }
    for(auto proton : missedProtons){
        momMissedProtons.push_back( XYZVector(proton->getMomentum()[0], proton->getMomentum()[1], proton->getMomentum()[2]) );
    }

    _nMissedTracks = missedMCs.size();
    _nMissedKaons = missedKaons.size();
    _nMissedProtons = missedProtons.size();
    _pMissedKaon = XYZVector();
    _pMissedProton = XYZVector();
    cout<<"N missed tracks: "<<_nMissedTracks<<endl;
    cout<<"N missed kaons: "<<_nMissedKaons<<endl;
    cout<<"N missed protons: "<<_nMissedProtons<<endl;
    if (_nMissedKaons != 0) _pMissedKaon = *std::min_element(momMissedKaons.begin(), momMissedKaons.end(), compVec);
    if (_nMissedProtons != 0) _pMissedProton = *std::min_element(momMissedProtons.begin(), momMissedProtons.end(), compVec);
    cout<<"Mom missed kaon: "<<_pMissedKaon.r()<<" GeV"<<endl;
    cout<<"Mom missed proton: "<<_pMissedProton.r()<<" GeV"<<endl;
    std::cout<<std::endl;


    vector<MCParticle*> missedDetectedMCs = missedMCs;
    vector<MCParticle*> missedDetectedKaons = missedKaons;
    vector<MCParticle*> missedDetectedProtons = missedProtons;
    vector<XYZVector> momMissedDetectedKaons;
    vector<XYZVector> momMissedDetectedProtons;

    for ( size_t j=0; j<missedMCs.size(); ++j){
        MCParticle* mc = missedMCs[j];
        if ( link.find(mc) == link.end() ){
            auto end = std::remove(missedDetectedMCs.begin(), missedDetectedMCs.end(), mc);
            missedDetectedMCs.erase(end, missedDetectedMCs.end());
        }
    }

    for ( size_t j=0; j<missedKaons.size(); ++j){
        MCParticle* mc = missedKaons[j];
        if ( link.find(mc) == link.end() ){
            auto end = std::remove(missedDetectedKaons.begin(), missedDetectedKaons.end(), mc);
            missedDetectedKaons.erase(end, missedDetectedKaons.end());
        }
    }

    for ( size_t j=0; j<missedProtons.size(); ++j){
        MCParticle* mc = missedProtons[j];
        if ( link.find(mc) == link.end() ){
            auto end = std::remove(missedDetectedProtons.begin(), missedDetectedProtons.end(), mc);
            missedDetectedProtons.erase(end, missedDetectedProtons.end());
        }
    }

    for(auto kaon : missedDetectedKaons){
        momMissedDetectedKaons.push_back( XYZVector(kaon->getMomentum()[0], kaon->getMomentum()[1], kaon->getMomentum()[2]) );
    }
    for(auto proton : missedDetectedProtons){
        momMissedDetectedProtons.push_back( XYZVector(proton->getMomentum()[0], proton->getMomentum()[1], proton->getMomentum()[2]) );
    }

    _nMissedDetectedTracks = missedDetectedMCs.size();
    _nMissedDetectedKaons = missedDetectedKaons.size();
    _nMissedDetectedProtons = missedDetectedProtons.size();
    _pMissedDetectedKaon = XYZVector();
    _pMissedDetectedProton = XYZVector();
    cout<<"N missed Detected tracks: "<<_nMissedDetectedTracks<<endl;
    cout<<"N missed Detected kaons: "<<_nMissedDetectedKaons<<endl;
    cout<<"N missed Detected protons: "<<_nMissedDetectedProtons<<endl;
    if (_nMissedDetectedKaons != 0) _pMissedDetectedKaon = *std::min_element(momMissedDetectedKaons.begin(), momMissedDetectedKaons.end(), compVec);
    if (_nMissedDetectedProtons != 0) _pMissedDetectedProton = *std::min_element(momMissedDetectedProtons.begin(), momMissedDetectedProtons.end(), compVec);
    cout<<"Mom missed Detected kaon: "<<_pMissedDetectedKaon.r()<<" GeV"<<endl;
    cout<<"Mom missed Detected proton: "<<_pMissedDetectedProton.r()<<" GeV"<<endl;
    std::cout<<std::endl;



    vector<MCParticle*> confusedMCs = recoMCs;
    vector<MCParticle*> confusedKaons = recoKaons;
    vector<MCParticle*> confusedProtons = recoProtons;
    vector<XYZVector> momConfusedKaons;
    vector<XYZVector> momConfusedProtons;

    for ( size_t j=0; j<trueMCs.size(); ++j){
        MCParticle* trueMC = trueMCs[j];
        if ( std::find(confusedMCs.begin(), confusedMCs.end(), trueMC) != confusedMCs.end() ){
            auto end = std::remove(confusedMCs.begin(), confusedMCs.end(), trueMC);
            confusedMCs.erase(end, confusedMCs.end());
        }
        if ( std::find(confusedKaons.begin(), confusedKaons.end(), trueMC) != confusedKaons.end() ){
            auto end = std::remove(confusedKaons.begin(), confusedKaons.end(), trueMC);
            confusedKaons.erase(end, confusedKaons.end());
        }
        if ( std::find(confusedProtons.begin(), confusedProtons.end(), trueMC) != confusedProtons.end() ){
            auto end = std::remove(confusedProtons.begin(), confusedProtons.end(), trueMC);
            confusedProtons.erase(end, confusedProtons.end());
        }
    }
    XYZVector confusedVtxPos = XYZVector();
    for(auto kaon : confusedKaons){
        momConfusedKaons.push_back( XYZVector(kaon->getMomentum()[0], kaon->getMomentum()[1], kaon->getMomentum()[2]) );
    }

    for(auto proton : confusedProtons){
        momConfusedProtons.push_back( XYZVector(proton->getMomentum()[0], proton->getMomentum()[1], proton->getMomentum()[2]) );
    }

    _nConfusedTracks = confusedMCs.size();
    _nConfusedKaons = confusedKaons.size();
    _nConfusedProtons = confusedProtons.size();
    _pConfusedKaon = XYZVector();
    _pConfusedProton = XYZVector();
    cout<<"N confused tracks: "<<_nConfusedTracks<<endl;
    cout<<"N confused kaons: "<<_nConfusedKaons<<endl;
    cout<<"N confused protons: "<<_nConfusedProtons<<endl;

    if (_nConfusedKaons != 0) _pConfusedKaon = *std::min_element(momConfusedKaons.begin(), momConfusedKaons.end(), compVec);
    if (_nConfusedProtons != 0) _pConfusedProton = *std::min_element(momConfusedProtons.begin(), momConfusedProtons.end(), compVec);
    cout<<"Mom confused kaon: "<<_pConfusedKaon.r()<<" GeV"<<endl;
    cout<<"Mom confused proton: "<<_pConfusedProton.r()<<" GeV"<<endl;
    _dToTrue = 0.;
