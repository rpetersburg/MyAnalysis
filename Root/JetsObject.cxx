#include <stdlib.h>
#include <string>
#include "MyAnalysis/JetsObject.h"
#include <iostream>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////
// Just The constructor
JetsObject::JetsObject(Int_t year)
{
	dataYear = year;
	eventNumber = 0;
}

// Fill the intial vector with muon objects, Sets the type as well
void JetsObject::FillJets(D3PDReader::JetD3PDObject * jets_branch, Int_t type, Bool_t isMC, Int_t teventNumber)
{
	eventNumber = teventNumber;
	for(Int_t i = 0; i < jets_branch->n(); i++)
	{
		ChargedLepton *temp = new ChargedLepton (&((*jets_branch)[i]), type, i, isMC);
		jetsInitEvent.push_back(temp);
	}
}
// Fill the intial vector with truth jets  objects, Sets the type as well
void JetsObject::FillJetsTruth(D3PDReader::JetD3PDObject * jets_branch, Int_t type, Bool_t isMC)
{
	if(!isMC) return;

	for(Int_t i = 0; i < jets_branch->n(); i++)
	{
		ChargedLepton *temp = new ChargedLepton (&((*jets_branch)[i]), type, i, isMC);
		jetsInitEvent.push_back(temp);
	}
}

// Perfroms the Jets cut by calling the right funtion for the type
Bool_t JetsObject::JetsCut(Int_t *cutJetsPass, Int_t run)
{
	jetsBfOverlap.clear();
	for(vector<ChargedLepton *>::iterator itr = jetsInitEvent.begin();
			itr != jetsInitEvent.end(); ++itr)
	{
		ChargedLepton *jets_lep_i = *itr;
		D3PDReader::JetD3PDObjectElement *jets_i = jets_lep_i->GetJets();
		// HSG2 Jets
		if(jets_lep_i->type == jetsType::AntiKt4TopoEM)
		{
			if(!CutAntiKt4TopoEM(jets_i, cutJetsPass, run)) continue;
			else {
				jetsBfOverlap.push_back(jets_lep_i); 
			}
		}
		//true hets
		else if(jets_lep_i->type == jetsType::AntiKt4TopoEMTruth)
		{
			if(!CutAntiKt4TopoEMTruth(jets_i)) continue;
			else {
				jetsBfOverlap.push_back(jets_lep_i); 
			}
		}
		// Jets for Fiducial HSG1 definition
		else if(jets_lep_i->type == jetsType::AntiKt4TopoEM_Fid)
		{
			if(!CutAntiKt4TopoEM_Fid(jets_i, run)) continue;
			else {
				jetsBfOverlap.push_back(jets_lep_i); 
			}
		}
		//true hets
		else if(jets_lep_i->type == jetsType::AntiKt4TopoEMTruth_Fid)
		{
			if(!CutAntiKt4TopoEMTruth_Fid(jets_i)) continue;
			else {
				jetsBfOverlap.push_back(jets_lep_i); 
			}
		}


		
	}
	if(jetsBfOverlap.size() > 0) return true;
	else return false;
}
// HSG2 Jets
Bool_t JetsObject::CutAntiKt4TopoEM(D3PDReader::JetD3PDObjectElement *jets_i, Int_t *cutJetsPass, Int_t run)
{
	Double_t eta = jets_i->emscale_eta();
	Double_t pT = jets_i->pt();

	// Pt Cut
	if(pT> 25000)
	{cutJetsPass[cutJetsFlow::Pt]++;}
	else return false;
	
	// Eta Cut
	if((pT > 25000 && fabs(eta) < 2.4) || (pT > 30000 && fabs(eta) > 2.4 && fabs(eta) < 4.5))
	{cutJetsPass[cutJetsFlow::Eta]++;
	//cout << "event number " << eventNumber << "\tjet Pt\t"  <<  pT << endl;
	}
	else return false;
	
	// Pileup Removal
	Double_t jvfCut = 0;
	if(dataYear == 2011) {jvfCut = 0.75;}
	else if(dataYear == 2012) {jvfCut = 0.5;}

	if((pT < 50*1000 && fabs(eta) < 2.4))
	{	
		if(fabs(jets_i->jvtxf()) > jvfCut){}
		else return false;
	}
	cutJetsPass[cutJetsFlow::Pileup]++;
	
	// Jets Cleaning
	if(jets_i->isBadLooseMinus() == 0){}
	//{cutJetsPass[cutJetsFlow::Clean]++;}
	else return false;

	Double_t j_fmax = jets_i->fracSamplingMax();
 	Double_t j_smax = jets_i->SamplingMax();
  	Double_t j_eta  = jets_i->eta();
  	Double_t j_phi  = jets_i->phi();
  	//Int_t    run    =  TreeObject->RunNumber;                                                                                                                                

  	Bool_t _etaphi28 = false;
  	if( j_eta>-0.2 && j_eta<-0.1 && j_phi>2.65 && j_phi< 2.75 ) _etaphi28 = true;
  	Bool_t affected = (run==202660||run==202668||run==202712||run==202740||run==202965||run==202987||run==202991||run==203027) ? true : false;
  	if ( j_fmax>0.6 && j_smax==13 && _etaphi28 && affected ){return false;}

	cutJetsPass[cutJetsFlow::Clean]++;

	return true;
}
// HSG2 Jets
Bool_t JetsObject::CutAntiKt4TopoEMTruth(D3PDReader::JetD3PDObjectElement *jets_i)
{
	Double_t eta = jets_i->eta();
	Double_t pT = jets_i->pt();

	// Pt Cut
	if(pT> 25000){}
	else return false;
	
	// Eta Cut
	if((pT > 25000 && fabs(eta) < 2.4) || (pT > 30000 && fabs(eta) > 2.4 && fabs(eta) < 4.5)){}
	else return false;
	

	return true;

}

// HSG1_Fid Jets
Bool_t JetsObject::CutAntiKt4TopoEM_Fid(D3PDReader::JetD3PDObjectElement *jets_i, Int_t run)
{
	Double_t eta = jets_i->emscale_eta();
	Double_t pT = jets_i->pt();
	Double_t phi = jets_i->phi();
	Double_t M = jets_i->m();
	
	TLorentzVector	m_momentum;
	m_momentum.SetPtEtaPhiM(pT, jets_i->eta(), phi, M);


	//cout<<"Fid Reco Jet eta: "<<eta<<" pT: "<<pT<<" rapidity: "<<m_momentum.Rapidity()<<endl;
	// Pt Cut
	if(pT> 30*1000){}
	else return false;
	
	// Eta Cut
	if( fabs(m_momentum.Rapidity()) < 4.4 ){}
	else return false;
	
	// Pileup Removal
	Double_t jvfCut = 0.25;

	if(fabs(eta) < 2.4 && pT < 50 * 1000)
	{	
		if(fabs(jets_i->jvtxf()) > jvfCut){}
		else return false;
	}
	
	// Jets Cleaning
	if(jets_i->isBadLooseMinus() == 0){}
	else return false;

	Double_t j_fmax = jets_i->fracSamplingMax();
 	Double_t j_smax = jets_i->SamplingMax();
  	Double_t j_eta  = jets_i->eta();
  	Double_t j_phi  = jets_i->phi();
  	//Int_t    run    =  TreeObject->RunNumber;                                                                                                                                

  	Bool_t _etaphi28 = false;
  	if( j_eta>-0.2 && j_eta<-0.1 && j_phi>2.65 && j_phi< 2.75 ) _etaphi28 = true;
  	Bool_t affected = (run==202660||run==202668||run==202712||run==202740||run==202965||run==202987||run==202991||run==203027) ? true : false;
  	if ( j_fmax>0.6 && j_smax==13 && _etaphi28 && affected ){return false;}

	//cout<<"Jet passed the cut"<<endl;
	return true;
}
// HSG2 Jets
Bool_t JetsObject::CutAntiKt4TopoEMTruth_Fid(D3PDReader::JetD3PDObjectElement *jets_i)
{
	Double_t eta = jets_i->eta();
	Double_t pT = jets_i->pt();
	//cout<<"Fid true Jet eta: "<<eta<<" pT: "<<pT<<endl;
	Double_t phi = jets_i->phi();
	Double_t M = jets_i->m();
	
	TLorentzVector	m_momentum;
	m_momentum.SetPtEtaPhiM(pT, eta, phi, M);
	// Pt Cut
	if(pT> 30*1000){}
	else return false;
	
	// Eta Cut
	if( fabs(m_momentum.Rapidity()) < 4.4 ){}
	else return false;
	
	//cout<<"Jet passed the cut"<<endl;

	return true;

}


// Clear Vars
void JetsObject::clearVars()
{
	while(!jetsInitEvent.empty()) delete jetsInitEvent.back(), jetsInitEvent.pop_back();
	jetsInitEvent.clear(); 
	jetsAuthor.clear();				
	jetsBfOverlap.clear();		
	jetsOverlapGoodEvent.clear(); 
	jetsEvent.clear();

}

