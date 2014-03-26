#ifndef MACRODEF_H
#define MACRODEF_H

// Definition of macro for period
#define period2011_BD "period2011_BD"
#define period2011_EH "period2011_EH"
#define period2011_IK "period2011_IK"
#define period2011_I "period2011_I"
#define period2011_J "period2011_J"
#define period2011_K "period2011_K"
#define period2011_LM "period2011_LM"

#define period2012_All "period2012_All"

// Z mass
#define pdgZMass  91187.6 
#define pdgMuMass 105.6583715
#define pdgElMass 0.510998928

namespace doAnalysis{
	enum{
		StdHZZllll,
		trigeff4l,
		Zllll
	};
}

namespace  MCCollection{
	enum {
		MC11c,
		MC11d,
		MC12a,
		MC12b,
		MC12c
	};
}
namespace cutFlow {
	enum {
		Total,
		DataPreselection,
		Preselection,
		Trigger,
		Trigger4Mu,
		Trigger4e
	};
}
namespace cutFlowCH {
	enum {
		Total,
		Trigger,
		Lepton,
		SFOS,
		Kinematics,
		TriggerMatch,
		Z1Mass,
		Z2Mass,
		DeltaR,
		TrackIso,
		CaloIso,
		D0Sig,
		Final
	};
}

namespace cutElFlow {
	enum {
		Total,
		DataPreselection,
		Preselection,
		Trigger,
		Author,
		Loose,
		Eta,
		Et,
		ObjectQuality,
		Z0,
		OverLapElEl,
		OverLapClElEl,
		OverLap
	};
}

namespace cutMuFlow {
	enum {
		Total,
		DataPreselection,
		Preselection,
		Trigger,
		Author,
		Pt,
		Eta,
		BLayer,
		Pix,
		SCT,
		Holes,
		TRT,
		D0,
		OverLap
	};
}

namespace cutJetsFlow {
	enum {
		Total,
		DataPreselection,
		Preselection,
		Trigger,
		Pt,
		Eta,
		Pileup,
		Clean,
		OverLap
	};
}

namespace flavor {
	enum {
		Muon,
		Electron,
		Jet
	};
}

namespace leptonType{
	enum {
		MuonStaco,
		MuonCalo,
		MuonStandAlone,
		ElectronGSF
	};
}
namespace jetsType{
	enum {
		AntiKt4TopoEM,
		AntiKt4TopoEMTruth,
		AntiKt4TopoEM_Fid,
		AntiKt4TopoEMTruth_Fid
	};
}

namespace quadType{
	enum {
		Mu4,
		El4,
		Mu2El2,
		El2Mu2
	};
}

namespace analysisType{
	enum {
		Mu4,
		El4,
		Mu2El2,
		El2Mu2
	};
}

namespace electronCollection{
	enum {
		LoosePlusPlus,
		MultiLepton,
		Likelihood
	};
}
namespace muonCollection{
	enum {
		Loose
	};
}

namespace productionChannel{
	enum {
		VBF,
		VHLep,
		VHHad,		
		ggF,
		VH		
	};
}

namespace sampleType{
	enum {
		ggF,
		VBF,
		WH,
		ZH,
		ttH,
		qqF,
		Background,
		ggF_ZpZp
	};
}

namespace streamContainer{
	enum {
		eGamma,
		Muon,
		Other
	};
}

namespace truthTypeQuad{
	enum {
		_4mu,
		_2mu2e,
		_4e,
		_4tau,
		_2tau2mu,
		_2tau2e,
		noStatus
	};
}

namespace leptonIDType{
	enum {
		// what we use
		mu_staco_cb,
		mu_staco_st,
		mu_staco_sa,
		// what we used to use
		mu_muid_cb,
		mu_muid_st,
		mu_muid_sa,
		// what we use
		mu_calomuon,
		// For future
		mu_muon_cb,
		mu_muon_st,
		mu_muon_sa,
		mu_muon_calomuon,
		// Electrons 2011
		el_loosepp_H4l,
		el_loosepp_H4l_Ep_comb,
		el_loosepp_H4l_relax,
		// Electrons 2012
		el_multilepton,
		el_multilepton_Ep_comb,
		el_multilepton_relax,
		el_likelihood_loose,
		el_likelihood_loose_Ep_comb,
		el_likelihood_loose_relax,
		unknown
	};
}

namespace MCGeneratorName{
	enum {
		Pythia,
		other
	};
}

namespace calibrationType{
	enum {
		stdCalib,
		stdCalibEp,
		MvaCalib,
		MvaCalibEp,
		noCalib
	};
}
namespace fsrType{
	enum {
		collFSRZ1mumu,
		farFSRZ1,
		farFSRZ2, 
		noFSR
	};
}

namespace dataCalibType{
	enum {
		y2011c,
		y2011d,
		y2012ab,
		y2012c
	};
}

namespace diLeptonType{
	enum {
		_2e, 
		_2mu
	};
}

namespace muonType{
	enum
	{
		CB, 
		MS,
		ID
	};
}

namespace doSys{
	enum
	{
		ZeeStatUp,
		ZeeStatDown, 
		ZeeSystUp,
		ZeeSystDown,
		ZeeAllUp,
		ZeeAllDown,
		PSUp,
		PSDown,
		S12Up,
		S12Down,
		MatIDUp,
		MatIDDown,
		MatCryoUp,
		MatCryoDown,
		MatCaloUp,
		MatCaloDown,
		LArCalibUp,
		LArCalibDown,
		LArUnconvCalibUp,
		LArUnconvCalibDown,
		LArElecUnconvUp,
		LArElecUnconvDown,
		LArElecCalibUp,
		LArElecCalibDown,
		GainUp,
		GainDown,
		G4Up,
		G4Down,
		MomentumUp,
		MomentumDown,
		ZSmearingUp,
		ZSmearingDown,
		SamplingTermUp,
		SamplingTermDown,
		MaterialIDUp,
		MaterialIDDown,
		MaterialCaloUp,
		MaterialCaloDown,
		MaterialGapUp,
		MaterialGapDown,
		MaterialCryoUp,
		MaterialCryoDown,
		PileUpUp,
		PileUpDown,
		Nom
	};
}

// For category
namespace VHLeptonType{
	enum {
		electronPlus,
		electronMinus,
		muonPlus,
		muonMinus,
		unknown
	};
}
namespace VHTruthType{
	enum {
		ZnoHiggs,
		WnoHiggs,
		Z1fromHiggs,
		Z2fromHiggs,
		unknown
	};
}
#endif

