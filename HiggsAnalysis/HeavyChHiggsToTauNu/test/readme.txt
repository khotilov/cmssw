10.10.2007/S.Lehti CMSSW_1_6_5	jet calibration data needed:
				cvs co -r jet_corrections_16X JetMETCorrections/MCJet/data

29.10.2007/S.Lehti CMSSW_1_6_7	Old SimpleTree version tagged with name "SimpleTree". 
				Contains the analysis, event filter, and storing 
				data using a simple root tree.

29.10.2007/S.Lehti CMSSW_1_6_7	Analysis changed to pTDR ORCA style using MyEvent 
				class for storing information in the root files. 
				In cvs head.

1.11.2007/S.Lehti  CMSSW_1_6_7	PFlow tags needed, co version V02-06-00, and compile.
				List of all cvs tags to be checked out:
				cvs co -r jet_corrections_16X JetMETCorrections/MCJet/data
				cvs co -r V02-06-00 DataFormats/ParticleFlowReco
				cvs co -r V02-06-00 DataFormats/ParticleFlowCandidate
				cvs co -r V02-06-00 RecoParticleFlow

5.11.2007/S.Lehti  CMSSW_1_6_7	Type1MET added in analysis, MET storing changed
				to contain rawMET and muon correction separately

8.11.2007/S.Lehti  CMSSW_1_6_7  Tau changed to use CaloTau and PFTau classes. Old code using
				IsolatedTauTagInfo tagged with name "IsolatedTauTagInfo". 
				New code in cvs head. 
				List of all cvs tags to be checked out:
                                cvs co -r V02-06-00 DataFormats/ParticleFlowReco
                                cvs co -r V02-06-00 DataFormats/ParticleFlowCandidate
                                cvs co -r V02-06-00 RecoParticleFlow
				cvs co -r V00-00-17 DataFormats/TauReco
				cvs co -r V00-00-37 RecoTauTag/RecoTau 
				cvs co -r V00-00-06 RecoTauTag/TauTagTools
				cvs co -r jet_corrections_16X JetMETCorrections/MCJet/data

22.11.2007/S.Lehti CMSSW_1_6_7  Added particleType in MyTrack class for storing
				PF ParticleType: X=0,h=1,e=2,mu=3,gamma=4,h0=5


	How to compile:
	- First make event dictionaries for the MyEvent class
		cd HiggsAnalysis/HeavyChHiggsToTauNu/src
		(chmod +x dictionary.csh if needed)
		dictionary.csh

		Rerun dictionary.csh every time after changing anything in MyEvent class
	-compile
		cd HiggsAnalysis/HeavyChHiggsToTauNu
		scramv1 b

	      * if you have taken any tags, then the compilation should be
                done in  CMSSW_1_6_X/src directory instead of 
                HiggsAnalysis/HeavyChHiggsToTauNu, so that everything gets
                compiled

	-run
	Take jet calibration data for 1_5_X samples:
	cvs co -r jet_corrections_16X JetMETCorrections/MCJet/data

	Example for reading and analyzing the produced root files can be found in:
	/afs/cern.ch/user/s/slehti/public/CMSSW/analysis_src

