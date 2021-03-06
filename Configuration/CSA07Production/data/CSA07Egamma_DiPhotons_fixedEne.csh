cat > CSA07Egamma_DiPhotons_Ene$1_GEN_SIM.cfg <<EOF
process Sim  = {

   untracked PSet maxEvents = {untracked int32 input = 25000}

   untracked PSet configurationMetadata = {
           untracked string version = "\$Revision\$"
           untracked string name = "\$Source\$"
           untracked string annotation = "back-to-back photons; fixed energy=$1 GeV"
   }


   include "FWCore/MessageService/data/MessageLogger.cfi"
   replace MessageLogger.cout.threshold = "ERROR"
   replace MessageLogger.cerr.default.limit = 10

   service = RandomNumberGeneratorService
  {
      untracked uint32 sourceSeed = 123456789
      PSet moduleSeeds =
      {
         untracked uint32 VtxSmeared = 98765432
         untracked uint32 g4SimHits = 9876
      }

  }

  include "SimGeneral/HepPDTESSource/data/pythiapdt.cfi"

  source = FlatRandomEGunSource 
  { 	
    untracked PSet PGunParameters =
    {
      untracked vint32 PartID = { 22 }
      untracked double MinEta = -2.7
      untracked double MaxEta =  2.7
      untracked double MinPhi = -3.14159265358979323846 # in radians
      untracked double MaxPhi =  3.14159265358979323846
      untracked double MinE  = $1.
      untracked double MaxE  = $1.
    }
    untracked bool AddAntiParticle = true # back-to-back particles
    untracked int32 Verbosity = 0 # set to 1 (or greater)  for printouts
  }
   
  include "Configuration/StandardSequences/data/Simulation.cff"
  include "Configuration/StandardSequences/data/VtxSmearedGauss.cff"
  include "Configuration/StandardSequences/data/MixingNoPileUp.cff"

  path p1 = {psim} 

  include "Configuration/EventContent/data/EventContent.cff"

   module GEN-SIM = PoolOutputModule
   {
        using FEVTSIMEventContent
        untracked string fileName = "file:DiPhotons_Ene$1.root"
	untracked PSet dataset = {
	    untracked string dataTier = "GEN-SIM"
      	}
   }
   endpath outpath = {GEN-SIM}
   schedule = {p1,outpath}
}
EOF
