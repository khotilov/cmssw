import FWCore.ParameterSet.Config as cms

source = cms.Source("PythiaSource",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(3.8e-05),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    crossSection = cms.untracked.double(37150000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(51)=7     ! structure function chosen', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=1     ! Defines the multi-parton model : 1 for gammajetwithbg', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.9409   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ', 
            'MSTP(91)=1     !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! '),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters'),
        processParameters = cms.vstring('MSEL=0       ', 
            'MSUB(14)=1   ', 
            'MSUB(29)=1   ', 
            'MSUB(114)=1  ', 
            'MSUB(115)=1  ', 
            'MSUB(11)=1   ', 
            'MSUB(12)=1   ', 
            'MSUB(13)=1   ', 
            'MSUB(15)=1   ', 
            'MSUB(16)=1   ', 
            'MSUB(18)=1   ', 
            'MSUB(19)=1   ', 
            'MSUB(20)=1   ', 
            'MSUB(28)=1   ', 
            'MSUB(30)=1   ', 
            'MSUB(31)=1   ', 
            'MSUB(53)=1   ', 
            'MSUB(68)=1   ', 
            'CKIN(3)=45          ! minimum pt hat for hard interactions', 
            'CKIN(4)=220          ! maximum pt hat for hard interactions')
    )
)


