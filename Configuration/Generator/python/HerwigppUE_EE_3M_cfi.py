import FWCore.ParameterSet.Config as cms

herwigppUESettingsBlock = cms.PSet(
	herwigppUE_EE_3M_Base = cms.vstring(
		'+pdfMRST2008LOss',

		'cd /Herwig',
		'create Herwig::O2AlphaS O2AlphaS',
		'set Model:QCD/RunningAlphaS O2AlphaS',

		# Energy-independent MPI parameters
		#   Colour reconnection settings
		'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes',
		'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.61',
		#   Colour Disrupt settings
		'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.75',
		#   inverse hadron radius
		'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.35',
		#   MPI model settings
		#'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0',
		'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes',
		'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes',
		'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2',
	),

	herwigppUE_EE_3M_900GeV = cms.vstring(
		'+herwigppUE_EE_3M_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 900.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 1.86',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 3.72',
	),

	herwigppUE_EE_3M_1800GeV = cms.vstring(
		'+herwigppUE_EE_3M_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 1800.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 2.55',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 5.1',
	),

	herwigppUE_EE_3M_2760GeV = cms.vstring(
		'+herwigppUE_EE_3M_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 2760.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 2.62',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 5.24',
	),

	herwigppUE_EE_3M_7000GeV = cms.vstring(
		'+herwigppUE_EE_3M_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 7000.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 3.06',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 6.12',
		'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.37*GeV',
	),

	herwigppUE_EE_3M_8000GeV = cms.vstring(
		'+herwigppUE_EE_3M_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 8000.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 3.21',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 6.42',
	),

	herwigppUE_EE_3M_14000GeV = cms.vstring(
		'+herwigppUE_EE_3M_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 14000.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 3.53',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 7.06',
	),
)
