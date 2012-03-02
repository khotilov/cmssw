import FWCore.ParameterSet.Config as cms

herwigppUESettingsBlock = cms.PSet(
	herwigppUE_EE_3C_Base = cms.vstring(
		'+pdfCTEQ6L',

		'cd /Herwig',
		'create Herwig::O2AlphaS O2AlphaS',
		'set Model:QCD/RunningAlphaS O2AlphaS',

		# Energy-independent MPI parameters
		#   Colour reconnection settings
		'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes',
		'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.54',
		#   Colour Disrupt settings
		'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80',
		#   inverse hadron radius
		'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.11',
		#   MPI model settings
		#'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0',
		'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes',
		'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes',
		'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2',
	),

	herwigppUE_EE_3C_900GeV = cms.vstring(
		'+herwigppUE_EE_3C_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 900.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 1.55',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 3.1',
		'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 1.81*GeV',
	),

	herwigppUE_EE_3C_1800GeV = cms.vstring(
		'+herwigppUE_EE_3C_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 1800.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 2.26',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 4.52',
	),

	herwigppUE_EE_3C_2760GeV = cms.vstring(
		'+herwigppUE_EE_3C_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 2760.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 2.33',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 4.66',
	),

	herwigppUE_EE_3C_7000GeV = cms.vstring(
		'+herwigppUE_EE_3C_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 7000.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 2.752',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 5.504',
		'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.34*GeV',
	),

	herwigppUE_EE_3C_8000GeV = cms.vstring(
		'+herwigppUE_EE_3C_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 8000.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 2.85',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 5.7',
	),

	herwigppUE_EE_3C_14000GeV = cms.vstring(
		'+herwigppUE_EE_3C_Base',
		'set LHCGenerator:EventHandler:LuminosityFunction:Energy 14000.0',
		'set /Herwig/UnderlyingEvent/KtCut:MinKT 3.53',
		'set /Herwig/UnderlyingEvent/UECuts:MHatMin 7.06',
	),
)
