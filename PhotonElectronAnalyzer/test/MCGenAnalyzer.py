import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
    "file:./DYphotospp_Tune4C_13TeV_pythia8_cfi_py_GEN_VALIDATION.root"
#     "file:./DYnorad_Tune4C_13TeV_pythia8_cfi_py_GEN_VALIDATION.root"
    )
)

process.demo = cms.EDAnalyzer('GenAnalyzer',
 outputfileName = cms.untracked.string("Gen_DY.root"),
 HepMCLabel = cms.untracked.InputTag('generatorSmeared','')
)

process.out_step = cms.EndPath(process.demo)
