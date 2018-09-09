import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
    #"file:./DYphotospp_Tune4C_13TeV_pythia8_cfi_py_GEN_VALIDATION.root"
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0EE4BD9B-7971-E011-AE11-1CC1DE1CF1BA.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/10A720F8-8B71-E011-9DA2-00237DA15C7C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1249CAEF-B471-E011-B206-0018FE29606E.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/12842656-9471-E011-A068-0017A477141C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/12B25738-F970-E011-ACCB-001F2969CD10.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/12F6D5D2-2371-E011-9BBA-1CC1DE1CDCAE.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/146F78E5-F870-E011-9D3D-1CC1DE1D014A.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1499C59E-8671-E011-AB6B-0017A4771014.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/160DDB46-3071-E011-98F9-0025B3E01FC2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/168F159D-9A71-E011-B043-0017A4770414.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/18103F29-A871-E011-A316-0017A4771030.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1870F388-2071-E011-AF5B-1CC1DE1D2028.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/18716A54-4671-E011-BD9C-1CC1DE1D036C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/18EBDA97-8171-E011-A25A-00237DA16692.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1AB6D2C4-1271-E011-83B4-0017A4770C20.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1AC63409-A271-E011-81E6-00237DA2F1DC.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1C3FBBE6-4071-E011-B478-0017A477080C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/1E382BB6-1F71-E011-8DB4-1CC1DE1CED1C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/20205AB5-A271-E011-8863-0017A4770410.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/20CEEF2C-7471-E011-9E04-1CC1DE1D16E0.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/24A1E017-3571-E011-9EDB-001E0B5FC57A.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/24E31CAF-8D71-E011-9582-0017A4770828.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/24F4E4E2-1771-E011-A01B-1CC1DE1CE170.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/26F5E436-8E71-E011-942D-00237DA41368.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/26FBC15E-4C71-E011-9520-00237DA13FC2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/280E23CB-8771-E011-ABDA-1CC1DE1D023A.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/28D480DC-9D71-E011-9976-001E0B5FD4A6.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2ABC281F-9870-E011-BDB7-78E7D1E4B772.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2AFE35FD-A071-E011-8924-1CC1DE1D16C8.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2C3FC99D-AC71-E011-9BB9-1CC1DE1CEDB2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2C56174E-5171-E011-B712-001E0BE922E2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2CACD981-EB71-E011-86DC-00237DA13C2E.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2E604CC3-0E71-E011-8B35-0017A4770010.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/2E97BCC4-AB71-E011-973F-1CC1DE046F78.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/3267531B-3B71-E011-8FE1-0017A4770030.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/326C908D-8B71-E011-8971-1CC1DE1CE01A.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/34B047AF-1C71-E011-8E9B-1CC1DE1CDDBC.root'
    )
)

##data Lumi filter
import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = './Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt').getVLuminosityBlockRange()
import FWCore.ParameterSet.Types as CfgTypes
myLumis = LumiList.LumiList(filename = './Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

process.demo = cms.EDAnalyzer('PhotonElectronAnalyzer',
 outputfileName = cms.untracked.string("PhotonElectron_Test.root"),
 triggerResultsTag = cms.untracked.InputTag("TriggerResults","","HLT"),
 vertexProducer  = cms.untracked.InputTag("offlinePrimaryVertices"),
 photonProducer  = cms.untracked.InputTag("photons"),                          
 electronProducer = cms.untracked.InputTag("gsfElectrons")
)

process.out_step = cms.EndPath(process.demo)
