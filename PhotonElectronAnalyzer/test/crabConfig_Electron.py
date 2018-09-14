from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Run2010BPhoton_EGAnalyzer_Electron'
config.General.workArea = 'Run2010BPhoton_EGAnalyzer_Electron'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PhotonElectronAnalyzer.py' 
config.JobType.outputFiles = ['PhotonElectron_Test.root']
config.Data.inputDataset = '/Electron/Run2010B-Apr21ReReco-v1/AOD'
#config.Data.userInputFiles = ['root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0000/0003D1EF-B071-E011-BA0A-78E7D1651098.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0000/001CE57A-F470-E011-8208-001CC4A934D8.root']
config.Data.publication = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =  10

config.Data.outLFNDirBase = "/store/user/jtao/CMSOpenData2010/"
config.Site.storageSite = 'T2_CH_CERN'  # you might need to change this to a site you have acces too
config.section_("Debug")

