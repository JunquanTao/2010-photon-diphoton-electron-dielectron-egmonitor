from CMSOpenDataAnalysis.DecaysToPhotons.Sample import Sample

data_files = [
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/004DDBA5-7471-E011-A381-0017A4770C08.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/00C1B689-F670-E011-B054-1CC1DE1CF1BA.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/063C79CD-3271-E011-934D-0025B3E022C2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/065409FD-0071-E011-B185-001F296B758E.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/06C1851C-8E71-E011-83BB-00237DA16C42.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/08118EBE-9B71-E011-915E-001F296B758E.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0A6EBDBB-8371-E011-A8EE-0017A477003C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0C59D16B-F670-E011-A3D4-1CC1DE1CDDBC.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0C9AE5C4-3371-E011-85B7-1CC1DE1CEDB2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0CC571FD-7571-E011-98CB-1CC1DE046F78.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0E1ED66B-3971-E011-9C33-00237DA1EDE0.root'
]

data = Sample('data', False, data_files, 1)

sources = [data]
