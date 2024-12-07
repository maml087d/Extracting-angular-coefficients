theApp.EvtMax = -1

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [ 'EVNT.Zjj_LO_pol.merged.pool.root' ] # eventfile

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from xAODEventInfoCnv.xAODEventInfoCnvConf import xAODMaker__EventInfoCnvAlg
job += xAODMaker__EventInfoCnvAlg()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
import os
rivet.AnalysisPath = os.environ['PWD'] # directory for analysis path

rivet.Analyses += [ 'ATLAS_2020_I1803608:TYPE=EW_ONLY:CUTS=YES' ]
rivet.RunName = ''
rivet.HistoFile = 'Zjj_ptbin_data_xsectest_cuts.yoda.gz'
# rivet.CrossSection = 1.0
# rivet.IgnoreBeamCheck = True
# rivet.SkipWeights=True
job += rivet

