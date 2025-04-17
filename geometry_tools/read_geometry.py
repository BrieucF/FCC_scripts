import dd4hep as dd4hepModule
from ROOT import dd4hep
import os
detector = dd4hep.Detector.getInstance()
detector.fromCompact(os.path.join(os.environ["K4GEO"], "FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"))
detector.detector("DCH_v2")
print(detector.constantAsDouble("DCH_alpha"))
