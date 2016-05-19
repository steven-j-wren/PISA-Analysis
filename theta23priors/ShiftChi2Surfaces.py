
import numpy as np
import scipy.interpolate

xNO2014, yNO2014 = np.loadtxt("Data/NuFit2014/sin2theta23dataNO.txt",unpack=True)
xIO2014, yIO2014 = np.loadtxt("Data/NuFit2014/sin2theta23dataIO.txt",unpack=True)
xLEMNO2015, yLEMNO2015 = np.loadtxt("Data/NuFit2015/sin2theta23dataLEMNO.txt",unpack=True)
xLEMIO2015, yLEMIO2015 = np.loadtxt("Data/NuFit2015/sin2theta23dataLEMIO.txt",unpack=True)
xLIDNO2015, yLIDNO2015 = np.loadtxt("Data/NuFit2015/sin2theta23dataLIDNO.txt",unpack=True)
xLIDIO2015, yLIDIO2015 = np.loadtxt("Data/NuFit2015/sin2theta23dataLIDIO.txt",unpack=True)

ShiftedyNO2014 = yNO2014 - np.min(yNO2014)
ShiftedyIO2014 = yIO2014 - np.min(yIO2014)
ShiftedyLEMNO2015 = yLEMNO2015 - np.min(yLEMNO2015)
ShiftedyLEMIO2015 = yLEMIO2015 - np.min(yLEMIO2015)
ShiftedyLIDNO2015 = yLIDNO2015 - np.min(yLIDNO2015)
ShiftedyLIDIO2015 = yLIDIO2015 - np.min(yLIDIO2015)

np.savetxt("Data/NuFit2014/sin2theta23dataShiftedNO.txt", np.c_[xNO2014,ShiftedyNO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2014/sin2theta23dataShiftedIO.txt", np.c_[xIO2014,ShiftedyIO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/sin2theta23dataShiftedLEMNO.txt", np.c_[xLEMNO2015,ShiftedyLEMNO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/sin2theta23dataShiftedLEMIO.txt", np.c_[xLEMIO2015,ShiftedyLEMIO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/sin2theta23dataShiftedLIDNO.txt", np.c_[xLIDNO2015,ShiftedyLIDNO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/sin2theta23dataShiftedLIDIO.txt", np.c_[xLIDIO2015,ShiftedyLIDIO2015], fmt=['%.3f', '%.7e'])

theta23xNO2014 = np.arcsin(np.sqrt(xNO2014))
theta23xIO2014 = np.arcsin(np.sqrt(xIO2014))
theta23xLEMNO2015 = np.arcsin(np.sqrt(xLEMNO2015))
theta23xLEMIO2015 = np.arcsin(np.sqrt(xLEMIO2015))
theta23xLIDNO2015 = np.arcsin(np.sqrt(xLIDNO2015))
theta23xLIDIO2015 = np.arcsin(np.sqrt(xLIDIO2015))

np.savetxt("Data/NuFit2014/theta23dataNO.txt", np.c_[theta23xNO2014,yNO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2014/theta23dataIO.txt", np.c_[theta23xIO2014,yIO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataLEMNO.txt", np.c_[theta23xLEMNO2015,yLEMNO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataLEMIO.txt", np.c_[theta23xLEMIO2015,yLEMIO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataLIDNO.txt", np.c_[theta23xLIDNO2015,yLIDNO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataLIDIO.txt", np.c_[theta23xLIDIO2015,yLIDIO2015], fmt=['%.3f', '%.7e'])

np.savetxt("Data/NuFit2014/theta23dataShiftedNO.txt", np.c_[theta23xNO2014,ShiftedyNO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2014/theta23dataShiftedIO.txt", np.c_[theta23xIO2014,ShiftedyIO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataShiftedLEMNO.txt", np.c_[theta23xLEMNO2015,ShiftedyLEMNO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataShiftedLEMIO.txt", np.c_[theta23xLEMIO2015,ShiftedyLEMIO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataShiftedLIDNO.txt", np.c_[theta23xLIDNO2015,ShiftedyLIDNO2015], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2015/theta23dataShiftedLIDIO.txt", np.c_[theta23xLIDIO2015,ShiftedyLIDIO2015], fmt=['%.3f', '%.7e'])
