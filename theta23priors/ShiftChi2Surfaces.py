
import numpy as np
import scipy.interpolate

xNO2014, yNO2014 = np.loadtxt("Data/NuFit2014/sin2theta23dataNO.txt",unpack=True)
xIO2014, yIO2014 = np.loadtxt("Data/NuFit2014/sin2theta23dataIO.txt",unpack=True)
xLEMNO2016, yLEMNO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLEMNO.txt",unpack=True)
xLEMIO2016, yLEMIO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLEMIO.txt",unpack=True)
xLIDNO2016, yLIDNO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLIDNO.txt",unpack=True)
xLIDIO2016, yLIDIO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLIDIO.txt",unpack=True)

ShiftedyNO2014 = yNO2014 - np.min(yNO2014)
ShiftedyIO2014 = yIO2014 - np.min(yIO2014)
ShiftedyLEMNO2016 = yLEMNO2016 - np.min(yLEMNO2016)
ShiftedyLEMIO2016 = yLEMIO2016 - np.min(yLEMIO2016)
ShiftedyLIDNO2016 = yLIDNO2016 - np.min(yLIDNO2016)
ShiftedyLIDIO2016 = yLIDIO2016 - np.min(yLIDIO2016)

print np.arcsin(np.sqrt(xNO2014[np.where(ShiftedyNO2014==np.min(ShiftedyNO2014))][0]))*180/np.pi
print np.arcsin(np.sqrt(xIO2014[np.where(ShiftedyIO2014==np.min(ShiftedyIO2014))][0]))*180/np.pi
print np.arcsin(np.sqrt(xLEMNO2016[np.where(ShiftedyLEMNO2016==np.min(ShiftedyLEMNO2016))][0]))*180/np.pi
print np.arcsin(np.sqrt(xLEMIO2016[np.where(ShiftedyLEMIO2016==np.min(ShiftedyLEMIO2016))][0]))*180/np.pi
print np.arcsin(np.sqrt(xLIDNO2016[np.where(ShiftedyLIDNO2016==np.min(ShiftedyLIDNO2016))][0]))*180/np.pi
print np.arcsin(np.sqrt(xLIDIO2016[np.where(ShiftedyLIDIO2016==np.min(ShiftedyLIDIO2016))][0]))*180/np.pi

np.savetxt("Data/NuFit2014/sin2theta23dataShiftedNO.txt", np.c_[xNO2014,ShiftedyNO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2014/sin2theta23dataShiftedIO.txt", np.c_[xIO2014,ShiftedyIO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/sin2theta23dataShiftedLEMNO.txt", np.c_[xLEMNO2016,ShiftedyLEMNO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/sin2theta23dataShiftedLEMIO.txt", np.c_[xLEMIO2016,ShiftedyLEMIO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/sin2theta23dataShiftedLIDNO.txt", np.c_[xLIDNO2016,ShiftedyLIDNO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/sin2theta23dataShiftedLIDIO.txt", np.c_[xLIDIO2016,ShiftedyLIDIO2016], fmt=['%.3f', '%.7e'])

theta23xNO2014 = np.arcsin(np.sqrt(xNO2014))
theta23xIO2014 = np.arcsin(np.sqrt(xIO2014))
theta23xLEMNO2016 = np.arcsin(np.sqrt(xLEMNO2016))
theta23xLEMIO2016 = np.arcsin(np.sqrt(xLEMIO2016))
theta23xLIDNO2016 = np.arcsin(np.sqrt(xLIDNO2016))
theta23xLIDIO2016 = np.arcsin(np.sqrt(xLIDIO2016))

np.savetxt("Data/NuFit2014/theta23dataNO.txt", np.c_[theta23xNO2014,yNO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2014/theta23dataIO.txt", np.c_[theta23xIO2014,yIO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataLEMNO.txt", np.c_[theta23xLEMNO2016,yLEMNO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataLEMIO.txt", np.c_[theta23xLEMIO2016,yLEMIO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataLIDNO.txt", np.c_[theta23xLIDNO2016,yLIDNO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataLIDIO.txt", np.c_[theta23xLIDIO2016,yLIDIO2016], fmt=['%.3f', '%.7e'])

np.savetxt("Data/NuFit2014/theta23dataShiftedNO.txt", np.c_[theta23xNO2014,ShiftedyNO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2014/theta23dataShiftedIO.txt", np.c_[theta23xIO2014,ShiftedyIO2014], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataShiftedLEMNO.txt", np.c_[theta23xLEMNO2016,ShiftedyLEMNO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataShiftedLEMIO.txt", np.c_[theta23xLEMIO2016,ShiftedyLEMIO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataShiftedLIDNO.txt", np.c_[theta23xLIDNO2016,ShiftedyLIDNO2016], fmt=['%.3f', '%.7e'])
np.savetxt("Data/NuFit2016/theta23dataShiftedLIDIO.txt", np.c_[theta23xLIDIO2016,ShiftedyLIDIO2016], fmt=['%.3f', '%.7e'])
