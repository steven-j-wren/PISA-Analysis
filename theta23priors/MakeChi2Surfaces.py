
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import scipy.interpolate

xNO2014, yNO2014 = np.loadtxt("Data/NuFit2014/sin2theta23dataNO.txt",unpack=True)
xIO2014, yIO2014 = np.loadtxt("Data/NuFit2014/sin2theta23dataIO.txt",unpack=True)

x2014min = 0.4525

fNO = scipy.interpolate.splrep(xNO2014,yNO2014,s=0)
fIO = scipy.interpolate.splrep(xIO2014,yIO2014,s=0)

plt.plot(xNO2014,yNO2014,color='r',label='Normal')
plt.plot(xIO2014,yIO2014,color='b',label='Inverted')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.title(r'Nu-Fit 2014 $\Delta\chi^2$ projection for $\theta_{23}$',size='x-large')
plt.legend(loc='upper left')
plt.axvline(x2014min,color='r',linestyle='--',ymax=0.19)
plt.figtext(0.40,0.15,r'$\Delta\chi^2$(NO-IO)$\sim1$',color='k',size='xx-large')
plt.subplots_adjust(bottom=0.12)
plt.grid()
plt.savefig('Plots/NuFit2014Th23Chi2Surfaces.png')
plt.close()

xLEMNO2016, yLEMNO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLEMNO.txt",unpack=True)
xLEMIO2016, yLEMIO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLEMIO.txt",unpack=True)

xLEM2016min = 0.580

fLEMNO = scipy.interpolate.splrep(xLEMNO2016,yLEMNO2016,s=0)
fLEMIO = scipy.interpolate.splrep(xLEMIO2016,yLEMIO2016,s=0)

plt.plot(xLEMNO2016,yLEMNO2016,color='r',label='Normal')
plt.plot(xLEMIO2016,yLEMIO2016,color='b',label='Inverted')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.title(r'Nu-Fit 2016 $\Delta\chi^2$ projection for $\theta_{23}$ (NO$\nu$A LEM Data)',size='x-large')
plt.legend(loc='upper left')
plt.axvline(xLEM2016min,color='b',linestyle='--',ymax=0.19)
plt.figtext(0.50,0.15,r'$\Delta\chi^2$(IO-NO)$\sim1$',color='k',size='x-large')
plt.subplots_adjust(bottom=0.12)
plt.grid()
plt.savefig('Plots/NuFit2016LEMTh23Chi2Surfaces.png')
plt.close()

xLIDNO2016, yLIDNO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLIDNO.txt",unpack=True)
xLIDIO2016, yLIDIO2016 = np.loadtxt("Data/NuFit2016/sin2theta23dataLIDIO.txt",unpack=True)

xLID2016min = 0.450

fLIDNO = scipy.interpolate.splrep(xLIDNO2016,yLIDNO2016,s=0)
fLIDIO = scipy.interpolate.splrep(xLIDIO2016,yLIDIO2016,s=0)

plt.plot(xLIDNO2016,yLIDNO2016,color='r',label='Normal')
plt.plot(xLIDIO2016,yLIDIO2016,color='b',label='Inverted')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.title(r'Nu-Fit 2016 $\Delta\chi^2$ projection for $\theta_{23}$ (NO$\nu$A LID Data)',size='x-large')
plt.legend(loc='upper left')
plt.axvline(xLID2016min,color='r',linestyle='--',ymax=0.11)
plt.figtext(0.15,0.15,r'$\Delta\chi^2$(NO-IO)$\sim0.5$',color='k',size='x-large')
plt.subplots_adjust(bottom=0.12)
plt.grid()
plt.savefig('Plots/NuFit2016LIDTh23Chi2Surfaces.png')
plt.close()

plt.plot(xNO2014,yNO2014,color='r',label='Normal, NuFit 2014')
plt.plot(xIO2014,yIO2014,color='b',label='Inverted, NuFit 2014')
plt.plot(xLIDNO2016,yLIDNO2016,color='r',linestyle='--',label='Normal, NuFit 2016 (LID)')
plt.plot(xLIDIO2016,yLIDIO2016,color='b',linestyle='--',label='Inverted, NuFit 2016 (LID)')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.subplots_adjust(bottom=0.12,top=0.8)
plt.title(r'Nu-Fit $\Delta\chi^2$ projections for $\theta_{23}$',size='x-large',x=0.5,y=1.13)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,fontsize='small')
plt.grid()
plt.savefig('Plots/NuFit20142016LIDTh23Chi2Surfaces.png')
plt.close()

plt.plot(xNO2014,yNO2014,color='r',label='Normal, NuFit 2014')
plt.plot(xIO2014,yIO2014,color='b',label='Inverted, NuFit 2014')
plt.plot(xLEMNO2016,yLEMNO2016,color='r',linestyle=':',label='Normal, NuFit 2016 (LEM)')
plt.plot(xLEMIO2016,yLEMIO2016,color='b',linestyle=':',label='Inverted, NuFit 2016 (LEM)')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.subplots_adjust(bottom=0.12,top=0.8)
plt.title(r'Nu-Fit $\Delta\chi^2$ projections for $\theta_{23}$',size='x-large',x=0.5,y=1.13)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,fontsize='small')
plt.grid()
plt.savefig('Plots/NuFit20142016LEMTh23Chi2Surfaces.png')
plt.close()

plt.plot(xLIDNO2016,yLIDNO2016,color='r',linestyle='--',label='Normal, NuFit 2016 (LID)')
plt.plot(xLIDIO2016,yLIDIO2016,color='b',linestyle='--',label='Inverted, NuFit 2016 (LID)')
plt.plot(xLEMNO2016,yLEMNO2016,color='r',linestyle=':',label='Normal, NuFit 2016 (LEM)')
plt.plot(xLEMIO2016,yLEMIO2016,color='b',linestyle=':',label='Inverted, NuFit 2016 (LEM)')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.subplots_adjust(bottom=0.12,top=0.8)
plt.title(r'Nu-Fit $\Delta\chi^2$ projections for $\theta_{23}$',size='x-large',x=0.5,y=1.13)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,fontsize='small')
plt.grid()
plt.savefig('Plots/NuFit2016LID2016LEMTh23Chi2Surfaces.png')
plt.close()

plt.plot(xNO2014,yNO2014,color='r',label='Normal, NuFit 2014')
plt.plot(xIO2014,yIO2014,color='b',label='Inverted, NuFit 2014')
plt.plot(xLIDNO2016,yLIDNO2016,color='r',linestyle='--',label='Normal, NuFit 2016 (LID)')
plt.plot(xLIDIO2016,yLIDIO2016,color='b',linestyle='--',label='Inverted, NuFit 2016 (LID)')
plt.plot(xLEMNO2016,yLEMNO2016,color='r',linestyle=':',label='Normal, NuFit 2016 (LEM)')
plt.plot(xLEMIO2016,yLEMIO2016,color='b',linestyle=':',label='Inverted, NuFit 2016 (LEM)')
plt.axis([0.35,0.65,0,5])
plt.xlabel(r'$\sin^2\theta_{23}$',fontsize=30, labelpad=10)
plt.ylabel(r'$\Delta\chi^2$',fontsize=30, labelpad=10)
plt.subplots_adjust(bottom=0.12,top=0.8)
plt.title(r'Nu-Fit $\Delta\chi^2$ projections for $\theta_{23}$',size='x-large',x=0.5,y=1.13)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.,fontsize='small')
plt.grid()
plt.savefig('Plots/NuFit20142016LID2016LEMTh23Chi2Surfaces.png')
plt.close()


