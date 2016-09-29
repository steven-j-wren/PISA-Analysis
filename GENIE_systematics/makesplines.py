
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import scipy.interpolate

from pisa.utils.jsons import from_json

indict = from_json('/Users/steven/IceCube/PISA/pisa/pisa/resources/aeff/relative-crosssections/genie-weigthed-crosssections.json')

axislabels = {}
axislabels['nue'] = r'$\nu_e$ CC'
axislabels['nue_bar'] = r'$\bar{\nu}_e$ CC'
axislabels['numu'] = r'$\nu_{\mu}$ CC'
axislabels['numu_bar'] = r'$\bar{\nu}_{\mu}$ CC'
axislabels['nutau'] = r'$\nu_{\tau}$ CC'
axislabels['nutau_bar'] = r'$\bar{\nu}_{\tau}$ CC'
axislabels['nuall'] = r'$\nu$ NC'
axislabels['nuallbar'] = r'$\bar{\nu}$ NC'

GENIElabels = {}
GENIElabels['MaRES'] = r'$M_A^{RES}$'
GENIElabels['MaCCQE'] = r'$M_A^{CCQE}$'
GENIElabels['AhtBY'] = r'$A_{HT}$'
GENIElabels['BhtBY'] = r'$B_{BT}$'
GENIElabels['CV1uBY'] = r'$C_{\nu1u}$'
GENIElabels['CV2uBY'] = r'$C_{\nu2u}$'

GENIEcolours = {}
GENIEcolours['MaRES'] = 'b'
GENIEcolours['MaCCQE'] = 'r'
GENIEcolours['AhtBY'] = 'g'
GENIEcolours['BhtBY'] = 'c'
GENIEcolours['CV1uBY'] = 'k'
GENIEcolours['CV2uBY'] = 'm'

TauResults = {}
MuResults = {}
EResults = {}
NCResults = {}
TauBarResults = {}
MuBarResults = {}
EBarResults = {}
NCBarResults = {}

for parameter in indict.keys():
    
    xvals = []
    yvals = []
    
    for energyval in indict[parameter].keys():
        xvals.append(float(energyval))
    xvals = sorted(xvals)
    for energyval in xvals:
        yvals.append(float(indict[parameter][str(energyval)]))

    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)

    paramlabel = ''

    for flavour in axislabels.keys():
        if flavour in parameter:
            if 'bar' in flavour and 'bar' in parameter:
                paramlabel += axislabels[flavour]
                flavourlabel = flavour
            if 'bar' not in flavour and 'bar' not in parameter:
                paramlabel += axislabels[flavour]
                flavourlabel = flavour

    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,color=GENIEcolour,label='Normal')
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty due to %s'%paramlabel)
    plt.savefig('%sUncertaintySpline.png'%parameter)

    plt.close()

    if flavourlabel == 'nutau':
        TauResults[parameter] = {}
        TauResults[parameter]['x'] = xvals
        TauResults[parameter]['y'] = yvals
        TauResults[parameter]['s'] = fParam

    if flavourlabel == 'numu':
        MuResults[parameter] = {}
        MuResults[parameter]['x'] = xvals
        MuResults[parameter]['y'] = yvals
        MuResults[parameter]['s'] = fParam

    if flavourlabel == 'nue':
        EResults[parameter] = {}
        EResults[parameter]['x'] = xvals
        EResults[parameter]['y'] = yvals
        EResults[parameter]['s'] = fParam

    if flavourlabel == 'nuall':
        NCResults[parameter] = {}
        NCResults[parameter]['x'] = xvals
        NCResults[parameter]['y'] = yvals
        NCResults[parameter]['s'] = fParam

    if flavourlabel == 'nutau_bar':
        TauBarResults[parameter] = {}
        TauBarResults[parameter]['x'] = xvals
        TauBarResults[parameter]['y'] = yvals
        TauBarResults[parameter]['s'] = fParam

    if flavourlabel == 'numu_bar':
        MuBarResults[parameter] = {}
        MuBarResults[parameter]['x'] = xvals
        MuBarResults[parameter]['y'] = yvals
        MuBarResults[parameter]['s'] = fParam

    if flavourlabel == 'nue_bar':
        EBarResults[parameter] = {}
        EBarResults[parameter]['x'] = xvals
        EBarResults[parameter]['y'] = yvals
        EBarResults[parameter]['s'] = fParam

    if flavourlabel == 'nuallbar':
        NCBarResults[parameter] = {}
        NCBarResults[parameter]['x'] = xvals
        NCBarResults[parameter]['y'] = yvals
        NCBarResults[parameter]['s'] = fParam

for parameter in TauResults.keys():
    xvals = TauResults[parameter]['x']
    yvals = TauResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['nutau']
    flavourlabel = axislabels['nutau']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]

    plt.plot(xvals,yvals,label=GENIElabel)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NuTauUncertaintySplines.png')
plt.close()

for parameter in MuResults.keys():
    xvals = MuResults[parameter]['x']
    yvals = MuResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['numu']
    flavourlabel = axislabels['numu']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NuMuUncertaintySplines.png')
plt.close()

for parameter in EResults.keys():
    xvals = EResults[parameter]['x']
    yvals = EResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['nue']
    flavourlabel = axislabels['nue']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NuEUncertaintySplines.png')
plt.close()

for parameter in NCResults.keys():
    xvals = NCResults[parameter]['x']
    yvals = NCResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['nuall']
    flavourlabel = axislabels['nuall']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NCUncertaintySplines.png')
plt.close()

for parameter in TauBarResults.keys():
    xvals = TauBarResults[parameter]['x']
    yvals = TauBarResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['nutau_bar']
    flavourlabel = axislabels['nutau_bar']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NuTauBarUncertaintySplines.png')
plt.close()

for parameter in MuBarResults.keys():
    xvals = MuBarResults[parameter]['x']
    yvals = MuBarResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['numu_bar']
    flavourlabel = axislabels['numu_bar']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NuMuBarUncertaintySplines.png')
plt.close()

for parameter in EBarResults.keys():
    xvals = EBarResults[parameter]['x']
    yvals = EBarResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['nue_bar']
    flavourlabel = axislabels['nue_bar']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NuEBarUncertaintySplines.png')
plt.close()

for parameter in NCBarResults.keys():
    xvals = NCBarResults[parameter]['x']
    yvals = NCBarResults[parameter]['y']
    fParam = scipy.interpolate.splrep(np.array(xvals),np.array(yvals),s=0)
    paramlabel = ''
    paramlabel += axislabels['nuallbar']
    flavourlabel = axislabels['nuallbar']
    for GEN in GENIElabels.keys():
        if GEN in parameter:
            paramlabel += ' ' + GENIElabels[GEN]
            GENIElabel = GENIElabels[GEN]
            GENIEcolour = GENIEcolours[GEN]

    plt.plot(xvals,yvals,label=GENIElabel,color=GENIEcolour)
    plt.xlabel(r'Energy [GeV]')
    plt.ylabel(r'Relative Uncertainty on Cross Section')
    plt.title(r'Energy-Dependent Uncertainty on %s'%flavourlabel)

plt.grid()
plt.legend(loc='best')
plt.savefig('NCBarUncertaintySplines.png')
plt.close()

