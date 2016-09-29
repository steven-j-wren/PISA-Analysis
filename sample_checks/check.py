#!/usr/bin/env python
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
from icecube import dataclasses, icetray
from icecube.dataio import I3File
from icecube.icetray import I3Frame, I3Units

from msu import passed_msu

#msu = "/data/user/hignight/nue/12600/"
#nbi = "/data/user/mlarson/tau_appearance_samples/level7_24Nov2015/12600/"

msu = "/data/user/hignight/numu/14600/"
nbi = "/data/user/mlarson/tau_appearance/scripts/level7/output/14600/"
#nbi = "/data/user/mlarson/tau_appearance_samples/level7_24Nov2015/14600/"

#msu = "/data/user/hignight/nutau/16600/"
#nbi = "/data/user/mlarson/tau_appearance_samples/level7_24Nov2015/16600/"

filenames = {}
eventids = {}
weight = {}
energy = {}
coszen = {}

reco_energy = {}
reco_coszen = {}

filenames[msu] = glob.glob(msu + "*")
filenames[msu].sort()
filenames[msu] = filenames[msu][:10]

filenames[nbi] = glob.glob(nbi + "*")
filenames[nbi].sort()
filenames[nbi] = filenames[nbi][:1000]

for x in filenames.keys():
    eventids[x] = {}
    weight[x] = {}
    energy[x] = {}
    coszen[x] = {}

    reco_energy[x] = {}
    reco_coszen[x] = {}

    for filename in filenames[x]:
        print filename

        input = I3File(filename, 'r')

        if x==msu:
            fnumber = int(os.path.basename(filename).split(".")[3])
        else:
            fnumber = int(os.path.basename(filename).split(".")[2])

        eventids[x][fnumber] = []
        weight[x][fnumber] = []
        energy[x][fnumber] = []
        coszen[x][fnumber] = []

        reco_energy[x][fnumber] = []
        reco_coszen[x][fnumber] = []
        
        while input.more():
            frame = input.pop_physics()

            if x==msu:                 
                if not passed_msu(frame): continue
                w = frame["NeutrinoWeights"]["UnoscillatedRate"]
                cascfit = frame['IC86_Dunkman_L6_PegLeg_MultiNest8D_HDCasc']
                trackfit = frame["IC86_Dunkman_L6_PegLeg_MultiNest8D_Track"]
            else:
                if not frame["TauL6_bool"].value: continue
                w = frame["OcelotWeight"]["UnoscillatedRate"]
                print f
                cascfit = frame["Pegleg_Fit_MNHDCasc"]
                trackfit = frame["Pegleg_Fit_MNTrack"]

            p = frame["I3MCTree"].most_energetic_neutrino
            eventid = frame["I3EventHeader"].event_id            
            
            isCC = frame["I3MCWeightDict"]["InteractionType"]==1

            if not isCC: continue

            weight[x][fnumber].append(w)
            energy[x][fnumber].append(p.energy)
            coszen[x][fnumber].append(numpy.cos(p.dir.zenith))
            eventids[x][fnumber].append(eventid)        

            reco_energy[x][fnumber].append(trackfit.energy + cascfit.energy)
            reco_coszen[x][fnumber].append(numpy.cos(trackfit.dir.zenith))

        input.close()
        eventids[x][fnumber] = numpy.array(eventids[x][fnumber])
        weight[x][fnumber] = numpy.array(weight[x][fnumber])
        energy[x][fnumber] = numpy.array(energy[x][fnumber])
        coszen[x][fnumber] = numpy.array(coszen[x][fnumber])

        reco_energy[x][fnumber] = numpy.array(reco_energy[x][fnumber])
        reco_coszen[x][fnumber] = numpy.array(reco_coszen[x][fnumber])
    



overlaps = {msu:[],
            nbi:[]}
energies = {msu:[],
            nbi:[]}
coszens = {msu:[],
           nbi:[]}
reco_energies = {msu:[],
                 nbi:[]}
reco_coszens = {msu:[],
                nbi:[]}

# do the comparisons
print "#\tNBI\t MSU\t Overlap"
for n in eventids[nbi].keys():
    if n in eventids[msu].keys():
        msu_overlap = []
        nbi_overlap = []
        for x in eventids[nbi][n]:
            if x in eventids[msu][n]: nbi_overlap.append(True)
            else: nbi_overlap.append(False)
        for x in eventids[msu][n]:
            if x in eventids[nbi][n]: msu_overlap.append(True)
            else: msu_overlap.append(False)

        overlap = numpy.array(nbi_overlap)
        overlaps[nbi].extend(nbi_overlap)
        overlaps[msu].extend(msu_overlap)
        
        energies[nbi].extend(energy[nbi][n])
        energies[msu].extend(energy[msu][n])
        coszens[nbi].extend(coszen[nbi][n])
        coszens[msu].extend(coszen[msu][n])
        reco_energies[nbi].extend(reco_energy[nbi][n])
        reco_energies[msu].extend(reco_energy[msu][n])
        reco_coszens[nbi].extend(reco_coszen[nbi][n])
        reco_coszens[msu].extend(reco_coszen[msu][n])
        


        #print "%i\t%i (%4.2e)\t %i (%4.2e)\t %i (%4.2e)" % (n, 
        #                                                    len(eventids[nbi][n]),
        #                                                    numpy.sum(weight[nbi][n]),
        #                                                    len(eventids[msu][n]),
        #                                                    numpy.sum(weight[msu][n]),
        #                                                    numpy.sum(overlap),
        #                                                    numpy.sum(weight[nbi][n][overlap]))
        

overlaps[nbi] = numpy.array(overlaps[nbi])
overlaps[msu] = numpy.array(overlaps[msu])

energies[nbi] = numpy.array(energies[nbi])
energies[msu] = numpy.array(energies[msu])

coszens[nbi] = numpy.array(coszens[nbi])
coszens[msu] = numpy.array(coszens[msu])

reco_coszens[nbi] = numpy.array(reco_coszens[nbi])
reco_coszens[msu] = numpy.array(reco_coszens[msu])

reco_energies[nbi] = numpy.array(reco_energies[nbi])
reco_energies[msu] = numpy.array(reco_energies[msu])

# Plot the reco
pyplot.figure()
pyplot.hist(numpy.log10(reco_energies[msu][numpy.logical_not(overlaps[msu])]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="MSU (Non-overlapping)",
            color='g',
            histtype='step',
            linewidth=2)
pyplot.hist(numpy.log10(reco_energies[msu][overlaps[msu]]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="MSU (Overlapping)",
            color='g',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.hist(numpy.log10(reco_energies[nbi][numpy.logical_not(overlaps[nbi])]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="NBI (Non-overlapping)",
            color='b',
            histtype='step',
            linewidth=2)
pyplot.hist(numpy.log10(reco_energies[nbi][overlaps[nbi]]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="NBI (Overlapping)",
            color='b',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.xlabel("Reco Energy")
pyplot.ylabel("A.U.")
pyplot.legend(loc='upper left')
pyplot.grid()
pyplot.savefig("reco_energy.pdf")

# Plot the truth
pyplot.figure()
pyplot.hist(numpy.log10(energies[msu][numpy.logical_not(overlaps[msu])]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="MSU (Non-overlapping)",
            color='g',
            histtype='step',
            linewidth=2)
pyplot.hist(numpy.log10(energies[msu][overlaps[msu]]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="MSU (Overlapping)",
            color='g',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.hist(numpy.log10(energies[nbi][numpy.logical_not(overlaps[nbi])]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="NBI (Non-overlapping)",
            color='b',
            histtype='step',
            linewidth=2)
pyplot.hist(numpy.log10(energies[nbi][overlaps[nbi]]), 
            bins = numpy.linspace(-0.5, 2, 30),
            label="NBI (Overlapping)",
            color='b',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.xlabel("True Energy")
pyplot.ylabel("A.U.")
pyplot.legend(loc='upper left')
pyplot.grid()
pyplot.savefig("true_energy.pdf")


# Plot the delta vs true energy
f, [[ax1, ax2], [ax3, ax4]] = pyplot.subplots(2,2)
h1, xedges, yedges, junk = ax1.hist2d(energies[msu][numpy.logical_not(overlaps[msu])],
                                      (reco_energies[msu]/energies[msu])[numpy.logical_not(overlaps[msu])],
                                      bins = (numpy.linspace(0, 100, 50), numpy.linspace(0, 2, 50)),
                                      )
h2, xedges, yedges, junk = ax2.hist2d(energies[msu][overlaps[msu]],
                                      (reco_energies[msu]/energies[msu])[overlaps[msu]],
                                      bins = (numpy.linspace(0, 100, 50), numpy.linspace(0, 2, 50)),
                                      )
h3, xedges, yedges, junk = ax3.hist2d(energies[nbi][numpy.logical_not(overlaps[nbi])],
                                      (reco_energies[nbi]/energies[nbi])[numpy.logical_not(overlaps[nbi])],
                                      bins = (numpy.linspace(0, 100, 50), numpy.linspace(0, 2, 50)),
                                      )
h4, xedges, yedges, junk = ax4.hist2d(energies[nbi][overlaps[nbi]],
                                      (reco_energies[nbi]/energies[nbi])[overlaps[nbi]],
                                      bins = (numpy.linspace(0, 100, 50), numpy.linspace(0, 2, 50)),
                                      )

h1_sigmas = []
h2_sigmas = []
h3_sigmas = []
h4_sigmas = []

for bin in range(h1.shape[0]):
    h1[bin,:] /= numpy.sum(h1[bin,:])
    h2[bin,:] /= numpy.sum(h2[bin,:])
    h3[bin,:] /= numpy.sum(h3[bin,:])
    h4[bin,:] /= numpy.sum(h4[bin,:])
    
    h1_cum = numpy.cumsum(h1[bin,:])
    h2_cum = numpy.cumsum(h2[bin,:])
    h3_cum = numpy.cumsum(h3[bin,:])
    h4_cum = numpy.cumsum(h4[bin,:])

    h1_sigmas.append( numpy.array([yedges[h1_cum > 0.16][0], yedges[h1_cum <= 0.84][-1]]) )
    h2_sigmas.append( numpy.array([yedges[h2_cum > 0.16][0], yedges[h2_cum <= 0.84][-1]]) )
    h3_sigmas.append( numpy.array([yedges[h3_cum > 0.16][0], yedges[h3_cum <= 0.84][-1]]) )
    h4_sigmas.append( numpy.array([yedges[h4_cum > 0.16][0], yedges[h4_cum <= 0.84][-1]]) )

h1_sigmas = numpy.array(h1_sigmas).T
h2_sigmas = numpy.array(h2_sigmas).T
h3_sigmas = numpy.array(h3_sigmas).T
h4_sigmas = numpy.array(h4_sigmas).T

f, [[ax1, ax2], [ax3, ax4]] = pyplot.subplots(2,2)
ax1.pcolormesh(xedges, yedges, h1.T)
ax2.pcolormesh(xedges, yedges, h2.T)
ax3.pcolormesh(xedges, yedges, h3.T)
ax4.pcolormesh(xedges, yedges, h4.T)

ax1.plot(xedges[:-1], h1_sigmas[0], color='w')
ax1.plot(xedges[:-1], h1_sigmas[1], color='w')
ax2.plot(xedges[:-1], h2_sigmas[0], color='w')
ax2.plot(xedges[:-1], h2_sigmas[1], color='w')
ax3.plot(xedges[:-1], h3_sigmas[0], color='w')
ax3.plot(xedges[:-1], h3_sigmas[1], color='w')
ax4.plot(xedges[:-1], h4_sigmas[0], color='w')
ax4.plot(xedges[:-1], h4_sigmas[1], color='w')

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax1.set_ylabel("Reco/True Energy")
ax3.set_ylabel("Reco/True Energy")
ax3.set_xlabel("True Energy")
ax4.set_xlabel("True Energy")
pyplot.tight_layout()
pyplot.savefig("energy_resolutions.pdf")

# Plot the truth
pyplot.figure()
pyplot.hist(coszens[msu][numpy.logical_not(overlaps[msu])], 
            bins = numpy.linspace(-1, 1, 30),
            label="MSU (Non-overlapping)",
            color='g',
            histtype='step',
            linewidth=2)
pyplot.hist(coszens[msu][overlaps[msu]], 
            bins = numpy.linspace(-1, 1, 30),
            label="MSU (Overlapping)",
            color='g',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.hist(coszens[nbi][numpy.logical_not(overlaps[nbi])], 
            bins = numpy.linspace(-1, 1, 30),
            label="NBI (Non-overlapping)",
            color='b',
            histtype='step',
            linewidth=2)
pyplot.hist(coszens[nbi][overlaps[nbi]], 
            bins = numpy.linspace(-1, 1, 30),
            label="NBI (Overlapping)",
            color='b',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.xlabel("True Cos(Zenith)")
pyplot.ylabel("A.U.")
pyplot.legend()
pyplot.grid()
pyplot.savefig("true_coszen.pdf")




# Plot the truth
pyplot.figure()
pyplot.hist(reco_coszens[msu][numpy.logical_not(overlaps[msu])], 
            bins = numpy.linspace(-1, 1, 30),
            label="MSU (Non-overlapping)",
            color='g',
            histtype='step',
            linewidth=2)
pyplot.hist(reco_coszens[msu][overlaps[msu]], 
            bins = numpy.linspace(-1, 1, 30),
            label="MSU (Overlapping)",
            color='g',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.hist(reco_coszens[nbi][numpy.logical_not(overlaps[nbi])], 
            bins = numpy.linspace(-1, 1, 30),
            label="NBI (Non-overlapping)",
            color='b',
            histtype='step',
            linewidth=2)
pyplot.hist(reco_coszens[nbi][overlaps[nbi]], 
            bins = numpy.linspace(-1, 1, 30),
            label="NBI (Overlapping)",
            color='b',
            histtype='step',
            linestyle='dashed',
            linewidth=2)
pyplot.xlabel("Reco Cos(Zenith)")
pyplot.ylabel("A.U.")
pyplot.legend()
pyplot.grid()
pyplot.savefig("reco_coszen.pdf")
