#!/usr/bin/env python
import os, sys, pickle, numpy

from icecube import dataclasses, icetray, MuonGun
from icecube.icetray import I3Frame, I3Units 
from icecube.dataio import I3File
from icecube.MuonGun import Cylinder

infiles = ["Level5_corsika.009036.100000files.i3",
           "Level5_corsika.009255.95000files.i3"]
	
weights = []
energy_of_primary = []
energy_at_depth = []
coszens = []
code_of_primary = []
type_of_primary = []

id_1000160320 = False
id_1000080160 = False
id_2212 = False
id_1000020040 = False
id_1000260560 = False
id_1000100200 = False
id_1000060120 = False
id_1000040090 = False
id_1000070140 = False
id_1000130270 = False

for filename in infiles:
    print filename
    input = I3File(filename, 'r')
    
    while input.more():
        frame = input.pop_frame()
        if (frame.Stop == I3Frame.DAQ):
            continue
        elif (frame.Stop == I3Frame.Physics):

            # Apply the L5 BDT cut
            if not frame.Has("L4BDT"): continue
            if not frame["L4BDT"]["BDTScore"] > 0.04: continue
            
            # Lets just weight to H3a for now
            w = frame["GaisserH3aWeight"].value
            if "9036" in filename: w /= 100000.0
            elif "9255" in filename: w /= 95000.0
            elif "10282" in filename w /= 11200.0
            elif "10309" in filename w /= 12200.0
            elif "10369" in filename w /= 400.0
            elif "11905" in filename w /= 90000.0
            elif "12268" in filename w /= 96000.0
            elif "12332" in filename w /= 99000.0

            # Now, we want the primary. Let's get it
            primary = dataclasses.get_most_energetic_primary(frame["I3MCTree"])
            muon = dataclasses.get_most_energetic_muon(frame["I3MCTree"])


            # Propagate the muons muons to the MuonGun generating surface
            outSurface = Cylinder(1600*I3Units.m, 800*I3Units.m)
            partVec = MuonGun.muons_at_surface(frame, outSurface)

            muongun_energy = 1e-9
            for p in partVec:
                if p.energy > muongun_energy:
                    muongun_energy = p.energy                    

            weights.append(w)
            energy_of_primary.append(primary.energy)
            energy_at_depth.append(muongun_energy)
            coszens.append(numpy.cos(muon.dir.zenith))
	    code_of_primary.append(primary.pdg_encoding)
	    if not id_1000160320 and primary.pdg_encoding == 1000160320:
		print '1000160320 = %s'%primary.type
		id_1000160320 = True
	    if not id_1000080160 and primary.pdg_encoding == 1000080160:
		print '1000080160 = %s'%primary.type
		id_1000080160 = True
	    if not id_2212 and primary.pdg_encoding == 2212:
		print '2212 = %s'%primary.type
		id_2212 = True
	    if not id_1000020040 and primary.pdg_encoding == 1000020040:
		print '1000020040 = %s'%primary.type
		id_1000020040 = True
	    if not id_1000260560 and primary.pdg_encoding == 1000260560:
		print '1000260560 = %s'%primary.type
		id_1000260560 = True
	    if not id_1000100200 and primary.pdg_encoding == 1000100200:
	        print '1000100200 = %s'%primary.type
		id_1000100200 = True
	    if not id_1000060120 and primary.pdg_encoding == 1000060120:
	        print '1000060120 = %s'%primary.type
		id_1000060120 = True
	    if not id_1000040090 and primary.pdg_encoding == 1000040090:
		print '1000040090 = %s'%primary.type
		id_1000040090 = True
	    if not id_1000070140 and primary.pdg_encoding == 1000070140:
		print '1000070140 = %s'%primary.type
		id_1000070140 = True
	    if not id_1000130270 and primary.pdg_encoding == 1000130270:
		print '1000130270 = %s'%primary.type
		id_1000130270 = True
	    type_of_primary.append(primary.type)
                
    input.close()

weights = numpy.array(weights)
energy_of_primary = numpy.array(energy_of_primary)
energy_at_depth = numpy.array(energy_at_depth)
coszens = numpy.array(coszens)
code_of_primary = numpy.array(code_of_primary)
type_of_primary = numpy.array(type_of_primary)

results = {"weights":weights,
           "energy_of_primary": energy_of_primary,
           "energy_at_depth": energy_at_depth,
           "coszens": coszens,
	   "code_of_primary": code_of_primary,
	   "type_of_primary": type_of_primary,
           }

pickle.dump(results, open("corsika.pckl","w"))

print "Found %i events with total weight %4.2e Hz" % (len(weights), numpy.sum(weights))
