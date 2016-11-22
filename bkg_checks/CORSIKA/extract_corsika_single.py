#!/usr/bin/env python
import os, sys, pickle, numpy

from icecube import dataclasses, icetray, MuonGun
from icecube.icetray import I3Frame, I3Units 
from icecube.dataio import I3File
from icecube.MuonGun import Cylinder

infiles = ["Level5_corsika.009036.100000files.i3",
           "Level5_corsika.009255.95000files.i3",
           "Level5_corsika.010282.11200files.i3",
           "Level5_corsika.010309.12200files.i3",
           "Level5_corsika.010369.400files.i3",
           "Level5_corsika.011905.90000files.i3",
           "Level5_corsika.012268.96000files.i3",
           "Level5_corsika.012332.99000files.i3"]

id_1000160320 = False
id_1000080160 = False
id_2212 = False
id_1000100200 = False
id_1000260560 = False
id_1000020040 = False
id_1000060120 = False
id_1000040090 = False
id_1000070140 = False
id_1000130270 = False
id_1000220480 = False
id_1000120240 = False
id_1000030070 = False
id_1000190390 = False
id_1000140280 = False
id_1000250550 = False
id_1000050110 = False
id_1000110230 = False
id_1000200400 = False
id_1000090190 = False
id_1000240520 = False
id_1000180400 = False
id_1000150310 = False
id_1000230510 = False

for filename in infiles:
	
    weights = []
    energy_of_primary = []
    energy_at_depth = []
    coszens = []
    code_of_primary = []
    type_of_primary = []

    i=0

    input = I3File(filename, 'r')
    
    while input.more():
        frame = input.pop_frame()
        if (frame.Stop == I3Frame.DAQ):
            continue
        elif (frame.Stop == I3Frame.Physics):

            if not frame["IC2013_LE_L3"].value: continue

            # Apply the L5 BDT cut
            if "9036" in filename: 
                BDT_Name = "L3BDT"
                Next_Name = "L4BDT"
            elif "9255" in filename:
                BDT_Name = "L3BDT"
                Next_Name = "L4BDT"
            else:
                BDT_Name = "L4BDT"
                Next_Name = "L5BDT"
            if not frame.Has(BDT_Name): continue
            if not frame[BDT_Name]["BDTScore"] > 0.04: continue
            if not frame.Has(Next_Name): continue
            
            # Lets just weight to H3a for now
            w = frame["GaisserH3aWeight"].value
            if numpy.isnan(w) or numpy.isinf(w): continue
            if "9036" in filename:
                w /= 100000.0
                name = "9036"
                if i==0:
                    print frame["CorsikaWeightMap"]["NEvents"]
                    i=1
            elif "9255" in filename:
                w /= 95000.0
                name = "9255"
            elif "10282" in filename:
                w /= 11200.0
                name = "10282"
            elif "10309" in filename:
                w /= 12200.0
                name = "10309"
            elif "10369" in filename:
                w /= 400.0
                name = "10369"
            elif "11905" in filename:
                w /= 90000.0
                name = "11905"
            elif "12268" in filename:
                w /= 96000.0
                name = "12268"
            elif "12332" in filename:
                w /= 99000.0
                name = "12332"
                if i==0:
                    print frame["CorsikaWeightMap"]["NEvents"]
                    i=1

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
            if not id_1000100200 and primary.pdg_encoding == 1000100200:
                print '1000100200 = %s'%primary.type
                id_1000100200 = True
            if not id_1000260560 and primary.pdg_encoding == 1000260560:
                print '1000260560 = %s'%primary.type
                id_1000260560 = True
            if not id_1000020040 and primary.pdg_encoding == 1000020040:
                print '1000020040 = %s'%primary.type
                id_1000020040 = True
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
            if not id_1000220480 and primary.pdg_encoding == 1000220480:
                print '1000220480 = %s'%primary.type
                id_1000220480 = True
            if not id_1000120240 and primary.pdg_encoding == 1000120240:
                print '1000120240 = %s'%primary.type
                id_1000120240 = True
            if not id_1000030070 and primary.pdg_encoding == 1000030070:
                print '1000030070 = %s'%primary.type
                id_1000030070 = True
            if not id_1000190390 and primary.pdg_encoding == 1000190390:
                print '1000190390 = %s'%primary.type
                id_1000190390 = True
            if not id_1000140280 and primary.pdg_encoding == 1000140280:
                print '1000140280 = %s'%primary.type
                id_1000140280 = True
            if not id_1000250550 and primary.pdg_encoding == 1000250550:
                print '1000250550 = %s'%primary.type
                id_1000250550 = True
            if not id_1000050110 and primary.pdg_encoding == 1000050110:
                print '1000050110 = %s'%primary.type
                id_1000050110 = True
            if not id_1000110230 and primary.pdg_encoding == 1000110230:
                print '1000110230 = %s'%primary.type
                id_1000110230 = True
            if not id_1000200400 and primary.pdg_encoding == 1000200400:
                print '1000200400 = %s'%primary.type
                id_1000200400 = True
            if not id_1000090190 and primary.pdg_encoding == 1000090190:
                print '1000090190 = %s'%primary.type
                id_1000090190 = True
            if not id_1000240520 and primary.pdg_encoding == 1000240520:
                print '1000240520 = %s'%primary.type
                id_1000240520 = True
            if not id_1000180400 and primary.pdg_encoding == 1000180400:
                print '1000180400 = %s'%primary.type
                id_1000180400 = True
            if not id_1000150310 and primary.pdg_encoding == 1000150310:
                print '1000150310 = %s'%primary.type
                id_1000150310 = True
            if not id_1000230510 and primary.pdg_encoding == 1000230510:
                print '1000230510 = %s'%primary.type
                id_1000230510 = True
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
               "type_of_primary": type_of_primary}

    pickle.dump(results, open("%scorsika.pckl"%name,"w"))

    print "Found %i events with total weight %4.2e Hz in file %s" % (len(weights), numpy.sum(weights), name)
