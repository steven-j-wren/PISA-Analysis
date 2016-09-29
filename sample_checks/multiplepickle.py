#!/usr/bin/env python
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True
from icecube import dataclasses, icetray
from icecube.dataio import I3File
from icecube.icetray import I3Frame, I3Units
from icecube import ocelot
from icecube import NuFlux

from msu import passed_msu
from nbi import passed_nbi
from prd import passed_prd

# Do it for one file just to see if it works

def read_in_files(indirsdict):

    filenames = {}

    for first_selection in indirsdict.keys():

        first_indirs = indirsdict[first_selection]
        filenames[first_selection] = []

        for first_indir in first_indirs:

            first_infiles = glob.glob(first_indir + "*")
            first_infiles.sort()
            if first_selection in ['prd','msu']:
                first_infiles = first_infiles
                first_filenums = [x.split('/')[-1].split('.')[-3] for x in first_infiles]
            elif first_selection in ['nbi']:
                first_infiles = first_infiles
                first_filenums = [x.split('/')[-1].split('.')[-4] for x in first_infiles]

            # Initialise a boolean array of the right length. All true to start.
            matching_filenums = numpy.in1d(first_filenums, first_filenums)

            for second_selection in indirsdict.keys():

                second_indirs = indirsdict[second_selection]
            
                for second_indir in second_indirs:
                    
                    second_infiles = glob.glob(second_indir+"*")
                    second_infiles.sort()
                    second_sizes = [os.path.getsize(x) for x in second_infiles]
                    second_sizes_mean = sum(second_sizes)/len(second_sizes)
                    second_sizes_keep = numpy.array([x > second_sizes_mean/5 for x in second_sizes])
                    if second_selection in ['prd','msu']:
                        second_filenums = [x.split('/')[-1].split('.')[-3] for x in second_infiles]
                    elif second_selection in ['nbi']:
                        second_filenums = [x.split('/')[-1].split('.')[-4] for x in second_infiles]

                    second_filenums = numpy.array(second_filenums)
                    this_matching_filenums = numpy.in1d(first_filenums,second_filenums[second_sizes_keep])
                    matching_filenums = numpy.logical_and(matching_filenums, this_matching_filenums)
            matching_files = numpy.array(first_infiles)[matching_filenums]

            first_sizes = [os.path.getsize(x) for x in matching_files]
            first_sizes_mean = sum(first_sizes)/len(first_sizes)
            first_sizes_keep = numpy.array([x > first_sizes_mean/5 for x in first_sizes])
            kept_files = matching_files[first_sizes_keep]
            for kept_file in kept_files:
                filenames[first_selection].append(kept_file)

    return filenames
                
        
def extract_run_num(filename, selection):

    if selection in ['msu']:
        run_num = int(os.path.basename(filename).split(".")[2])
    elif selection in ['nbi','prd']:
        run_num = int(os.path.basename(filename).split(".")[1])
    else:
        raise ValueError("I didn't understand that selection choice - %s"%selection)

    return run_num

                    
if __name__ == '__main__':

    Honda2014SPLFluxIP = NuFlux.makeFlux('IPhonda2014_spl_solmin')

    indirs_dict = {}
    indirs_dict['msu'] = ["/data/user/hignight/nue/12600/",
                          "/data/user/hignight/numu/14600/",
                          "/data/user/hignight/nutau/16600/"]
    indirs_dict['nbi'] = ["/data/user/mlarson/tau_appearance/scripts/level7/output/12600/",
                          "/data/user/mlarson/tau_appearance/scripts/level7/output/14600/",
                          "/data/user/mlarson/tau_appearance/scripts/level7/output/16600/"]
    indirs_dict['prd'] = ["/data/user/terliuk/nu_L5/12600/",
                          "/data/user/terliuk/nu_L5/14600/",
                          "/data/user/terliuk/nu_L5/16600/"]

    filenames = read_in_files(indirs_dict)

    eventids = {}
    unosc_weight = {}
    osc_weight = {}
    energy = {}
    coszen = {}
    isnu = {}
    reco_energy = {}
    reco_coszen = {}

    cases = ['nue_cc', 'numu_cc', 'nutau_cc', 'nuall_nc']
    requirements = {}
    requirements['nue_cc'] = {}
    requirements['nue_cc']['run_num'] = 12600
    requirements['numu_cc'] = {}
    requirements['numu_cc']['run_num'] = 14600
    requirements['nutau_cc'] = {}
    requirements['nutau_cc']['run_num'] = 16600
    requirements['nuall_nc'] = {}
    requirements['nuall_nc']['run_num'] = [12600, 14600, 16600]

    for selection in filenames.keys():
        
        eventids[selection] = {}
        unosc_weight[selection] = {}
        osc_weight[selection] = {}
        energy[selection] = {}
        coszen[selection] = {}
        reco_energy[selection] = {}
        reco_coszen[selection] = {}
        isnu[selection] = {}

        for case in cases:
            eventids[selection][case] = []
            unosc_weight[selection][case] = []
            osc_weight[selection][case] = []
            energy[selection][case] = []
            coszen[selection][case] = []
            reco_energy[selection][case] = []
            reco_coszen[selection][case] = []
            isnu[selection][case] = []

        for filename in filenames[selection]:

            print filename
            input = I3File(filename, 'r')
            run_num = extract_run_num(filename = filename,
                                      selection = selection)
        
            while input.more():

                try:
                    frame = input.pop_physics()
                except:
                    continue

                if selection == 'msu':                 
                    if not passed_msu(frame): continue
                    
                    cascfit = frame['IC86_Dunkman_L6_PegLeg_MultiNest8D_HDCasc']
                    trackfit = frame["IC86_Dunkman_L6_PegLeg_MultiNest8D_Track"]
                    just_file = filename.split('/')[-1]
                    just_name, ext = os.path.splitext(os.path.splitext(just_file)[0])
                    file_num = str(run_num)+'.'+just_name.split('.')[-1]
                elif selection == 'nbi':
                    if not passed_nbi(frame): continue
                    
                    cascfit = frame["Pegleg_Fit_MNHDCasc"]
                    trackfit = frame["Pegleg_Fit_MNTrack"]
                    just_file = filename.split('/')[-1]
                    just_name, ext = os.path.splitext(os.path.splitext(just_file)[0])
                    file_num = str(run_num)+'.'+just_name.split('.')[-2]
                    
                elif selection == 'prd':
                    if not passed_prd(frame): continue
                    
                    cascfit = frame["SANTA_Fit_CascadeHad"]
                    trackfit = frame["SANTA_Fit_Muon"]
                    just_file = filename.split('/')[-1]
                    just_name, ext = os.path.splitext(os.path.splitext(just_file)[0])
                    file_num = str(run_num)+'.'+just_name.split('.')[-1]
                else:
                    continue

                ocelot.Ocelot(frame,
                              nue_flux_service = Honda2014SPLFluxIP,
                              numu_flux_service = Honda2014SPLFluxIP,
                              delta_m21_squared = 7.6e-5,             #< eV^2                
                              delta_m31_squared = 2.426e-3,           #< eV^2                 
                              theta_12 = math.asin(math.sqrt(0.312)), #< rad                    
                              theta_23 = math.asin(math.sqrt(0.42)),  #< rad
                              theta_13 = math.asin(math.sqrt(0.025)), #< rad
                              delta_cp = 0.,                          #< rad
                              cache_base = "Ocelot",
                              output_name = "DummyNeutrinoWeights")

                osc_w = frame["DummyNeutrinoWeights"]["OscillatedRate"] / float(len(filenames['msu'])/3.0)
                unosc_w = frame["DummyNeutrinoWeights"]["UnoscillatedRate"] / float(len(filenames['msu'])/3.0)

                # fix the nu/nubar junk from ocelot, which is just wrong
                p = dataclasses.get_most_energetic_neutrino(frame['I3MCTree'])
                if p.pdg_encoding > 0:
                    osc_w *= 0.7/2.0
                    unosc_w *= 0.7/2.0
                else:
                    osc_w *= 0.3/2.0
                    unosc_w *= 0.3/2.0

                eventid = file_num + str(frame["I3EventHeader"].event_id)
                uniqueid = str(eventid)+'_'+str(p.energy)+'_'+str(numpy.cos(p.dir.zenith))
            
                isCC = frame["I3MCWeightDict"]["InteractionType"]==1
                isneutrino = p.pdg_encoding > 0.0

                if isCC:
                    for case in requirements.keys():
                        if requirements[case]['run_num'] == run_num:
                            unosc_weight[selection][case].append(unosc_w)
                            osc_weight[selection][case].append(osc_w)
                            energy[selection][case].append(p.energy)
                            coszen[selection][case].append(numpy.cos(p.dir.zenith))
                            eventids[selection][case].append(uniqueid)        
                            reco_energy[selection][case].append(trackfit.energy + cascfit.energy)
                            reco_coszen[selection][case].append(numpy.cos(trackfit.dir.zenith))
                            isnu[selection][case].append(isneutrino)
                                
                else:
                    unosc_weight[selection]['nuall_nc'].append(unosc_w)
                    osc_weight[selection]['nuall_nc'].append(osc_w)
                    energy[selection]['nuall_nc'].append(p.energy)
                    coszen[selection]['nuall_nc'].append(numpy.cos(p.dir.zenith))
                    eventids[selection]['nuall_nc'].append(uniqueid)     
                    reco_energy[selection]['nuall_nc'].append(trackfit.energy + cascfit.energy)
                    reco_coszen[selection]['nuall_nc'].append(numpy.cos(trackfit.dir.zenith))
                    isnu[selection]['nuall_nc'].append(isneutrino)

            input.close()

        for case in cases:
        
            eventids[selection][case] = numpy.array(eventids[selection][case])
            unosc_weight[selection][case] = numpy.array(unosc_weight[selection][case])
            osc_weight[selection][case] = numpy.array(osc_weight[selection][case])
            energy[selection][case] = numpy.array(energy[selection][case])
            coszen[selection][case] = numpy.array(coszen[selection][case])    
            reco_energy[selection][case] = numpy.array(reco_energy[selection][case])
            reco_coszen[selection][case] = numpy.array(reco_coszen[selection][case])
            isnu[selection][case] = numpy.array(isnu[selection][case])

    results = {}
    results['eventids'] = eventids
    results['unosc_weight'] = unosc_weight
    results['osc_weight'] = osc_weight
    results['energy'] = energy
    results['coszen'] = coszen
    results['isnu'] = isnu
    results['reco_energy'] = reco_energy
    results['reco_coszen'] = reco_coszen

    pickle.dump(results, open("samplecomparisons.pckl","w"))
