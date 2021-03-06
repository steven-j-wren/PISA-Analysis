#! /usr/bin/env python
#
# get_DH_slopes.py
#
# Get the slopes for DomEff and HoleIce fits.
#
# author: Feifei Huang - fxh140@psu.edu
#         Tim Arlen - tca3@psu.edu
#
# date:   02-July-2015
#
import matplotlib as mpl
mpl.use('Agg')
import copy
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pisa.utils.utils as utils
from pisa.utils.log import logging, tprofile, physics
from pisa.utils.jsons import from_json,to_json
from pisa.analysis.GetMCError import GetMCError
from pisa.analysis.TemplateMakerNoDiscSyst import TemplateMaker
from pisa.utils.params import get_values, select_hierarchy
from pisa.utils.plot import show_map
import os

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from scipy.optimize import curve_fit

parser = ArgumentParser(description='''Get the slopes for DOMeff and HoleIce fits. ''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-t','--template_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the template generation and systematics.''')
parser.add_argument('-s','--sim',type=str,
                    metavar='simu', required = True,
                    help='''Which simulation, can only be 4digit, 5digit, or dima''')
parser.add_argument('--templ_already_saved',action='store_true',default=False,
                    help="Read templates from already saved file; saves time when only need plotting.")
parser.add_argument('--plot',action='store_true',default=False,
                    help="Plot the fits of DOM efficiency and hole ice for each bin.")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in plot titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in plot titles")
parser.add_argument('-o','--outdir',type=str,
                    metavar='DIR', required = True, help='''Output directory''')
args = parser.parse_args()

#Read in the settings
detector = args.detector
selection = args.selection
outdir = args.outdir
utils.mkdir(outdir)
utils.mkdir(outdir+'/plots/')
template_settings = from_json(args.template_settings)
czbin_edges = template_settings['binning']['czbins']
ebin_edges = template_settings['binning']['ebins']
channel = template_settings['params']['channel']['value']
x_steps = 0.0001

if args.sim == '4digit':
    MC_name = '1XXX'
elif args.sim == '5digit':
    MC_name = '1XXXX'
elif args.sim == 'dima':
    MC_name = 'Dima'
else:
    MC_name = 'Other'

params = get_values(select_hierarchy(template_settings['params'],normal_hierarchy=True))
run_list = params['run_list']
run_nominal = params['run_nominal']
run_dict = params['run_dict']

for norm in [False]:

    if norm:
        HierarchyPrefix = 'NH'
        DMRange = np.linspace(template_settings['params']['deltam31_nh']['range'][0],template_settings['params']['deltam31_nh']['range'][1],21)
        THRange = np.linspace(0.35,1.25,10)
    else:
        HierarchyPrefix = 'IH'
        DMRange = np.linspace(template_settings['params']['deltam31_ih']['range'][0],template_settings['params']['deltam31_ih']['range'][1],21)
        THRange = np.linspace(0.35,1.25,10)

    for DMval in DMRange:

        DMPrefix = HierarchyPrefix + '_DM31_%.2f'%(DMval*1000.0)

        for THval in THRange:

            OutputFilePrefix = DMPrefix + '_TH23_%.2f'%THval

            templates = {}
            MCmaps = {}
            fits_DOMEff = {'trck':{'slopes':{}, 'fixed_ratios':{}},
                           'cscd':{'slopes':{}, 'fixed_ratios':{}},
                           'nominal_value': params['dom_eff']}
            fits_HoleIce = {'trck':{'slopes':{}, 'fixed_ratios':{}},
                            'cscd':{'slopes':{}, 'fixed_ratios':{}},
                            'nominal_value': params['hole_ice']}

            # Get templates and MC events 
            if not args.templ_already_saved:
                for run_num in run_list:
                    DH_template_settings = copy.deepcopy(template_settings)
                    MCmaps[str(run_num)] = {'trck':{}, 'cscd':{}}
                    templates[str(run_num)] = {'trck':{}, 'cscd':{}}
                    print "run_num = ", run_num
        
                    if args.sim == '5digit':
                        assert(DH_template_settings['pid_mode']['value']=='mc') # right now, only use MC mode for PID for the 5-digit sets 
                        aeff_mc_file = 'aeff/events__deepcore__IC86__runs_12%s-16%s:20000__proc_v5digit__unjoined.hdf5' % (run_num,run_num)
                        reco_mc_file = 'aeff/events__deepcore__IC86__runs_12%s-16%s:20000__proc_v5digit__joined_G_nue_cc+nuebar_cc_G_numu_cc+numubar_cc_G_nutau_cc+nutaubar_cc_G_nuall_nc+nuallbar_nc.hdf5' % (run_num, run_num)
            
                    elif args.sim == '4digit':
                        aeff_mc_file = run_dict[run_num]["aeff_file"]
                        reco_mc_file = run_dict[run_num]["reco_file"]
                        pid_param_file = run_dict[run_num]["pid_file"]
                        DH_template_settings['params']['pid_paramfile']['value'] = pid_param_file
            
                    else:
                        #TODO
                        print "to do, dima sets"
                    DH_template_settings['params']['aeff_weight_file']['value'] = aeff_mc_file
                    DH_template_settings['params']['reco_mc_wt_file']['value'] = reco_mc_file
                    if 'atmos_mu_scale' in DH_template_settings['params']:
                        DH_template_settings['params']['atmos_mu_scale']['value'] = 0.0

                    DH_params = select_hierarchy(DH_template_settings['params'],
                                                 normal_hierarchy=norm)
                    DH_params['deltam31']['value'] = DMval
                    DH_params['theta23']['value'] = THval
                    DH_template_maker = TemplateMaker(get_values(DH_params), **DH_template_settings['binning'])
                    template = DH_template_maker.get_template(get_values(DH_params))
    
                    templates[str(run_num)]['trck'] = template['trck']['map']
                    templates[str(run_num)]['cscd'] = template['cscd']['map']
                    templates[str(run_num)]['dom_eff'] = run_dict[str(run_num)]['dom_eff']
                    templates[str(run_num)]['hole_ice'] = run_dict[str(run_num)]['hole_ice']

                    # Get MC events map 
                    MCMap = GetMCError(get_values(DH_template_settings['params']), DH_template_settings['binning']['ebins'], DH_template_settings['binning']['czbins'], reco_mc_file)
                    tmap_MC = MCMap.get_mc_events_map(True, get_values(DH_template_settings['params']), reco_mc_file)
                    MCmaps[str(run_num)]['trck'] = tmap_MC['trck']['map']
                    MCmaps[str(run_num)]['cscd'] = tmap_MC['cscd']['map']
    
                #Assemble output dict
                output_template = {'templates' : templates,
                                   'MCmaps': MCmaps,
                                   'template_settings' : template_settings}
                to_json(output_template, outdir+'%s_%s_DomEff_HoleIce_templates_%i_by_%i.json'%(OutputFilePrefix, MC_name,len(ebin_edges)-1,len(czbin_edges)-1))

            # if templates already saved
            else:
                output_template = from_json(outdir+ '%s_DomEff_HoleIce_templates_10_by_16.json'%args.sim) 
                templates = output_template['templates']
                MCmaps = output_template['MCmaps']

            # Do fits (linear and quadratic) for each bin
            for flav in ['trck','cscd']:
                templ_list = []
                templ_err_list = []
                dom_eff_list = []
                hole_ice_list = []
                k_DE_linear = np.empty(np.shape(templates[run_nominal][flav])) 
                k_DE_quad = np.empty(np.shape(templates[run_nominal][flav])) 
                p_DE_quad = np.empty(np.shape(templates[run_nominal][flav])) 
                if args.sim == '4digit':
                    fixed_ratio = np.empty(np.shape(templates[run_nominal][flav])) 
                k_HI_quad = np.empty(np.shape(templates[run_nominal][flav])) 
                p_HI_quad = np.empty(np.shape(templates[run_nominal][flav])) 
                k_HI_linear = np.empty(np.shape(templates[run_nominal][flav])) 
                for run_num in run_list: 
                    dom_eff_list.append(templates[str(run_num)]['dom_eff'])
                    hole_ice_list.append(templates[str(run_num)]['hole_ice'])
                    templ_list.append(templates[str(run_num)][flav])
                    templ_err_list.append(templates[str(run_num)][flav]/np.sqrt(MCmaps[str(run_num)][flav]))
   
                dom_eff = np.array(dom_eff_list)
                hole_ice = np.array(hole_ice_list)      # unit: cm-1
                templ = np.array(templ_list)
                templ_err = np.array(templ_err_list)

                y_val_max = np.max(np.divide(templ, templates[run_nominal][flav]))
                y_val_min = np.min(np.divide(templ, templates[run_nominal][flav]))

                tml_shape = np.shape(templ)
                n_ebins = tml_shape[1] 
                n_czbins = tml_shape[2] 

############################# DOM efficiency ################################
    
                for i in range(0,n_ebins):
                    for j in range(0,n_czbins):
    
                        ########### Get Data ############
                        cut = hole_ice==0.02        # select elements when hole ice = 0.02
                        dom_eff_values = dom_eff[cut]
                        bin_counts = templ[cut,i,j]
                        bin_err = templ_err[cut,i,j]
                        nominal_bin_counts = bin_counts[dom_eff_values==1.0]
                        nominal_bin_err = bin_err[dom_eff_values==1.0]
                        #print "dom_eff_values = ", dom_eff_values
                        bin_ratio_values = bin_counts/nominal_bin_counts  #divide by the nominal value
                        bin_ratio_err_values = bin_ratio_values * np.sqrt(np.square(nominal_bin_err/nominal_bin_counts)+np.square(bin_err/bin_counts))

                        # line goes through point (0.02, fixed_r_val), fixed_r_val is the value for dom_eff = 0.91 and hole ice = 0.02
                        if args.sim == '4digit':
                            # line goes through point (0.02, fixed_r_val), fixed_r_val is the value for dom_eff = 0.91 and hole ice = 0.02
                            fixed_r_val = bin_ratio_values[dom_eff_values==0.91]
                            fixed_ratio[i][j]= fixed_r_val
                            exec('def dom_eff_linear_through_point(x, k): return k*x + %s - k*0.91'%fixed_r_val)
                            exec('def dom_eff_quadratic_through_point(x, k, p): return k*(x- 0.91) + p*(x- 0.91)**2 + %s'%fixed_r_val)
                        elif args.sim == '5digit':
                            exec('def dom_eff_linear_through_point(x, k): return k* (x - 1.0) + 1.0')
                            exec('def dom_eff_quadratic_through_point(x, k, p): return k*(x- 1.0) + p*(x-1.0)**2 + 1.0')
                        else:
                            #TODO
                            print "dima sets, to do"

                        ########### DOM efficiency #############

                        popt_1, pcov_1 = curve_fit(dom_eff_linear_through_point, dom_eff_values, bin_ratio_values)
                        k1 = popt_1[0]
                        k_DE_linear[i][j]= k1
   
                        popt_2, pcov_2 = curve_fit(dom_eff_quadratic_through_point, dom_eff_values, bin_ratio_values)
                        k2 = popt_2[0]
                        p2 = popt_2[1]
                        k_DE_quad[i][j]= k2
                        p_DE_quad[i][j]= p2

                        if args.plot:
                            fig_num = i * n_czbins+ j
                            if (fig_num == 0 or fig_num == n_czbins * n_ebins):
                                fig = plt.figure(num=1, figsize=( 4*n_czbins, 4*n_ebins))
                            subplot_idx = n_czbins*(n_ebins-1-i)+ j + 1
                            #print 'subplot_idx = ', subplot_idx
                            plt.subplot(n_ebins, n_czbins, subplot_idx)
                            #plt.title("CZ:[%s, %s] E:[%.1f, %.1f]"% (czbin_edges[j], czbin_edges[j+1], ebin_edges[i], ebin_edges[i+1]))
                            plt.scatter(dom_eff_values, bin_ratio_values, color='blue')
                            plt.errorbar(dom_eff_values, bin_ratio_values, yerr=bin_ratio_err_values,fmt='none')
                            plt.xlim(0.7,1.2)
                            plt.ylim(y_val_min-0.1,y_val_max+0.1)

                            dom_func_plot_x = np.arange(0.8 - x_steps, 1.2 + x_steps, x_steps)
                            dom_func_plot_y_linear = dom_eff_linear_through_point(dom_func_plot_x, k1)
                            dom_func_plot_linear, = plt.plot(dom_func_plot_x, dom_func_plot_y_linear, 'k-')
                            dom_func_plot_y_quad = dom_eff_quadratic_through_point(dom_func_plot_x, k2,p2)
                            dom_func_plot_quad, = plt.plot(dom_func_plot_x, dom_func_plot_y_quad, 'r-')
                            if j > 0:
                                plt.setp(plt.gca().get_yticklabels(), visible=False)
                            if i > 0:
                                plt.setp(plt.gca().get_xticklabels(), visible=False)
                            if(fig_num == n_czbins * n_ebins-1):
                                plt.figtext(0.5, 0.04, r'$\cos\theta_Z$',fontsize=60,ha='center')
                                plt.figtext(0.6, 0.06, r'x-axis - DOM Efficiency Value',fontsize=30,ha='left')
                                plt.figtext(0.6, 0.04, r'y-axis - Relative Bin Content (to nominal)',fontsize=30,ha='left') 
                                plt.figtext(0.04, 0.5, 'Energy [GeV]',rotation=90,fontsize=60,ha='center')
                                for i,enval in enumerate(ebin_edges):
                                    plt.figtext(0.085, 0.097+0.089*i, '%.2f'%enval,fontsize=25,ha='center')
                                for i,czval in enumerate(czbin_edges):
                                    plt.figtext(0.123+0.155*i, 0.080, '%.2f'%czval,fontsize=25,ha='center')
                                plt.figtext(0.5, 0.95, r'DOM Efficiency Slopes', fontsize=60,ha='center')
                                if flav == 'cscd':
                                    plt.figtext(0.5, 0.93, r'%s %s Cascade Channel'%(detector,selection), fontsize=60,ha='center')
                                elif flav == 'trck':
                                    plt.figtext(0.5, 0.93, r'%s %s Track Channel'%(detector,selection), fontsize=60,ha='center')
                                fig.subplots_adjust(hspace=0)
                                fig.subplots_adjust(wspace=0)
                                plt.savefig(outdir+ 'plots/'+'%s_%s_fits_domeff_%s.png'%(OutputFilePrefix, args.sim, flav))
                                plt.savefig(outdir+ 'plots/'+'%s_%s_fits_domeff_%s.pdf'%(OutputFilePrefix, args.sim, flav))
                                plt.clf()

    
                ########### Hole Ice #############
                
                for i in range(0,n_ebins):
                    for j in range(0,n_czbins):
    
                        ########### Get Data ############
                        if args.sim == '4digit':
                            cut = dom_eff== 0.91       # select elements when dom_eff = 0.91 
                        elif args.sim == '5digit':
                            cut = dom_eff== 1.0        # select elements when dom_eff = 1
                        else:
                            #TODO
                            print "Dima sets, to do"
                        hole_ice_values = hole_ice[cut]
                        bin_counts = templ[cut,i,j]
                        bin_err = templ_err[cut,i,j]
                        nominal_bin_counts = bin_counts[hole_ice_values==0.02]
                        nominal_bin_err = bin_err[hole_ice_values==0.02]
                        #print "hole_ice_values = ", hole_ice_values
                        bin_ratio_values = bin_counts/nominal_bin_counts  #divide by the nominal value
                        bin_ratio_err_values = bin_ratio_values * np.sqrt(np.square(nominal_bin_err/nominal_bin_counts)+np.square(bin_err/bin_counts))

                        exec('def hole_ice_linear_through_point(x, k): return k* (x - 0.02) + 1.0')
                        if args.sim == '4digit':
                            fixed_r_val = bin_ratio_values[hole_ice_values==0.02]
                            # line goes through point (0.02, fixed_r_val), fixed_r_val is the value for dom_eff = 0.91 and hole ice = 0.02
                            exec('def hole_ice_quadratic_through_point(x, k, p): return k*(x-0.02) + p*(x-0.02)**2 + %s'%fixed_r_val)

                        elif args.sim == '5digit':
                            exec('def hole_ice_quadratic_through_point(x, k, p): return k*(x-0.02) + p*(x-0.02)**2 + 1.0')
                        else:
                            #TODO
                            print "to do"

                        popt_1, pcov_1 = curve_fit(hole_ice_linear_through_point, hole_ice_values, bin_ratio_values)
                        k1 = popt_1[0]
                        k_HI_linear[i][j]= k1

                        popt_2, pcov_2 = curve_fit(hole_ice_quadratic_through_point, hole_ice_values, bin_ratio_values)
                        k2 = popt_2[0]
                        p2 = popt_2[1]
                        k_HI_quad[i][j]= k2
                        p_HI_quad[i][j]= p2

                        if args.plot:
                            fig_num = i * n_czbins+ j
                            if (fig_num == 0 or fig_num == n_czbins * n_ebins):
                                fig = plt.figure(num=2, figsize=( 4*n_czbins, 4*n_ebins))
                            subplot_idx = n_czbins*(n_ebins-1-i)+ j + 1
                            plt.subplot(n_ebins, n_czbins, subplot_idx)
                            #plt.title("CZ:[%s, %s] E:[%.1f, %.1f]"% (czbin_edges[j], czbin_edges[j+1], ebin_edges[i], ebin_edges[i+1]))
                            plt.scatter(hole_ice_values, bin_ratio_values, color='blue')
                            plt.errorbar(hole_ice_values, bin_ratio_values, yerr=bin_ratio_err_values,fmt='none')
                            plt.ylim(y_val_min-0.1,y_val_max+0.1)

                            if args.sim == '4digit':
                                ice_func_plot_x = np.arange(-0.02, 0.06 + x_steps, x_steps)
                                plt.xlim(-0.02,0.06)
                            elif args.sim == '5digit':
                                ice_func_plot_x = np.arange(0.005, 0.04 + x_steps, x_steps)
                                plt.xlim(0.005,0.04+x_steps)
                            else:
                                #TODO
                                print "dima sets, to do"
                            ice_func_plot_y_linear = hole_ice_linear_through_point(ice_func_plot_x, k1)
                            ice_func_plot_linear, = plt.plot(ice_func_plot_x, ice_func_plot_y_linear, 'k-')
                            ice_func_plot_y_quad = hole_ice_quadratic_through_point(ice_func_plot_x, k2,p2)
                            ice_func_plot_quad, = plt.plot(ice_func_plot_x, ice_func_plot_y_quad, 'r-')
                            if j > 0:
                                plt.setp(plt.gca().get_yticklabels(), visible=False)
                            if i > 0:
                                plt.setp(plt.gca().get_xticklabels(), visible=False)

                            if(fig_num==n_czbins * n_ebins-1):
                                fig.subplots_adjust(hspace=0)
                                fig.subplots_adjust(wspace=0)
                                plt.figtext(0.5, 0.04, r'$\cos\theta_Z$',fontsize=60,ha='center')
                                plt.figtext(0.6, 0.06, r'x-axis - Hole Ice Value',fontsize=30,ha='left')
                                plt.figtext(0.6, 0.04, r'y-axis - Relative Bin Content (to nominal)',fontsize=30,ha='left') 
                                plt.figtext(0.04, 0.5, 'Energy [GeV]',rotation=90,fontsize=60,ha='center')
                                for i,enval in enumerate(ebin_edges):
                                    plt.figtext(0.085, 0.097+0.089*i, '%.2f'%enval,fontsize=25,ha='center')
                                for i,czval in enumerate(czbin_edges):
                                    plt.figtext(0.123+0.155*i, 0.080, '%.2f'%czval,fontsize=25,ha='center')
                                plt.figtext(0.5, 0.95, r'Hole Ice Slopes', fontsize=60,ha='center')
                                if flav == 'cscd':
                                    plt.figtext(0.5, 0.93, r'%s %s Cascade Channel'%(detector,selection), fontsize=60,ha='center')
                                elif flav == 'trck':
                                    plt.figtext(0.5, 0.93, r'%s %s Track Channel'%(detector,selection), fontsize=60,ha='center')
                                plt.savefig(outdir+ 'plots/'+'%s_%s_fits_holeice_%s.png'%(OutputFilePrefix, args.sim, flav))
                                plt.savefig(outdir+ 'plots/'+'%s_%s_fits_holeice_%s.pdf'%(OutputFilePrefix, args.sim, flav))
                                plt.clf()

                fits_DOMEff[flav]['slopes'] = k_DE_linear
                fits_DOMEff[flav]['linear'] = k_DE_quad
                fits_DOMEff[flav]['quadratic'] = p_DE_quad

                fits_HoleIce[flav]['slopes'] = k_HI_linear
                fits_HoleIce[flav]['linear'] = k_HI_quad
                fits_HoleIce[flav]['quadratic'] = p_HI_quad

                if args.sim == '4digit':
                    fits_DOMEff[flav]['fixed_ratios'] = fixed_ratio 
                    fits_HoleIce[flav]['fixed_ratios'] = fixed_ratio

            #And write to file

            to_json(fits_DOMEff,outdir+'%s_%s_MC_DomEff_fits_%i_by_%i.json'%(OutputFilePrefix, MC_name,len(ebin_edges)-1,len(czbin_edges)-1))
            to_json(fits_HoleIce,outdir+'%s_%s_MC_HoleIce_fits_%i_by_%i.json'%(OutputFilePrefix, MC_name,len(ebin_edges)-1,len(czbin_edges)-1))
