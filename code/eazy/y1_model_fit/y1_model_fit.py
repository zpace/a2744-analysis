# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# make a bunch of spectra.param files (each with one template included), for templates with epoch less than 700My (.7Gy)
# change output directory specific to each template
# get chi2 for every fit (constrained within 99% bounds), plot each fit w/ transparency according to chi2 (chi2 has to be less than (1, 4, 9) above chi2 minimum)

# <codecell>

import numpy
import astropy.table as table
import matplotlib.pyplot as plt
import os
import re
import glob
from progressbar import ProgressBar
pbar = ProgressBar()

# <codecell>

#first, read in summary file, and figure out the redshift limits, minimum chi2, and all that jazz.
obj_summary = table.Table.read('../y1_K_newACS_BC03/OUTPUT_1-TEMPLATE/y1.zout', format = 'ascii')
u99 = obj_summary['u99'][0]
l99 = obj_summary['l99'][0]
minchi2 = obj_summary['chi_1'][0]
l99, u99, minchi2

# <codecell>

def read_template(fname):
    f = open(fname, 'r')
    lambda_conv = 1.0
    epoch = 0
    for i, line in enumerate(f):
        if 'Epoch (yr)' in line:
            match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
            num = [float(x) for x in re.findall(match_number, line)]
            epoch = float(num[0]*10**num[1]) / 1e9
    templatefile = ' '.join(('1', fname, str(lambda_conv), str(epoch), str(1.0)))
    return epoch, templatefile

# <codecell>

template_list = glob.glob('templates/BC03Templates/BC*')

# <codecell>

fig = plt.figure(figsize = (8, 16))
ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

ax1.semilogy()
ax1.set_title('photo-z fits (with K-band)')
ax1.set_xlabel('z')
ax1.set_ylabel('probability')

ax2.set_title('SED fits (K)')
ax2.set_xlabel('$\lambda (\AA)$')
ax2.set_ylabel('$F_{\\nu}$')

ax3.semilogy()
ax3.set_title('photo-z fits (no K-band)')
ax3.set_xlabel('z')
ax3.set_ylabel('probability')

ax4.set_title('SED fits (no K)')
ax4.set_xlabel('$\lambda (\AA)$')
ax4.set_ylabel('$F_{\\nu}$')

#plt.tight_layout()
#plt.show()
for template_name in pbar(template_list):
    epoch, templatefiletext = read_template(template_name)
    #print templatefiletext
    if epoch < 0.7:
        #print 'Running eazy for template', template_name
        #first make the OUTPUT directory, keeping in mind that it will later be renamed to reflect the actual template name
        
        #prepare to run eazy with K band detections
        if os.path.isdir('OUTPUT') == False: 
            os.mkdir('OUTPUT')
        paramfile = open('templates/BC03Templates.spectra.param', 'a')
        paramfile.truncate(0)
        paramfile.write(templatefiletext)
        paramfile.close()
    
        #now go ahead and run eazy, then rename K-band OUTPUT folder to OUTPUT-<template>
        new_output_dir = 'OUTPUT-' + template_name[24:]
        if os.path.isdir('K/' + new_output_dir) == False: #this saves some time by only regenerating as necessary
            sys_cmd = './eazy -p zphot.param.y1 -z zphot.zeropoint.a2744'
            os.system(sys_cmd)
            #os.system('pwd')
            new_output_dir = 'OUTPUT-' + template_name[24:]
            os.rename('OUTPUT', new_output_dir)
            os.system('mv ' + new_output_dir + ' K/' + new_output_dir)
            
        # =====
            
        #prepare to run eazy **without** K band detections
        if os.path.isdir('OUTPUT') == False: 
            os.mkdir('OUTPUT')
        paramfile = open('templates/BC03Templates.spectra.param', 'a')
        paramfile.truncate(0)
        paramfile.write(templatefiletext)
        paramfile.close()
        
        wait = raw_input('Change the zphot.translate file and hit ENTER to continue')
    
        #now go ahead and run eazy, then rename no K-band OUTPUT folder to no-K/OUTPUT-<template>
        new_output_dir = 'OUTPUT-' + template_name[24:]
        if os.path.isdir('no-K/' + new_output_dir) == False: #this saves some time by only regenerating as necessary
            print new_output_dir
            sys_cmd = './eazy -p zphot.param.y1 -z zphot.zeropoint.a2744'
            os.system(sys_cmd)
            #os.system('pwd')
            new_output_dir = 'OUTPUT-' + template_name[24:]
            os.rename('OUTPUT', new_output_dir)
            os.system('mv ' + new_output_dir + ' no-K/' + new_output_dir)
            
        # =====
            
        #now test whether the K-band run has a chi2 below the threshold established above
        if os.stat('K/' + new_output_dir + '/y1.zout')[6] != 0:
            template_fit_summary = table.Table.read('K/' + new_output_dir + '/y1.zout', format = 'ascii')
            template_fit_chi2 = template_fit_summary['chi_1'][0]
            z_peak = template_fit_summary['z_peak'][0]
            #print template_fit_chi2
        
        print new_output_dir
        
        if os.stat('K/' + new_output_dir + '/y1.zout')[6] != 0:
            #now test whether the no-K-band run has a chi2 below the threshold established above
            template_fit_summary_noK = table.Table.read('no-K/' + new_output_dir + '/y1.zout', format = 'ascii')
            template_fit_chi2_noK = template_fit_summary['chi_1'][0]
            z_peak_noK = template_fit_summary['z_peak'][0]
            #print template_fit_chi2
        
        # =====
        
        sigma = 2
        maxchi2 = minchi2 + minchi2 * sigma**2.
        
        #now plot the fits
        if (template_fit_chi2 <= maxchi2) and (os.path.isfile('K/' + new_output_dir + '/1.pz') == True) and (os.path.isfile('K/' + new_output_dir + '/1.temp_sed') == True):
            #print 'plotting!'
            
            #first do all the plots with the K band data
            pz = table.Table.read('K/' + new_output_dir + '/1.pz', format = 'ascii')
            z = pz['z']
            prob = pz['pz']/np.trapz(pz['pz'], z)
            ax1.plot(z, prob, color = 'b', alpha = 0.5 * (template_fit_chi2/maxchi2), linewidth = 0.25)
            
            temp_sed = table.Table.read('K/' + new_output_dir + '/1.temp_sed', format = 'ascii')
            f_n = temp_sed['tempflux']
            f_l = f_n * (l/5000.)**2.
            ax2.plot(l[l < 46000.], f_l[l < 46000.], color = 'r', alpha = 0.5 * (template_fit_chi2/maxchi2), linewidth = 0.25)
        
        if (template_fit_chi2 <= maxchi2) and (os.path.isfile('no-K/' + new_output_dir + '/1.pz') == True) and (os.path.isfile('no-K/' + new_output_dir + '/1.temp_sed') == True):
            #now do all the plots with no K-band data
            pz = table.Table.read('no-K/' + new_output_dir + '/1.pz', format = 'ascii')
            z = pz['z']
            prob = pz['pz']/np.trapz(pz['pz'], z)
            
            temp_sed = table.Table.read('no-K/' + new_output_dir + '/1.temp_sed', format = 'ascii')
            ax3.plot(z, prob, color = 'b', alpha = 0.5 * (template_fit_chi2_noK/maxchi2), linewidth = 0.25)
            f_n = temp_sed['tempflux']
            f_l = f_n * (l/5000.)**2.
            ax4.plot(l[l < 46000.], f_l[l < 46000.], color = 'r', alpha = 0.5 * (template_fit_chi2_noK/maxchi2), linewidth = 0.25)
            
plt.show()

# <codecell>


