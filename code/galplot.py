import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from astropy import table
from astropy.io import fits
from progressbar import ProgressBar
import os

'''
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
'''

bands = {'105': 1055.2, '125': 1248.6, '140': 1392.3, '160': 1536.9, '435': 429.7, '606': 590.7, '814': 833.3, 'K': 2200}
bands_wid = {'105': 265/2., '125': 284.5/2., '140': 384/2., '160': 268.3/2., '435': 103.8/2., '606': 234.2/2., '814': 251.1/2., 'K': 400./2.}

def get_z(dir, row):
	'''
	given an integer row number for a particular galaxy, return its z, its image coordinates, its SED array and its probability distribution in z
	'''
	row = str(row)
	#read in its output table produced by eazy and get its 
	f = table.Table.read(dir + '/OUTPUT/' + row + '.pz', format = 'ascii')
	z = f['z']
	logchi2 = np.exp(-f['chi2']/2.)
	cumul = np.trapz(logchi2, z)
	prob = (1./cumul) * logchi2
	
	prob_prior = f['pz']
	cumul2 = np.trapz(prob_prior, z)
	prob_prior *= (1./cumul2)
	
	#now get the SED data for the object from the full catalog
	obj = table.Table.read(dir+ '/A2744_cat.dat', format = 'ascii')[int(row)-1]
	#print obj
	x_image = obj['X_IMAGE']
	y_image = obj['Y_IMAGE']
	
	F160Jy = obj['FLUX_AUTO_160']
	F160Mag = np.round(23.9 - 2.5*np.log10(F160Jy), 2)
	
	fluxnames = [i for i in obj.colnames if 'FLUX_APER' in i]
	'''
	#the SED will have a col for the band wavelength, a col for the half-width of the band, a col for the flux in that band, and a col for the flux error in that band
	SED = []
	
	for item in fluxnames:
		l = bands[item[10:]]
		lw = bands_wid[item[10:]]
		f = obj[item]
		df = obj['FLUXERR_APER_' + item[10:]]
		SED.append([l, lw, f, df])
	
	SED = np.asarray(SED)
	'''
	#finally, read in and return the fit SEDs
	obs_SED = table.Table.read(dir + '/OUTPUT/' + row + '.obs_sed', format = 'ascii')['lambda', 'flux_cat', 'err_cat', 'tempa_z']
	temp_SED = table.Table.read(dir + '/OUTPUT/' + row + '.temp_sed', format = 'ascii')['lambda', 'tempflux']
	
	return z, prob, prob_prior, x_image, y_image, obs_SED, temp_SED, F160Mag

def galplot(dir, row, save = False, plot_dir = 'figs'):	
	'''
	wraps around get_z(), and plots the following panels: 
	SED; redshift dist. (computed from chi2); image cutout in 105, 125, 140, 160, 435, 606, 814, K
	'''
	plt.close('all')
	if save == False:
		print 'Object #' + str(row) + '...' 
	z, prob, prob_prior, x_image, y_image, obs_SED, temp_SED, F160Mag = get_z(dir, row)
	zpeak = table.Table.read(dir + '/OUTPUT/eazy-cat.zout', format = 'ascii')['z_peak'][int(row)-1]
	
	fig = plt.figure(figsize = (8, 12))
	
	ax1 = plt.subplot2grid((4,2), (0,0), colspan = 2) #F_lambda SED
	ax2 = plt.subplot2grid((4,2), (1,0)) #z prob dist
	ax3 = plt.subplot2grid((4,2), (1,1)) #F_nu SED
	ax5 = plt.subplot2grid((4,2), (2,0), rowspan = 2, colspan = 2) #images
	
	len_a = np.sum(temp_SED['lambda'] < 2.5e4)
	F_l_obs = obs_SED['flux_cat'] * ((obs_SED['lambda']/5500.)**-2.)
	F_l_obs_e = F_l_obs - ((obs_SED['flux_cat'] + obs_SED['err_cat']) * ((obs_SED['lambda']/5500.)**-2.))
	
	F_l_temp = temp_SED['tempflux'] * ((temp_SED['lambda']/5500.)**-2.)
	
	ax1.errorbar(obs_SED['lambda'], F_l_obs, yerr = F_l_obs_e, marker = 'x', linestyle = 'None', color = 'b', label = 'Obs. SED')
	ax1.scatter(obs_SED['lambda'], obs_SED['tempa_z'] * ((obs_SED['lambda']/5500.)**-2.), marker = 's', s = 40, color = 'g', label = 'Filt. Flux', alpha = 0.3)
	ax1.plot(temp_SED['lambda'][:len_a], F_l_temp[:len_a], color = 'r', label = 'SED')
	ax1.legend(loc = 'best', prop = {'size': 8})	
	ax1.set_xlim([0., 25000.])
	
	ax1.set_xlabel('wavelength ($\AA$)')
	ax1.set_ylabel('Flux ($F_{\lambda}$)')
	
	prob_tol = .0001
	#find lowest & highest integer values of z s.t. prob_prior and prob are > prob_tol
	plotmin = np.floor( np.min( z[ (prob > prob_tol) * (prob_prior > prob_tol)]) )
	plotmax = np.ceil( np.max( z[ (prob > prob_tol) * (prob_prior > prob_tol)]) )
	#print plotmin, plotmax
	
	ax2.plot(z, prob, label = 'Raw prob.')
	ax2.plot(z, prob_prior, linestyle = '--', label = 'Prob. with prior')
	ax2.text(plotmin + .1*(plotmax-plotmin), 0.55 * np.max(prob), ' $z_{peak}$ = ' + str(np.round(zpeak, 3)))
	ax2.legend(loc = 'best', prop = {'size': 8})
	
	ax2.set_ylabel('Prob. dist.')
	ax2.set_xlabel('z')
	ax2.set_xlim([np.max([plotmin, 0.]), plotmax])
	
	ax3.errorbar(obs_SED['lambda'], obs_SED['flux_cat'], yerr = obs_SED['err_cat'], marker = 'x', linestyle = 'None', color = 'b', label = 'Obs. SED')
	ax3.scatter(obs_SED['lambda'], obs_SED['tempa_z'], marker = 's', s = 40, color = 'g', label = 'Filt. Flux', alpha = 0.3)
	ax3.plot(temp_SED['lambda'][:len_a], temp_SED['tempflux'][:len_a], color = 'r', label = 'SED')
	ax3.legend(loc = 'best', prop = {'size': 8})
	ax3.set_xlabel('wavelength ($\AA$)')
	ax3.set_ylabel('Flux ($F_{\\nu}$)')
	ax3.set_xlim([0., 25000.])
	
	if save == False: 
		print 'Getting images of object at', str(y_image) + ',', str(x_image), '...'
	
	im105 = fits.open('105.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	im125 = fits.open('125.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]	
	im140 = fits.open('140.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	im160 = fits.open('160.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	im435 = fits.open('435.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	im606 = fits.open('606.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	im814 = fits.open('814.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	imK = fits.open('K.fits')[0].data[y_image - 50:y_image + 50, x_image - 50:x_image + 50]
	
	im1 = np.hstack((im105, im125, im140, im160))
	im2 = np.hstack((im435, im606, im814, imK))
	
	#im3 = np.hstack((im105, im125))
	#im4 = np.hstack((im140, im160))
	#im5 = np.hstack((im435, im606))
	#im6 = np.hstack((im814, imK))
	
	im = np.vstack((im2, im1)) # this ordering puts im1 on the top, strangely...
	
	ax5.imshow(im + np.abs(np.min(im)) + .01, cmap = 'gray', norm = LogNorm(), interpolation = 'None') # I futzed with the logarithmic color mapping
	#add some identifying text
	textpos = [ [5, 105], [105, 105], [205, 105], [305, 105], [5, 5], [105, 5], [205, 5], [305, 5] ]
	texttext = ['F105W', 'F125W', 'F140W', 'F160W', 'F435', 'F606', 'F814', 'K-band']
	for i in range(len(texttext)):
		#print textpos[i]
		ax5.text(x = textpos[i][0], y = textpos[i][1], s = texttext[i], color = 'w')
	plt.axis('off')
	
	#now add flux circles of r = .5 arcseconds at the center of each image (couldn't figure out how to make things draw automatically, so here it is the old-fashioned way)
	imsize_y, imsize_x = np.shape(im)
	num_x = 4
	num_y = 2
	r = .5 / .060 + .060
	theta = np.linspace(0, 2.*np.pi, 100)
	circ_ctrs_x = np.linspace(0.5 * (imsize_x/num_x), imsize_x - 0.5 * (imsize_x/num_x), num = num_x, endpoint = True)
	circ_ctrs_y = np.linspace(0.5 * (imsize_y/num_y), imsize_y - 0.5 * (imsize_y/num_y), num = num_y, endpoint = True)
	#print circ_ctrs_x
	#print circ_ctrs_y
	for x in circ_ctrs_x:
		for y in circ_ctrs_y:
			ax5.plot(x + r*np.cos(theta), y + r*np.sin(theta), c = 'b')

	plt.suptitle('A2744 catalog - object ' + str(row) + ': F160W $M_{AB}$ = ' + str(F160Mag))
	plt.tight_layout()
	plt.subplots_adjust(top=0.95)
	if save == True:
		if not os.path.exists(dir + '/' + plot_dir): os.makedirs(dir + '/' + plot_dir)
		plt.savefig(dir + '/' + plot_dir + '/' + str(row) + '.png')
	else:
		plt.show()
		
pbar = ProgressBar() #don't move this

def find_highz(dir, zoutname, catname, zllim, zulim):
	print 'Filtering catalog...'
	f = table.Table.read(dir + '/' + zoutname, format = 'ascii')
	
	cat = table.Table.read(dir + '/' + catname, format = 'ascii')['id', 'X_IMAGE', 'Y_IMAGE']
	wht_105 = fits.open('105_wht.fits')[0].data
	wht_125 = fits.open('125_wht.fits')[0].data
	wht_140 = fits.open('140_wht.fits')[0].data
	wht_160 = fits.open('160_wht.fits')[0].data
	wht_435 = fits.open('435_wht.fits')[0].data
	wht_606 = fits.open('606_wht.fits')[0].data
	wht_814 = fits.open('814_wht.fits')[0].data
	
	fracmax = 0.1
	
	max_wht_105 = fracmax * np.max(wht_105.flatten())
	max_wht_125 = fracmax * np.max(wht_125.flatten())
	max_wht_140 = fracmax * np.max(wht_140.flatten())
	max_wht_160 = fracmax * np.max(wht_160.flatten())
	max_wht_435 = fracmax * np.max(wht_435.flatten())
	max_wht_606 = fracmax * np.max(wht_606.flatten())
	max_wht_814 = fracmax * np.max(wht_814.flatten())
	
	#print max_wht_105, max_wht_125, max_wht_140, max_wht_160, max_wht_435, max_wht_606, max_wht_814 
	
	f['wht_err'] = 0
	#print len(f['id']), len(cat['id'])
	#now test whether weights are good
	for obj in cat['id']:
		#print 'Testing object #' + str(obj)
		x_image = cat['X_IMAGE'][obj-1]
		y_image = cat['Y_IMAGE'][obj-1]
		#test each filter's wht file individually, except for K-band
		if np.any(wht_105[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_105):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '105 bad'
		if np.any(wht_125[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_125):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '125 bad'
		if np.any(wht_140[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_140):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '140 bad'
		if np.any(wht_160[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_160):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '160 bad'
		if np.any(wht_435[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_435):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '435 bad'
		if np.any(wht_606[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_606):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '606 bad'
		if np.any(wht_814[y_image - 9:y_image + 9, x_image - 9:x_image + 9] < max_wht_814):
			f['wht_err'][obj-1] += 1
			if obj == 302: print '814 bad'
	#print f['wht_err']
	'''
	plt.hist(f['wht_err'], bins = np.max(f['wht_err']))
	plt.title('fracmax: ' + str(fracmax))
	plt.show()
	'''	
	print f[f['id'] == 302]
	highztable = f[ (f['z_p'] > zllim) * (f['z_p'] < zulim) * (f['wht_err'] < 2)] #sort out redshifts we don't want, and then sort out images with more than 1 band error
	
	return highztable
	
def obj_lookup(dir, obj):
	import pysao
	cat = table.Table.read(dir + '/A2744_cat.dat', format = 'ascii')['id', 'X_IMAGE', 'Y_IMAGE']
	obj = cat[cat['id'] == obj]
	X_IMAGE = str(np.rint(obj['X_IMAGE'][0]))
	Y_IMAGE = str(np.rint(obj['Y_IMAGE'][0]))
	#print obj
	ds9 = pysao.ds9()
	ds9.set('file 160.fits')
	ds9.set('regions file zheng2014.reg')
	ds9.set('scale log 99.5')
	ds9.set('scale limits -.1 10')
	ds9.set('cmap value 10.0889 0.358543')
	#ds9.xpa_help('pan')
	panstring = 'pan to ' + X_IMAGE + ' ' + Y_IMAGE + ' image'
	ds9.set(panstring)
	a = raw_input('Press ENTER to exit...')
	#print ds9.get('cmap value')
	ds9.set('exit')

dir = 'eazy-cat-EA'
'''
highztable = find_highz(dir, 'OUTPUT/eazy-cat.zout', 'A2744_cat.dat', 3., 10.)
#print highztable
print 'Number of candidates:', len(highztable['id'])
'''

def filter_sortby_z(dir, zllim = 7., zulim = 10.):
	f = table.Table.read(dir + '/OUTPUT/eazy-cat.zout', format = 'ascii')
	f = f[ (f['z_peak'] > zllim) * (f['z_peak'] < zulim) ]
	return f['id'].tolist()

#galplot(dir, 462)
#galplot(dir, 187)
#galplot(dir, 788)
#obj_lookup(dir, 788)
#galplot(dir, 302)

#galplot(dir, 440)
#obj_lookup(dir, 440) #440 is genuine, looks like prior exerts influence

galplot(dir, 457)
obj_lookup(dir, 457) #440 is genuine, looks like prior exerts influence

'''
obj_list = filter_sortby_z(dir)
for obj in pbar(obj_list):
	galplot(dir, str(obj), save = True, plot_dir = 'z_7_10_figs')
'''
'''
for obj in pbar(highztable['id']):
	galplot(dir, obj, save = True)
	#if obj > 100: break
'''