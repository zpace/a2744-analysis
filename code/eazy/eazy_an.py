import numpy as np
import astropy.io.ascii as apa
import matplotlib.pyplot as plt

dirnames = ['y1_K','y1_K_newACS','y1_newACS','y1_old']
fnames = [dirname + '/OUTPUT/1.pz' for dirname in dirnames]

#build DF of fit data
zdata = apa.read(fnames[0], include_names = ['z', 'prior'])
for i in range(0, len(dirnames)):
	tempdata = apa.read(fnames[i])
	zdata[dirnames[i]+'_chi2'] = tempdata['chi2']
	zdata[dirnames[i]+'_pz'] = tempdata['pz']


#print zdata
plt.plot(zdata['z'], zdata['prior'], label = 'prior')
for dirname in dirnames:
	plt.plot(zdata['z'], zdata[dirname + '_chi2']/np.max(zdata[dirname + '_chi2']), label = dirname + '_chi2')
	print dirname + ': min err', np.min(zdata[dirname + '_chi2']), 'at z =', zdata['z'][np.argmin(zdata[dirname + '_chi2'])]

plt.xlabel('redshift')
plt.xlabel('prior & relative chi2 error')
plt.title('propagation of prior')
plt.legend(loc = 'best')
plt.show()

print '\n=====\n'

for dirname in dirnames:
	print dirname
	logchi2 = np.exp(-zdata[dirname + '_chi2']/2.)
	cumul = np.trapz(logchi2, zdata['z'])
	logchi2 *= 1./cumul
	plt.plot(zdata['z'], logchi2, label = dirname)
	#probability that z is more than 4
	p = np.trapz(logchi2[zdata['z'] > 4.], zdata['z'][zdata['z'] > 4.])
	print 'Probability of z > 4:', p
	
	#probability that z is between 7.5 and 8.5
	y = logchi2[(zdata['z'] > 7.5) * (zdata['z'] < 8.5)]
	x = zdata['z'][(zdata['z'] > 7.5) * (zdata['z'] < 8.5)]
	p = np.trapz(y, x)
	print 'Probability of 7.5 < z < 8.5:', p
	
	#probability that z is between 8 and 8.2
	y = logchi2[(zdata['z'] > 8) * (zdata['z'] < 8.2)]
	x = zdata['z'][(zdata['z'] > 8) * (zdata['z'] < 8.2)]
	p = np.trapz(y, x)
	print 'Probability of 8.0 < z < 8.2:', p
	
	#probability that z is between 8.05 and 8.13
	y = logchi2[(zdata['z'] > 8.05) * (zdata['z'] < 8.13)]
	x = zdata['z'][(zdata['z'] > 8.05) * (zdata['z'] < 8.13)]
	p = np.trapz(y, x)
	print 'Probability of 8.05 < z < 8.13:', p
	
plt.legend(loc = 'best')
plt.xlabel('z')
plt.ylabel('Probability')
plt.title('Redshift of Abell 2744-Y1')
#plt.semilogy()
plt.show()