import numpy as np
import glob as glob
import re

beginstring = ('# Template definition file'
'\n#'
'\n# No blank lines allowed (for now).'
'\n#'
'\n# Column definitions:'
'\n#   1. Template number'
'\n#   2. Template file name'
'\n#   3. Lambda_conv (multiplicative factor to correct wavelength units)'
'\n#   4. Age of template model in Gyr (0 means template is always used)'
'\n#   5. Template error amplitude (for INDIVIDUAL template fits)'
'\n#   6. Comma/space separated list of template numbers to be combined'
'\n#      with current template if combined fits are enabled.'
'\n#'
'\n# Sample entry:'
'\n# 1 [path_to_file]/template1.sed 1.0 14.7 0.2 2,3,5')

relpath = 'BC03Templates/'
template_list = glob.glob('templates/BC03Templates/BC*')
num_templates = len(template_list)
#print num_templates

def read_template(fname, template_number, num_templates):
	f = open(fname, 'r')
	lambda_conv = 1.0
	epoch = 0
	for i, line in enumerate(f):
		if 'Epoch (yr)' in line:
			match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
			num = [float(x) for x in re.findall(match_number, line)]
			epoch = float(num[0]*10**num[1]) / 1e9
			#print f, epoch
	usewith = ( str(range(template_number+1, num_templates + 1))[1:-1] ).replace(', ', ',')
	#print template_number, epoch
	print ' '.join((str(template_number), fname, str(lambda_conv), str(epoch), str(1.0), usewith))

for i, item in enumerate(template_list):
	#print item
	read_template(item, i+1, num_templates)