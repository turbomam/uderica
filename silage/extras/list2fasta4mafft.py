import fileinput
import re

for line in fileinput.input():
	fields = line.split('\t')
	fullpath = fields[0]
	p = re.compile('^.*\/')
	basename = p.sub('', fields[0])
	print ">{}{}".format(basename,basename)
