#!/usr/bin/env python
from sv2_backend import check_PAR,errFH,format_chrom,tokenize 
from pybedtools import BedTool
import sys
__help_url__='https://github.com/dantaki/SV2/wiki/input#sv-input' 
class Structural_Variant():
	def __init__(self,sv=None,variant_id=None,gen=None):
		self.raw=[(str(x[0]),int(x[1]),int(x[2]),str(x[3])) for x in BedTool(list(set(sv))).sort()]
		self.variant_id=variant_id
		self.par={} # [locus]=True/False; true if intersects PAR 
		for locus in self.raw:
			chrom = format_chrom(locus[0])
			if chrom != 'chrX' and chrom != 'chrY': continue
			key = (chrom,int(locus[1]),int(locus[2]),str(locus[3]))
			self.par[key]=check_PAR('{} {} {}'.format(chrom,locus[1],locus[2]),gen) 
def checkCall(cl,l,fh):
	keep=True
	if 'DEL' not in cl and 'DUP' not in cl: 
		sys.stderr.write('WARNING: Skipping bad entry {} in {}. {} not an accepted SV type. Accepted types: DEL or DUP: {}\n'.format(l.rstrip(),fh,cl,__help_url__))
		keep=False
	return keep
def Bed(fh,sv):
	errFH(fh)
	with open(fh,'r') as f:
		for l in f:
			if l.startswith('#'): continue
			r = tokenize(l)
			if r==0 or len(r)<4: 
				sys.stderr.write('WARNING: Skipping bad entry {} in BED file {}. BED format requirements: {}\n'.format(l.rstrip(),fh,__help_url__))
			elif r!=0 and checkCall(str(r[3]),l,fh)!=False: sv.append((str(r[0]),int(r[1]),int(r[2]),str(r[3])))
def Vcf(fh,sv,variant_id):
	errFH(fh)
	vcf_fh=None
	if fh.endswith('.gz'):
		import gzip
		vcf_fh = gzip.open(fh,'r')
	else: vcf_fh = open(fh,'r')
	with vcf_fh as f:
		for l in f:
			if l.startswith('#'): continue
			s,e,r=-9,-9,tokenize(l)
			if r!=0:
				c,s=str(r[0]),int(r[1])
				s-=1
				if len(r)<8: sys.stderr.write('WARNING: Skipping bad entry {} in VCF file {}. INFO column is required: {}\n'.format(l.rstrip(),fh,__help_url__))
				else:
					for i in r[7].split(';'):
						if 'SVTYPE=' in i: cl=i.replace('SVTYPE=','')
						if i.startswith('END=')and 'CI' not in i: e=int(i.replace('END=',''))
					if s==-9 or e==-9: sys.stderr.write('WARNING: Skipping bad entry {} in VCF file {}. SV entry missing positions: {}\n'.format(l.rstrip(),fh,__help_url__))
					else:
						if checkCall(cl,l,fh)==True:
							sv.append((c,s,e,cl))
							variant_id[(str(c),int(s),int(e),str(cl))]=str(r[2])
	vcf_fh.close()
def skip_sv(c,s,e,cl,limit):
	skip=False
	if e<=s or s>=e: 
		skip=True
		sys.stderr.write('WARNING: end position is less than or equal to start position. Error found here: {}:{}-{} {}. Skipping ...\n'.format(c,s,e,cl))
	if s == limit or e >= limit:
		skip=True
		sys.stderr.write('WARNING: positions surpass reference limit. Error found here: {}:{}-{} {}. Skipping ...\n'.format(c,s,e,cl))
	if s<0 or e<0:
		skip=True
		sys.stderr.write('WARNING: positions are less than 0. Error found here: {}:{}-{} {}. Skipping ...\n'.format(c,s,e,cl))
	return skip
def sv_init(beds=None,vcfs=None,gen=None):
	sv,variant_id=[],{}
	if beds!=None:
		for fh in beds: Bed(fh,sv)
	if vcfs!=None:
		for fh in vcfs: Vcf(fh,sv,variant_id)
	if len(sv)<1: 
		print 'FATAL ERROR: SVs were not formatted correctly. See https://github.com/dantaki/SV2/wiki/input#sv-input for details'
		sys.stderr.write('FATAL ERROR: SVs were not formatted correctly. See https://github.com/dantaki/SV2/wiki/input#sv-input for details\n')
		sys.exit(1)
	SV = Structural_Variant(sv,variant_id,gen)
	return SV