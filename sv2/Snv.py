#!/usr/bin/env python2
from Backend import errFH
import pysam,sys
class Snv():
	def __init__(self,fh=None):
		isHC,isFB=False,False
		self.fh=fh
		self.id={} # id[IID]=ROW INDEX
		self.chr_flag=False # True == 'chrN'; False == 'N'
		self.ad_format=None
		tbx = pysam.TabixFile(fh)
		header = list(tbx.header)
		if [x for x in header if 'FORMAT=<ID=AD' in x]: isHC=True
		if [x for x in header if 'FORMAT=<ID=DPR' in x]: isFB=True
		if isHC==isFB: sys.stderr.write('WARNING:{} SNV VCF does not contain Allele Depth. See details for SNV VCF format: github.com/dantaki/SV2/wiki/input#snv-vcf\nSkipping {} ...'.format(fh,fh))
		else:
			if isHC==True:self.ad_format='AD'
			if isFB==True:self.ad_format='DPR'
			head = [x for x in header if x.startswith('#CHROM')]
			if len(head)<1:
				sys.stderr.write('WARNING:{} SNV VCF does not contain a proper header. See details for SNV VCF format: github.com/dantaki/SV2/wiki/input#snv-vcf\nSkipping {} ...'.format(fh,fh))
				self.ad_format=None
			else:
				head = head[0].split('\t')
				for index in range(9,len(head)): self.id[head[index]]=index
				if str(tbx.contigs[0]).startswith('chr'): self.chr_flag=True
		tbx.close()		
def snv_init(l=None):
	snv=[]
	for fh in l:
		errFH(fh)
		Vcf = Snv(fh)
		if Vcf.ad_format!=None: snv.append(Vcf)
	return snv
