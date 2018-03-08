#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
from sv2_backend import errFH,format_chrom,match_chrom_key,match_chrom_prefix
import numpy as np
import os,pysam,sys
from pybedtools import BedTool
cimport numpy as np
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
def extract_snv_features(Bam,svs,master_sv,Pre,gen): return c_extract_snv_features(Bam,svs,master_sv,Pre,gen)
def preprocess_snv(Bam,masked): return c_preprocess_snv(Bam,masked)
def snv_init(l=None):
	snv=[]
	for fh in l:
		errFH(fh)
		try: pysam.TabixFile(fh)
		except IOError:
			sys.stderr.write('WARNING: tabix index file for {} not found. SNV files must be bgzipped and indexed with tabix. Skipping ...\n'.format(fh))
			continue
		Vcf = Snv(fh)
		if Vcf.ad_format!=None: snv.append(Vcf)
	return snv
cdef c_extract_snv_features(Bam,svs,master_sv,Pre,gen):
	cdef unsigned int s=0
	cdef unsigned int e=0
	cdef unsigned int depth=0
	cdef double ratio=0.
	cdef double norm_cov=0.
	snv_feats,het_feats={},{}
	Itr = pysam.TabixFile(Bam.Snv.fh)
	for call in svs:
		# foreach SV
		if master_sv.get(call)==None: continue
		sv_span,flank_span,windows = master_sv[call]
		locus_depth, locus_had=[],[]
		chrom = format_chrom(call[0])
		chrom_key = match_chrom_key(call,Bam,gen)
		if Pre.snv_cov[(Bam.id,chrom_key)]==None: continue
		else: norm_cov=Pre.snv_cov[(Bam.id,chrom_key)]
		#####################
		for locus in sv_span:
			# foreach masked region within SV
			c,s,e = str(locus[0]),int(locus[1]),int(locus[2])
			c=match_chrom_prefix(c,Bam.Snv.chr_flag)
			for Variant in Itr.fetch(region='{}:{}-{}'.format(c,s+1,e)):
				data= tokenize_vcf(Bam,Variant)
				if data==0:continue
				depth,ratio= data
				locus_depth.append(depth)
				if ratio!=-9:locus_had.append(ratio)
		#####################
		if len(locus_depth)==0: snv_feats[call]=(float('nan'),0)
		else: snv_feats[call]=(np.nanmedian(locus_depth)/norm_cov, len(locus_depth)) # normalized coverage, # snvs
		if len(locus_had)==0: het_feats[call]=(float('nan'),0)
		else: het_feats[call]=(np.nanmedian(locus_had),len(locus_had)) # median heterozygous allele depth, # heterozygous snvs
	Itr.close()
	return snv_feats,het_feats
cdef c_preprocess_snv(Bam,masked):
	cdef int genome_lens=0
	genome_depths=[]
	depths,snv_depth,genome={},{},[] # depth[ref]=[], snv_depth[ref]= median depth, # snvs
	Itr = pysam.TabixFile(Bam.Snv.fh)
	for entry in masked:
		c = str(entry[0])
		if Bam.chr_flag==False: c=c.replace('chr','')
		chrom = format_chrom(c)
		if depths.get(chrom)==None: depths[chrom]=[]
		region = '{}:{}-{}'.format(c,int(entry[1])+1,int(entry[2]))
		try:
			for Variant in Itr.fetch(region=region):
				data = tokenize_vcf(Bam,Variant,False) # False=only process depth; output is (depth, allele depth)
				if data == 0: continue
				depths[chrom].append(data[0])
		except ValueError: sys.stderr.write('WARNING: skipping {}. Region not found in {}\n'.format(region,Bam.Snv.fh))
	Itr.close()
	for chrom in depths:
		if len(depths[chrom])>0:
			snv_depth[chrom]= np.nanmedian(depths[chrom]),len(depths[chrom])
			genome_depths+=depths[chrom]
			genome_lens+=len(depths[chrom])
		else: snv_depth[chrom]=0,0
	if len(genome_depths)>0: snv_depth['GENOME']=np.nanmedian(genome_depths),genome_lens
	else: snv_depth['GENOME']=0,0
	return snv_depth	
cdef format_genotype(x): return x.replace('|','/')
cdef tokenize_vcf(Bam,Variant,HAD=True):
	cdef double ratio=-9
	cdef unsigned int depth
	cdef short depth_index=-9
	cdef short genotype_index=-9
	cdef short allele_depth_index=-9
	record=Variant.split('\t')
	try: record[Bam.snv_index]
	except IndexError:
		sys.stderr.write('WARNING: {} SNV VCF entry is malformed. Skipping ...\n'.format('\t'.join(record)))
		return 0
	if len(record) < Bam.snv_index: return 0
	c,s  = format_chrom(record[0]),int(record[1])
	_format = record[8].split(':')
	if HAD == True: # process allele depth and depth
		for ind in xrange(len(_format)):
			if _format[ind]=='DP': depth_index=ind
			if _format[ind]=='GT': genotype_index=ind
			if _format[ind]== Bam.Snv.ad_format: allele_depth_index=ind
		if depth_index==-9 or genotype_index==-9 or allele_depth_index == -9: return 0
		sample_entry = record[Bam.snv_index].split(':')
		if '.' in sample_entry[genotype_index] or '.' in sample_entry[depth_index]: return 0
		genotype=[int(x) for x in format_genotype(sample_entry[genotype_index]).split('/')]
		if len(genotype)!=2: return 0
		depth = int(sample_entry[depth_index])
		if genotype[0]!=genotype[1]:
			allele_depth = [float(x) for x in sample_entry[allele_depth_index].split(',')]
			a,b = allele_depth[genotype[0]],allele_depth[genotype[1]]
			if a== 0.0 or b == 0.0: return 0
			ratio = b/a
			if b > a: ratio=a/b
	else: # only process depth
		for ind in xrange(len(_format)):
			if _format[ind]=='DP': depth_index=ind
			if _format[ind]=='GT': genotype_index=ind
		if depth_index==-9: return 0
		sample_entry = record[Bam.snv_index].split(':')
		if '.' in sample_entry[genotype_index] or '.' in sample_entry[depth_index]:return 0
		depth=int(sample_entry[depth_index])
	return (depth,ratio)