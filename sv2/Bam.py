#!/usr/bin/env python
from sv2_backend import accepted_chrom,errFH,format_chrom,make_dir
from collections import OrderedDict
import pysam,sys
class Bam():
	def __init__(self,fh=None,Ped=None,gen=None):
		self.id = None
		self.chr_flag=False # True == 'chrN'; False == 'N'
		self.sex=None
		self.char='b'
		self.refs=OrderedDict()
		self.fh=fh
		self.Snv=None
		self.snv_index=None
		sm=[]
		bam = pysam.AlignmentFile(fh)
		if bam.is_cram==True: self.char='c'
		header = bam.header
		if header.get('RG')==None: sys.stderr.write('WARNING: {} lacks Read Group (@RG) entry in the header. Skipping ...\n'.format(fh))
		else:
			for entry in header['RG']:
				if entry.get('SM')!=None: sm.append(entry['SM'])
			sm = list(set(sm))
			if len(sm)>1: sys.stderr.write('WARNING: {} contains two sample entries ({}) in the header. Skipping {} ...\n'.format(fh,sm,fh))
			else: self.id=sm[0]
		if self.id!=None:
			if Ped.sex.get(self.id)==None: sys.stderr.write('WARNING: {} is not found in the PED file. Skipping {} ...\n'.format(self.id,fh))
			else:
				self.sex=Ped.sex[self.id]
				if str(bam.references[0]).startswith('chr'): self.chr_flag=True
				chroms = accepted_chrom(self.chr_flag,gen)
				for i in range(len(bam.references)):
					chrom,leng = str(bam.references[i]),bam.lengths[i]
					if chrom in chroms: self.refs[chrom]=int(leng)
		bam.close()
	def tmp_chrom_file(self,tmp_dir=None,genome=True,chrom=None):
		make_dir(tmp_dir)
		tmp_genome = tmp_dir+'{}.genome'.format(self.id)
		tmpfh = open(tmp_genome,'w')
		if genome==True:
			if chrom==None:
				for ref in self.refs: tmpfh.write('{}\t{}\n'.format(format_chrom(ref),self.refs[ref]))
			else: tmpfh.write('{}\t{}\n'.format(format_chrom(chrom),self.refs[chrom]))
		if genome==False:
			for ref in self.refs:tmpfh.write('{}\t0\t{}\n'.format(format_chrom(ref),self.refs[ref]))
		tmpfh.close()
		return tmp_genome
def bam_init(args=None,Ped=None,Snv=None,gen=None):
	bams=[]
	for f in args:
		errFH(f)
		try: pysam.AlignmentFile(f).check_index()
		except ValueError:
			sys.stderr.write('WARNING: {} is not indexed with samtools index. Skipping ...\n'.format(f))
			continue
		except AttributeError:
			sys.stderr.write('WARNING: {} appears to be in SAM format. Convert to BAM with `samtools view -bh {}` and index with samtools index\n'.format(f,f))
			continue
		bam = Bam(f,Ped,gen)
		if len(Snv)<1:
			print 'FATAL ERROR: SNV file(s) were not formatted correctly. See https://github.com/dantaki/SV2/wiki/input#snv-vcf for details'
			sys.stderr.write('FATAL ERROR: SNV file(s) were not formatted correctly. See https://github.com/dantaki/SV2/wiki/input#snv-vcf for details\n')
			sys.exit(1)
		for snv in Snv:
			if snv.id.get(bam.id)!=None:
				bam.Snv,bam.snv_index=snv,snv.id[bam.id]
		if bam.id!=None and bam.Snv==None: sys.stderr.write('WARNING: BAM file {} sample name (@RG  SM:<sample_id>):{} not found in SNV VCFs. Skipping {} ...\n'.format(f,bam.id,f))
		if bam.id==None: sys.stderr.write('WARNING: Skipping BAM file {}. No sample name (@RG SM:<sample_id>). See https://github.com/dantaki/SV2/wiki/input#bam for details')
		if bam.id!=None and bam.Snv!=None: bams.append(bam)
	if len(bams)<1:
		print 'FATAL ERROR: BAM file(s) were not formatted correctly. See https://github.com/dantaki/SV2/wiki/input#bam for details'
		sys.stderr.write('FATAL ERROR: BAM file(s) were not formatted correctly. See https://github.com/dantaki/SV2/wiki/input#bam for details\n')
		sys.exit(1)
	return bams