#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
from sv2_backend import format_chrom,get_path,make_dir,mask_bed,match_chrom_prefix
from Snv import preprocess_snv
import numpy as np
import os,pysam,sys
from pybedtools import BedTool
cimport numpy as np
np.import_array()
class Preprocess():
	def __init__(self,pre=None):
		chr_cov,insert_size,mad,read_length,snv_cov={},{},{},{},{}
		with open(pre) as f:
			for l in f:
				if l.startswith('#'):continue
				r = l.rstrip('\n').split('\t')
				chr_cov[(str(r[0]),str(r[1]))]=float(r[2])
				snv_cov[(str(r[0]),str(r[1]))]=float(r[7])
				if r[1] == "GENOME":
					mad[r[0]]=float(r[5])
					insert_size[r[0]]=float(r[4])
					read_length[r[0]]=float(r[3])
		self.chr_cov=chr_cov
		self.insert_size=insert_size
		self.insert_mad=mad
		self.read_len=read_length
		self.snv_cov=snv_cov
def preprocess(Bam,ofh,seed,gen,tmp_dir):
	outfh = open(ofh,'w')
	outfh.write('\t'.join(("#id","chrom","coverage_median","read_length_median","insert_size_median","insert_size_mad","bam_bp_parsed","snv_depth_median","snv_parsed"))+'\n')
	outfh.close()
	sv2_preprocess(Bam,ofh,seed,gen,tmp_dir)
cdef MAD(a, double c=0.6745): return np.around(np.median(np.fabs(a-np.median(a))/c),decimals=3)
cdef normalize_chr_cov(double read_count, double chr_size, read_length): return np.around( (read_count/chr_size)*np.median(read_length),decimals=2)
cdef random_int(): return ''.join(map(str,np.random.randint(9,size=6)))
cdef sv2_preprocess(Bam,ofh,seed,gen,tmp_dir):
	cdef unsigned int chr_size
	cdef unsigned int read_count
	cdef unsigned int genome_size
	genome_cov,genome_read_length,genome_insert_size,genome_mad,genome_size=[],[],[],[],0
	tmp_bed = Bam.tmp_chrom_file(tmp_dir,False) # tmp file chrom 0 chrom_len
	snv_depth=preprocess_snv(Bam,mask_bed(tmp_bed,gen))
	os.remove(tmp_bed)
	out = open(ofh,'a')
	Itr = pysam.AlignmentFile(Bam.fh,'r{}'.format(Bam.char))  
	from sv2Config import Config
	for chrom in Bam.refs:
		tmp_genome = Bam.tmp_chrom_file(tmp_dir,True,chrom)
		chrom=format_chrom(chrom)
		chr_size,read_count,read_stats,read_length,insert_size=0,0,{},[],[]
		rand_bed = BedTool()
		rand_bed = rand_bed.random(l=100000,n=100,seed=seed,g=tmp_genome)
		rand_bed = rand_bed.subtract(BedTool('{}{}_excluded.bed.gz'.format(Config().resource_path(),gen))).sort().merge()
		for entry in rand_bed:
			entry=tuple(entry)
			c,s,e = str(entry[0]),int(entry[1]),int(entry[2])
			if Bam.chr_flag==False: c=c.replace('chr','')
			region = '{}:{}-{}'.format(c,s+1,e)
			chr_size+=e-s
			for Aln in Itr.fetch(region=region):
				if (Aln.is_proper_pair == False or Aln.is_qcfail== True or Aln.is_duplicate == True or Aln.mapq < 40 or read_stats.get(str(Aln.qname)+str(Aln.is_read1)) != None): continue
				read_stats[str(Aln.qname)+str(Aln.is_read1)] = ((Aln.qlen,abs(Aln.isize)))
		read_count = len(read_stats)
		if read_count == 0: out.write('\t'.join(map(str,(Bam.id,chrom,'0','0','0','0',chr_size,snv_depth[chrom][0],snv_depth[chrom][1])))+'\n')
		else :
			for name in read_stats:
				read_length.append(read_stats[name][0])
				insert_size.append(read_stats[name][1])
			norm = normalize_chr_cov(read_count,chr_size,np.array(read_length))
			out.write('\t'.join(map(str,(Bam.id,chrom,norm,np.median(read_length),np.median(insert_size),MAD(np.array(insert_size)),chr_size,snv_depth[chrom][0],snv_depth[chrom][1])))+'\n')
			genome_cov.append(norm)
			genome_read_length.append(np.median(read_length))
			genome_insert_size.append(np.median(insert_size))
			genome_mad.append(MAD(np.array(insert_size)))
			genome_size += chr_size
		os.remove(tmp_genome)
	out.write('\t'.join(map(str,(Bam.id,'GENOME',np.median(genome_cov),np.median(genome_read_length),np.median(genome_insert_size),np.median(genome_mad),genome_size,snv_depth['GENOME'][0],snv_depth['GENOME'][1])))+'\n')
	out.close()
	Itr.close()