__version__='1.2'
import sys,os,argparse
from .core import writeConfig,checkConfig,check_in,Bed,check_cnv,errFH,reportTime,preprocess,extract_feats,genotype,annotate
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
from glob import glob
from time import time
def main():
	init_time = int(time())
	splash='\n                       ____\n  _____________   ___ |___ \\\n /   _____/\   \ /   // ___/\n \_____  \  \   Y   //_____)\n /        \  \     /\n/_________/   \___/\nSupport Vector Structural Variation Genotyper\nVersion 1.2        Author: Danny Antaki <dantaki at ucsd dot edu>\n'
	parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
	genoArgs = parser.add_argument_group('genotype arguments')
	configArgs = parser.add_argument_group('configure arguments')
	genoArgs.add_argument('-i','-in', help='Tab delimited input [ ID, BAM-PATH, VCF-PATH, M/F ]')
	genoArgs.add_argument('-r','-cnv', help='SV to genotype. Either in BED or VCF format',type=str)
	genoArgs.add_argument('-c','-cpu', help='Parallelize sample-wise. 1 per cpu',required=False,default=1,type=int)
	genoArgs.add_argument('-g','-genome',  help='Reference genome build [ hg19, hg38 ]',required=False,default='hg19',type=str)
	genoArgs.add_argument('-pcrfree',  help='GC content normalization for PCR free libraries',required=False,default=False,action="store_true")
	genoArgs.add_argument('-s','-seed', help='Preprocessing: integer seed for genome shuffling',required=False,default=42,type=int)
	genoArgs.add_argument('-o','-out', help='output',required=False,default="sv2_genotypes.vcf",type=str)
	genoArgs.add_argument('-pre', help='Preprocessing output directory',required=False,default=None)
	genoArgs.add_argument('-feats', help='Feature extraction output directory',required=False,default=None)
	configArgs.add_argument('-hg19',default=None,help='hg19 FASTA',required=False)
	configArgs.add_argument('-hg38',default=None,help='hg38 FASTA',required=False)	
	args = parser.parse_args()
	infh = args.i
	bed = args.r
	cores = args.c
	gen = args.g
	pcrfree = args.pcrfree
	ofh = args.o
	seed = args.s
	predir= args.pre
	featsdir = args.feats
	conf_hg19 = args.hg19
	conf_hg38 = args.hg38
	preprocess_files={}
	feats_files={}
	gens = ['hg19','hg38']
	if gen not in gens: 
		sys.stderr.write('ERROR -g must be either hg19 or hg38. NOT {}\n'.format(gen))
		sys.exit(1)
	if conf_hg19 != None or conf_hg38 != None: 
		writeConfig(conf_hg19,conf_hg38)
		sys.exit(0)
	checkConfig()
	bam_dict,vcf_dict,gender_dict=check_in(infh)
	raw,cnv=check_cnv(Bed(bed),gen)
	ofh = ofh.replace('.txt','.vcf').replace('.out','.vcf')
	if not ofh.endswith('.vcf'): ofh=ofh+'.vcf'
	"""
	PREPROCESSING
	"""
	if predir == None:
		if cores > 1 :
			pool = Pool(processes=cores)
	       		for bam_id in bam_dict:
				preofh = bam_id+'_sv2_preprocessing.txt'
	      			preprocess_files[bam_id]=os.getcwd()+'/sv2_preprocessing/'+preofh
				pool.apply_async(preprocess, args=(bam_id,bam_dict[bam_id],vcf_dict[bam_id],preofh,gen,seed) )
			pool.close()
			pool.join()
		else:
			for bam_id in bam_dict:
				preofh = bam_id+'_sv2_preprocessing.txt'
				preprocess_files[bam_id]=os.getcwd()+'/sv2_preprocessing/'+preofh
				preprocess(bam_id,bam_dict[bam_id],vcf_dict[bam_id],preofh,gen,seed)	
	else: 
		if not predir.endswith('/'): predir = predir+'/'
		if not os.path.isdir(predir): errFH(predir)
		for fh in glob(predir+'*sv2_preprocessing.txt'):
			f = open(fh)
			if sum(1 for l in open(fh)) <= 1: continue
			else:
				preids=[]
				head=f.next().rstrip('\n')
				for l in f: preids.append(l.rstrip('\n').split('\t').pop(0))
			f.close()
			for iid in set(preids): 
				if gender_dict.get(iid) != None: preprocess_files[iid]=fh
	reportTime(init_time,'PREPROCESSING COMPLETE')
	""""
	FEATURE EXTRACTION
	"""
	if bed == None: errFH(bed) 
	if featsdir == None:
		if cores > 1 :
			pool = Pool(processes=cores)
      			for bam_id in bam_dict:
				if preprocess_files.get(bam_id) == None:
					print 'WARNING: preprocessing iid {} does not match inlist iid'.format(bam_id)
					continue
				prefh = preprocess_files[bam_id]
				gtofh = bam_id+'_sv2_features.txt'
				feats_files[bam_id]=os.getcwd()+'/sv2_features/'+gtofh
				pool.apply_async(extract_feats, args=(bam_id,bam_dict[bam_id],vcf_dict[bam_id],cnv,prefh,gender_dict[bam_id],gtofh,gen,pcrfree) )
			pool.close()
			pool.join()
		else:
			for bam_id in bam_dict:
				if preprocess_files.get(bam_id) == None:
					print 'WARNING: preprocessing iid {} does not match inlist iid'.format(bam_id)
					continue
				prefh = preprocess_files[bam_id]
				gtofh = bam_id+'_sv2_features.txt'
				feats_files[bam_id]=os.getcwd()+'/sv2_features/'+gtofh
				extract_feats(bam_id,bam_dict[bam_id],vcf_dict[bam_id],cnv,prefh,gender_dict[bam_id],gtofh,gen,pcrfree)
	else: 
		if not featsdir.endswith('/'): featsdir=featsdir+'/'
		if not os.path.isdir(featsdir): errFH(featsdir)
		for fh in glob(featsdir+'*sv2_features.txt'):
			f = open(fh)
			if sum(1 for l in open(fh)) <= 1: continue
			else:
				featsid=[]
				head=f.next().rstrip('\n')
				for l in f: featsid.append(l.rstrip('\n').split('\t').pop(5))
			f.close()
			for iid in set(featsid): 
				if gender_dict.get(iid) != None: feats_files[iid]=fh
	reportTime(init_time,'FEATURE EXTRACTION COMPLETE')	
	"""
	GENOTYPING
	"""
	feats=[]
	for iid in feats_files:
		with open(feats_files[iid]) as f:
			for l in f: feats.append(tuple(l.rstrip('\n').split('\t')))
	genos,REF,NON,GQ,HEMI,STD_FILT,DNM_FILT = genotype(raw,feats,gender_dict,gen,ofh.replace('.vcf','.txt'))
	annotate(raw,genos,gen,REF,NON,GQ,ofh,gender_dict,HEMI,STD_FILT,DNM_FILT)
	reportTime(init_time,'GENOTYPING COMPLETE')	
