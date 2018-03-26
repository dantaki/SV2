#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
from sv2_backend import check_PAR,format_chrom
from Svm import biallelic_del_svm,biallelic_dup_svm,male_sex_chrom_svm
from Variant import Variant
import numpy as np
cimport numpy as np
import pandas as pd
import os,sys
np.import_array()
def init_dataframe(obj,header=['chr','start','end','type','size','id','old_covr','covr_GC','dpe','sr','SNP_coverage','HET_ratio','SNPs','HET_SNPs']):
	df = pd.DataFrame(obj)
	df.columns=header
	return df
def genotype(features,out,classifier_name,Ped,Conf,gen):
	SVs={} # SVs[SV.locus]=Variant object
	Conf.get_clf(classifier_name)
	biallelic_dels,biallelic_dups,male_sex_chrom_dels,male_sex_chrom_dups=partition_svs(features,Ped,gen)
	pd.options.mode.chained_assignment = None
	genos=[]
	if not out.endswith('.txt'): out = out+'.txt'
	ofh = open(out,'w')
	ofh.write('\t'.join(('#CHROM','START','END','TYPE','LENGTH','ID','COVERAGE','DISCORDANT_PAIRED-END','SPLIT_READS','SNV_COVERAGE','SNVs','HETEROZYGOUS_ALLELE_RATIO','HETEROZYGOUS_SNVs','COPY_NUMBER_GENOTYPE','REF_GENOTYPE_LIKELIHOOD','HET_GENOTYPE_LIKELIHOOD','HOM_GENOTYPE_LIKELIHOOD','ALT_GENOTYPE_LIKELIHOOD','REF_QUAL','ALT_QUAL','CLASSIFIER'))+'\n')
	if len(biallelic_dels)>0:
		c_biallele_dels(biallelic_dels,Ped,Conf,SVs,ofh)
	if len(biallelic_dups)>0:
		c_biallele_dups(biallelic_dups,Ped,Conf,SVs,ofh)
	if len(male_sex_chrom_dels)>0:
		c_male_sex_chrom_dels(male_sex_chrom_dels,Ped,Conf,SVs,ofh)
	if len(male_sex_chrom_dups)>0:
		c_male_sex_chrom_dups(male_sex_chrom_dups,Ped,Conf,SVs,ofh)
	ofh.close()
	for locus in SVs:
		Variant = SVs[locus]
		Variant.median_likelihood()
		Variant.add_standard_filters()
		Variant.add_denovo_filters()
	return SVs
def partition_svs(features,Ped,gen):
	biallelic_dels,biallelic_dups,male_sex_chrom_dels,male_sex_chrom_dups = [],[],[],[]
	for x in features:
		chrom,svtype,sample_id=format_chrom(x[0]),str(x[3]),str(x[5])
		if 'DEL' in svtype:
			if (chrom=='chrX' or chrom=='chrY') and Ped.males.get(sample_id)!=None and check_PAR('{} {} {}'.format(chrom,x[1],x[2]),gen)==False: male_sex_chrom_dels.append(x)
			else: biallelic_dels.append(x)
		elif 'DUP' in svtype:
			if (chrom=='chrX' or chrom=='chrY') and Ped.males.get(sample_id)!=None and check_PAR('{} {} {}'.format(chrom,x[1],x[2]),gen)==False: male_sex_chrom_dups.append(x)
			else: biallelic_dups.append(x)
	return biallelic_dels,biallelic_dups,male_sex_chrom_dels,male_sex_chrom_dups
cdef biallele(df,Ped,SVs,ofh):
	for x in df.values:
		ofh.write('\t'.join(map(str,x))+'\n')
		SV = Variant(x)
		if SVs.get(SV.locus)!=None: SV=SVs[SV.locus]
		clf=str(x[20])
		SV.clf.append(clf)
		cn = format(float(x[6])*2,'.2f')
		svgt = int(x[13])
		skip=False
		if svgt==2:
			lik=np.rint(float(x[18]))
			skip=set_biallelic_genotype(SV,x[0],x,Ped,'0/0',cn,lik)
			if skip==False:
				if SV.ref.get(clf)==None: SV.ref[clf]=[lik]
				else: SV.ref[clf].append(lik)
		else:
			lik=np.rint(float(x[19]))
			if svgt==1 or svgt==3: skip=set_biallelic_genotype(SV,x[0],x,Ped,'0/1',cn,lik)
			elif svgt==0 or svgt==4: skip=set_biallelic_genotype(SV,x[0],x,Ped,'1/1',cn,lik)
			if skip==False:
				
				if SV.alt.get(clf)==None: SV.alt[clf]=[lik]
				else: SV.alt[clf].append(lik)
		SVs[SV.locus]=SV
cdef c_biallele_dels(biallelic_dels,Ped,Conf,SVs,ofh):		
	biallele(biallelic_del_svm(init_dataframe(biallelic_dels),Conf.gtclf),Ped,SVs,ofh)	
cdef c_biallele_dups(biallelic_dups,Ped,Conf,SVs,ofh):
	biallele(biallelic_dup_svm(init_dataframe(biallelic_dups),Conf.gtclf),Ped,SVs,ofh)
cdef c_male_sex_chrom_dels(male_sex_chrom_dels,Ped,Conf,SVs,ofh):
	male_sex_chrom(male_sex_chrom_svm(init_dataframe(male_sex_chrom_dels),Conf.gtclf,'delmsc'),Ped,SVs,ofh)
cdef c_male_sex_chrom_dups(male_sex_chrom_dups,Ped,Conf,SVs,ofh):
	male_sex_chrom(male_sex_chrom_svm(init_dataframe(male_sex_chrom_dups),Conf.gtclf,'dupmsc'),Ped,SVs,ofh)
cdef male_sex_chrom(df,Ped,SVs,ofh):
	for x in df.values:
		ofh.write('\t'.join(map(str,x))+'\n')
		SV = Variant(x)
		if SVs.get(SV.locus)!=None: SV=SVs[SV.locus]
		clf=str(x[20])
		SV.clf.append(clf)
		cn = format(float(x[6]),'.2f')
		svgt = int(x[13])
		if svgt==1:
			lik=np.rint(float(x[18]))
			SV.genotype(x,'0',cn,lik,False)
			if SV.ref.get(clf)==None: SV.ref[clf]=[lik]
			else: SV.ref[clf].append(lik)
		else:
			lik=np.rint(float(x[19]))
			SV.genotype(x,'1',cn,lik,False)
			if SV.alt.get(clf)==None: SV.alt[clf]=[lik]
			else: SV.alt[clf].append(lik)
		SVs[SV.locus]=SV
cdef set_biallelic_genotype(SV,chrom,data,Ped,genotype,copy_number,lik):
	skip=False
	SV.genotype(data,genotype,copy_number,lik)
	sample_id=str(data[5])
	if Ped.males.get(sample_id)==None and format_chrom(chrom) == 'chrY':
		skip=True
		SV.genotype(data,'.',copy_number,lik)
	return skip