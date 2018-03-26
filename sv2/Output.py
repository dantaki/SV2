#!/usr/bin/env python2
from sv2_backend import format_chrom
from Genotype import init_dataframe,partition_svs
from Svm import training_features
from Vcf import VCF
import datetime
def output(Structural_Variant,SVs,Ped,ids,gen,ofh,anno_flag,tmp_dir):
	Ofh = VCF(ofh)
	Ofh.init_header(datetime.date.today(),ids,Structural_Variant,gen)
	Ofh.load_genotypes(Structural_Variant,SVs,Ped,ids,gen,anno_flag,tmp_dir)
	if not ofh.endswith('.vcf'): ofh = ofh+'.vcf'
	vcf_ofh = open(ofh,'w')
	vcf_ofh.write('\n'.join(Ofh.head)+'\n')
	chroms={}
	for entry in Structural_Variant.raw:
		chroms[entry[0]]=1
		variant_id = '.'
		if Structural_Variant.variant_id.get(entry)!=None: variant_id=Structural_Variant.variant_id[entry]
		locus=(format_chrom(entry[0]),int(entry[1]),int(entry[2]),str(entry[3]))
		genotypes=[]
		if Ofh.genotypes.get(locus)==None:
			for sample_id in ids:
				gt='./.'
				if (locus[0]=='chrX' or locus[0]=='chrY') and Ped.males.get(sample_id)!=None and Structural_Variant.par[locus]==False: gt='.'
				genotypes.append(gt)
		else: genotypes=Ofh.genotypes[locus]
		genotypes = '\t'.join(genotypes)
		out=Ofh.init_row(entry)
		vcf_ofh.write('{}\t{}\t{}\t.\t{}\t{}\n'.format(entry[0],int(entry[1])+1,variant_id,out,genotypes))
	vcf_ofh.close()
def sv2_train_output(feats,Ped,gen,opre):
	biallelic_dels,biallelic_dups,male_sex_chrom_dels,male_sex_chrom_dups = partition_svs(feats,Ped,gen)
	if len(biallelic_dels)>0: training_features(init_dataframe(biallelic_dels),0,opre)
	if len(biallelic_dups)>0: training_features(init_dataframe(biallelic_dups),1,opre)
	if len(male_sex_chrom_dels)>0: training_features(init_dataframe(male_sex_chrom_dels),2,opre)
	if len(male_sex_chrom_dups)>0: training_features(init_dataframe(male_sex_chrom_dups),3,opre)