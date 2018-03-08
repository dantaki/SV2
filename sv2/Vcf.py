#!/usr/bin/env python
from Annotation import Annotation,tokenize_sv
from sv2_backend import format_chrom
from collections import OrderedDict
from pybedtools import BedTool
import sys
class VCF():
	def __init__(self,fh=None):
		self.fh=fh
		self.Annotations=None # Annotation object
		self.genotypes={} # [(c,s,e,cl)]=[genotypes]
		self.allele_freq={} # [locus]=[#alt alleles, total alleles]
		self.quals={} # [locus]=[ref,alt]
		self.filters={} # [locus]=[std,dnm]
		self.head=None
	def init_header(self,date=None,ids=None,Structural_Variant=None,gen=None):
		chroms,refs,contigs={},OrderedDict(),[]
		from sv2Config import Config
		with open('{}{}.genome'.format(Config().resource_path(),gen),'r') as f:
			for l in f:
				chrom, leng = l.rstrip().split('\t')
				chroms[format_chrom(chrom)]=leng
		for x in Structural_Variant.raw: 
			if chroms.get(format_chrom(x[0]))!=None: refs[x[0]]=chroms[format_chrom(x[0])]
		for chrom in refs: contigs.append('##contig=<ID={},length={}>'.format(chrom,refs[chrom]))
		self.head=[
			'##fileformat=VCFv4.1',
			'##fileDate={}'.format(date),
			'##SV2_CMD="{}"'.format(' '.join(map(str,sys.argv[:]))),
			'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">',
			'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
			'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
			'##INFO=<ID=DENOVO_FILTER,Number=1,Type=String,Description="Stringent filter status, recommended for de novo mutation discovery">',
			'##INFO=<ID=REF_GTL,Number=1,Type=Float,Description="Median Phred-adjusted REF genotype likelihood">',
			'##INFO=<ID=AF,Number=1,Type=Float,Description="Alternate allele frequency,in the range (0,1)">',
			'##INFO=<ID=CYTOBAND,Number=.,Type=String,Description="Cytoband(s) overlapping the variant">',
			'##INFO=<ID=REPEATMASKER,Number=2,Type=String,Description="Name and reciprocal overlap of RepeatMasker variant">',
			'##INFO=<ID=1000G_ID,Number=1,Type=String,Description="1000 Genomes Phase 3 integrated SV callset variant identifier">',
			'##INFO=<ID=1000G_OVERLAP,Number=1,Type=Float,Description="Overlap to 1000 Genomes Phase 3 variant, in the range (0,1)">',
			'##INFO=<ID=DESCRIPTION,Number=1,Type=String,Description="Verbose description of SV, 1-based coordinates"',
			'##INFO=<ID=GENES,Number=1,Type=String,Description="Genes overlapping the variant, pipe-separated by transcripts>"',
			'##INFO=<ID=ABPARTS,Number=1,Type=Float,Description="Overlap to antibody parts, in the range (0,1)">',
			'##INFO=<ID=CENTROMERE,Number=1,Type=Float,Description="Centromere overlap, in the range (0,1)">',
			'##INFO=<ID=GAP,Number=1,Type=Float,Description="Overlap to gaps in the reference, in the range (0,1)">',
			'##INFO=<ID=SEGDUP,Number=1,Type=Float,Description="Segmental duplication overlap, in the range (0,1)">',
			'##INFO=<ID=STR,Number=1,Type=Float,Description="Short tandem repeat overlap, in the range (0,1)">',
			'##INFO=<ID=UNMAPPABLE,Number=1,Type=Float,Description="Overlap to DAC Blacklisted Regions, in the range (0,1)">',
			'##FILTER=<ID=ABPARTS,Description="Variant overlaps to antibody parts >50%">',
			'##FILTER=<ID=CENTROMERE,Description="Variant overlaps to centromere >50%">',
			'##FILTER=<ID=GAP>,Description="Variant overlaps to reference gaps >50%">',
			'##FILTER=<ID=GENOTYPEFAIL,Description="Variant was unable to be genotyped">',
			'##FILTER=<ID=NOALT,Description="No alternate allele detected">',
			'##FILTER=<ID=SEGDUP,Description="Variant overlaps to segmental duplications >50%">',
			'##FILTER=<ID=STR,Description="Variant overlaps to short tandem repeats >50%">',
			'##FILTER=<ID=UNMAPPABLE,Description="Variant overlaps to DAC Blacklisted Regions regions >50%">',
			'##FILTER=<ID=FAIL,Description="Variant failed standard filters">',
			'##FILTER=<ID=PASS,Description="Variant passed standard filters">',
			'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
			'##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number estimate">',
			'##FORMAT=<ID=PE,Number=1,Type=Float,Description="Normalized discordant paired-end count">',
			'##FORMAT=<ID=SR,Number=1,Type=Float,Description="Normalized split-read count">',
			'##FORMAT=<ID=SC,Number=1,Type=Float,Description="SNV normalized coverage">',
			'##FORMAT=<ID=NS,Number=1,Type=Integer,Description="Number of SNVs within locus">',
			'##FORMAT=<ID=HA,Number=1,Type=Float,Description="Heterozygous allele ratio">',
			'##FORMAT=<ID=NH,Number=1,Type=Integer,Description="Number of heterozygous SNVs">',
			'##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled genotype likelihood">',
			'##FORMAT=<ID=GL,Number=2|3,Type=Float,Description="Phred-scaled genotype likelihoods in the order, REF:(0/0), HET:(0/1), HOM:(1/1)">',
			'##ALT=<ID=DEL,Description="Deletion, if 80% reciprocal overlap with RepeatMasker element, the class, name, and family are given separated by colons">',
			'##ALT=<ID=DUP,Description="Duplication, if 80% reciprocal overlap with RepeatMasker element, the class, name, and family are given separated by colons">',
			'{}'.format('\n'.join(contigs)),
			'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format('\t'.join(ids)),
		]
	def init_row(self,data=None):
		svtype, _1kgp_id, _1kgp_ovr,mei_id,mei_ovr= str(data[3]),'NA','NA','NA',0
		gene,cytoband,median_ref,qual,allele,std_flt,dnm_flt = 'intergenic','NA','NA','NA','NA','PASS','NA'
		locus = tokenize_sv(data)
		svlen=locus[2]-locus[1]
		if self.quals.get(locus+(svtype,))!=None: median_ref,qual = self.quals[locus+(svtype,)]
		if self.allele_freq.get(locus+(svtype,))!=None: allele=float(format(float(self.allele_freq[locus+(svtype,)][0])/float(self.allele_freq[locus+(svtype,)][1]),'.3f'))
		if self.filters.get(locus+(svtype,))!=None: std_flt,dnm_flt=self.filters[locus+(svtype,)]
		description='{}:{}-{}_{}'.format(data[0],int(data[1])+1,data[2],svtype.replace('DEL','deletion').replace('DUP','duplication'))
		fail=[]
		if self.Annotations!=None:
			if self.Annotations.repeatmasker.get(locus)!=None:
				mei_id,mei_ovr= self.Annotations.repeatmasker[locus]
				mei_ovr= format(mei_ovr,'.2f')
				svtype=svtype+':'+mei_id
			if self.Annotations.cytoband.get(locus)!=None: cytoband=self.Annotations.cytoband[locus]
			abpart_ovr = check_dict(self.Annotations.excluded,locus+('abparts',))
			centromere_ovr = check_dict(self.Annotations.excluded,locus+('centromere',))
			gap_ovr = check_dict(self.Annotations.excluded,locus+('gap',))
			segdup_ovr = check_dict(self.Annotations.excluded,locus+('segDup',))
			str_ovr = check_dict(self.Annotations.excluded,locus+('STR',))
			unmap_ovr = check_dict(self.Annotations.excluded,locus+('unmappable',))
			if self.Annotations._1kgp.get(locus+(str(data[3]),))!=None: _1kgp_id,_1kgp_ovr=self.Annotations._1kgp[locus+(str(data[3]),)]
			if self.Annotations.genes.get(locus)!=None: gene=self.Annotations.genes[locus]
			# filters
			if float(abpart_ovr) >= 0.5: fail.append('ABPARTS')
			if float(centromere_ovr) >= 0.5: fail.append('CENTROMERE')
			if float(gap_ovr)>=0.5: fail.append('GAP')
			if float(segdup_ovr) >= 0.5: fail.append('SEGDUP')
			if float(str_ovr) >= 0.5: fail.append('STR')
			if float(unmap_ovr) >= 0.5: fail.append('UNMAPPABLE')
		if median_ref=='NA' and qual=='NA': fail.append('GENOTYPEFAIL')
		if allele==0: fail.append('NOALT')
		if std_flt=='FAIL': fail.append('FAIL')
		elif std_flt=='PASS' and ('GENOTYPEFAIL' not in fail and 'NOALT' not in fail): fail.append('PASS')
		if len(fail) > 0: filt = ','.join(fail)
		if std_flt=='FAIL':dnm_flt='FAIL'
		if ('NOALT' in filt or 'GENOTYPEFAIL' in filt) and dnm_flt=='PASS':dnm_flt='NA' 
		if self.Annotations!=None:
			info = 'END={};SVTYPE={};SVLEN={};DENOVO_FILTER={};REF_GTL={};AF={};CYTOBAND={};REPEATMASKER={},{};1000G_ID={};1000G_OVERLAP={};DESCRIPTION={};GENES={};ABPARTS={};CENTROMERE={};GAP={};SEGDUP={};STR={};UNMAPPABLE={}'.format(data[2],svtype,svlen,dnm_flt,median_ref,allele,cytoband,mei_id,mei_ovr,_1kgp_id,_1kgp_ovr,description,gene,abpart_ovr,centromere_ovr,gap_ovr,segdup_ovr,str_ovr,unmap_ovr)
		else:
			info = 'END={};SVTYPE={};SVLEN={};DENOVO_FILTER={};REF_GTL={};AF={};'.format(data[2],svtype,svlen,dnm_flt,median_ref,allele)
		return '<{}>\t{}\t{}\t{}\tGT:CN:PE:SR:SC:NS:HA:NH:SQ:GL'.format(svtype,qual,filt,info)
	def load_genotypes(self,Structural_Variant=None,SVs=None,Ped=None,ids=None,gen=None,no_anno=None,tmp_dir=None):
		svs=BedTool(list(set([(format_chrom(x[0]),x[1],x[2]) for x in Structural_Variant.raw]))).sort()
		if no_anno==False:
			Annot = Annotation()
			Annot.check_overlap(svs,Structural_Variant.raw,gen,tmp_dir)
			self.Annotations=Annot
		for locus in SVs:
			Variant = SVs[locus]
			self.quals[locus]=Variant.med_ref,Variant.med_alt
			self.filters[locus]=Variant.standard_filter,Variant.denovo_filter
			for sample_id in ids:
				gt='./.'
				if (locus[0]=='chrX' or locus[0]=='chrY') and Ped.males.get(sample_id)!=None and Structural_Variant.par[locus]==False: gt='.'
				if Variant.gt.get(locus+(sample_id,))!=None:
					gt=Variant.gt[locus+(sample_id,)]
					for allele in gt.split(':').pop(0).split('/'):
						if allele=='.': continue
						if self.allele_freq.get(locus)==None: self.allele_freq[locus]=[int(allele),1]
						else: self.allele_freq[locus]=[self.allele_freq[locus][0]+int(allele),self.allele_freq[locus][1]+1]
				if self.genotypes.get(locus)==None: self.genotypes[locus]=[gt]
				else: self.genotypes[locus].append(gt)
def check_dict(d,key):
	if d.get(key)!=None: return float(format(d[key],'.3f'))
	else: return 0