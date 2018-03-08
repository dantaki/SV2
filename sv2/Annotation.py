#!/usr/bin/env python
from sv2_backend import format_chrom,reciprocal_overlap
from sv2Config import Config
import os,pybedtools
from pybedtools import BedTool
class Annotation():
	def __init__(self):
		self.cytoband={} # [(c,s,e)]=cytoband
		self.excluded={}
		self.repeatmasker={} #[(c,s,e)]= id,ovr
		self._1kgp={} #[(c,s,e,cl)]= id,ovr
		self.genes={}
	def check_overlap(self,svs=None,raw=None,gen=None,tmp_dir=None):
		# overlap cytobands
		try:
			for entry in svs.intersect('{}annotation_files/{}_cytoband.bed'.format(Config().resource_path(),gen),wao=True):
				locus = tokenize_sv(entry)
				if self.cytoband.get(locus)==None: self.cytoband[locus]=str(entry[3]).replace('chr','')+str(entry[6])
				else: self.cytoband[locus]=self.cytoband[locus]+','+str(entry[3]).replace('chr','')+str(entry[6])
		except pybedtools.cbedtools.MalformedBedLineError: self.cytoband[locus]='NA'
		# overlap excluded elements
		tmp_bed=tmp_dir+'tmp_anno.bed'
		svs.intersect('{}annotation_files/{}_excluded.bed.gz'.format(Config().resource_path(),gen),wao=True,output=tmp_bed)
		with open(tmp_bed,'r') as f:
			for l in f:
				entry = tuple(l.rstrip().split('\t'))
				locus, ovr = tokenize_sv(entry),int(entry[-1])
				if ovr==0: continue
				key = locus+(str(entry[6]),)
				if self.excluded.get(key)==None: self.excluded[key]=ovr
				else: self.excluded[key]+=ovr
		for locus in self.excluded: self.excluded[locus]=float(self.excluded[locus])/(int(locus[2])-int(locus[1]))
		os.remove(tmp_bed)
		# overlap repeatmasker
		svs.intersect('{}annotation_files/{}_repeatmasker.bed.gz'.format(Config().resource_path(),gen), f=0.8, F=0.8, wa=True, wb=True,output=tmp_bed)
		with open(tmp_bed,'r') as f:
			for l in f:
				entry = tuple(l.rstrip().split('\t'))
				locus, rovr=tokenize_sv(entry), reciprocal_overlap(map(int,(entry[1],entry[4],entry[2],entry[5])))
				name = '{}:{}:{}'.format(entry[7],entry[8],entry[6])
				if self.repeatmasker.get(locus)==None:
					self.repeatmasker[locus]=(name, rovr)
				elif self.repeatmasker.get(locus)!=None and rovr > self.repeatmasker[locus][1]:
					self.repeatmasker[locus]=(name, rovr)
				else: continue
		os.remove(tmp_bed)
		# overlap 1KGP phase 3 DEL/DUP
		if 'hg' in gen:
			self.load_1kgp(raw,'DEL',gen,tmp_bed)
			self.load_1kgp(raw,'DUP',gen,tmp_bed)
		# overlap genes
		genes={}
		svs.intersect('{}annotation_files/{}_genes.bed.gz'.format(Config().resource_path(),gen), wa=True,wb=True,output=tmp_bed)
		with open(tmp_bed,'r') as f:
			for l in f:
				entry = tuple(l.rstrip().split('\t'))
				locus, gene = tokenize_sv(entry), str(entry[6])
				if genes.get(locus)==None: genes[locus]=[gene]
				elif gene not in genes[locus]: genes[locus].append(gene)
				else: continue
		os.remove(tmp_bed)
		for x in genes:
			ori,exons,introns,trx,genlist = {},{},{},{},[]
			exon_total,exon_num, intron_total,intron_num,exoncnt,exontot,introncnt,introntot=0,0,0,0,0,0,0,0
			for y in genes[x]:
				gene = y.split(',')
				ori[gene[0]]=gene[1]
				if 'exon' in y:
					trx[gene[0]]=gene[2]
					exon_total = int(gene[len(gene)-1].split('/').pop())
					exon_num = int(gene[len(gene)-1].split('/').pop(0).replace('exon_',''))
					exons[(gene[0],'TOTAL')]=exon_total
					if exons.get(gene[0])==None:
						exons[gene[0]]=1
						exons[(gene[0],exon_num)]=1
					elif exons.get((gene[0],exon_num)) == None:
						exons[gene[0]]+=1
						exons[(gene[0],exon_num)]=1
					else: continue
				elif 'UTR3' in y:
					trx[gene[0]]=gene[2]
					exons[(gene[0],'UTR3')]=1
				elif 'UTR5' in y:
					trx[gene[0]]=gene[2]
					exons[(gene[0],'UTR5')]=1
				elif 'intron' in y:
					trx[gene[0]]=gene[2]
					intron_total = int(gene[len(gene)-1].split('/').pop())
					intron_num = int(gene[len(gene)-1].split('/').pop(0).replace('intron_',''))
					introns[(gene[0],'TOTAL')]=intron_total
					if introns.get(gene[0])==None:
						introns[gene[0]]=1
						introns[(gene[0],intron_num)]=1
					elif introns.get((gene[0],intron_num)) == None:
						introns[gene[0]]+=1
						introns[(gene[0],intron_num)]=1
					else: continue
				elif 'stream' in y:
					trx[gene[0]]=gene[2]
					exons[(gene[0],gene[3])]=1
				else: continue
			for y in trx:
				if ori.get(y) == None: continue
				orient=ori[y]
				exoncnt=0
				exontot=0
				if exons.get(y) != None: exoncnt=exons[y]
				if exons.get((y,'TOTAL')) != None: exontot=exons[(y,'TOTAL')]
				if orient == '+':
					if exons.get((y,'upstream_1kb')) == 1 and exons.get((y,'downstream_1kb')) != 1:
						genlist.append(','.join(map(str,(y,trx[y],'upstream_1kb'))))
					if exons.get((y,'UTR3'))!=1 and exons.get((y,'UTR5')) == 1:
						genlist.append(','.join(map(str,(y,trx[y],'UTR5'))))
					if exoncnt != 0: genlist.append(','.join(map(str,(y,trx[y],'exon_{}/{}'.format(exoncnt,exontot)))))
					if exons.get(y) == None or exoncnt != exontot:
						introncnt,introntot=0,0
						if introns.get(y) != None: introncnt=introns[y]
						if introns.get((y,'TOTAL'))!=None: introntot=introns[(y,'TOTAL')]
						if introncnt !=0: genlist.append(','.join(map(str,(y,trx[y],'intron_{}/{}'.format(introncnt,introntot)))))
					if exons.get((y,'UTR3'))==1 and exons.get((y,'UTR5')) != 1:
						genlist.append(','.join(map(str,(y,trx[y],'UTR3'))))
					if exons.get((y,'upstream_1kb')) != 1 and exons.get((y,'downstream_1kb')) == 1:
						genlist.append(','.join(map(str,(y,trx[y],'downstream_1kb'))))
				else:
					if exons.get((y,'upstream_1kb')) != 1 and exons.get((y,'downstream_1kb')) == 1:
						genlist.append(','.join(map(str,(y,trx[y],'downstream_1kb'))))
					if exons.get((y,'UTR3'))==1 and exons.get((y,'UTR5')) != 1:
						genlist.append(','.join(map(str,(y,trx[y],'UTR3'))))
					if exoncnt != 0: genlist.append(','.join(map(str,(y,trx[y],'exon_{}/{}'.format(exoncnt,exontot)))))
					if exons.get(y) == None or exoncnt != exontot:
						introncnt,introntot=0,0
						if introns.get(y) != None: introncnt=introns[y]
						if introns.get((y,'TOTAL'))!=None: introntot=introns[(y,'TOTAL')]
						if introncnt !=0: genlist.append(','.join(map(str,(y,trx[y],'intron_{}/{}'.format(introncnt,introntot)))))
					if exons.get((y,'UTR3'))!=1 and exons.get((y,'UTR5')) == 1:
						genlist.append(','.join(map(str,(y,trx[y],'UTR5'))))
					if exons.get((y,'upstream_1kb')) == 1 and exons.get((y,'downstream_1kb')) != 1:
						genlist.append(','.join(map(str,(y,trx[y],'upstream_1kb'))))
			if len(genlist)>=1 : self.genes[x]='|'.join(genlist)
	def load_1kgp(self,raw=None,svtype=None,gen=None,tmp_bed=None):
		sv = BedTool([(format_chrom(x[0]),x[1],x[2],x[3]) for x in raw if svtype in str(x[3])]).sort()
		sv.intersect('{}annotation_files/{}_1000Genomes_{}.bed'.format(Config().resource_path(),gen,svtype), f=0.8, F=0.8, wao=True,output=tmp_bed)
		with open(tmp_bed,'r') as f:
			for l in f:
				x = tuple(l.rstrip().split('\t'))
				locus = tokenize_sv(x)+(str(x[3]),)
				ovr = int(x[-1])
				if ovr==0: continue
				ovr = format(float(x[len(x)-1])/(int(x[2])-int(x[1])),'.2f')
				if self._1kgp.get(locus)==None:
					self._1kgp[locus]=(x[len(x)-2],ovr)
				elif self._1kgp.get(locus)!=None and float(ovr) > float(self._1kgp[locus][1]):
					self._1kgp[locus]=(x[len(x)-2],ovr)
				else: continue
		os.remove(tmp_bed)
def tokenize_sv(x): return (format_chrom(str(x[0])),int(x[1]),int(x[2]))