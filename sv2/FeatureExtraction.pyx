#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
from sv2_backend import accepted_chrom,format_chrom,match_chrom_key,match_chrom_prefix
from sv2Config import Config
from Preprocess import Preprocess
from Snv import extract_snv_features
from SV import skip_sv
import os,pysam,sys
import numpy as np
cimport numpy as np
np.import_array()
from pybedtools import BedTool
def extract_feats(Bam,sv,prefh,out,gen,pcrfree,legacyM,Conf,tmp_dir):
	master_sv = filter_sv(sv,Bam,gen,tmp_dir)
	Conf.gc_norm(pcrfree)
	sv_gc=gc_sv(sv,master_sv,gen,Conf)
	Pre = Preprocess(prefh)
	snv_feats,het_feats=extract_snv_features(Bam,sv,master_sv,Pre,gen)
	Itr = pysam.AlignmentFile(Bam.fh,'r{}'.format(Bam.char))
	insert_size, insert_mad = Pre.insert_size[Bam.id],Pre.insert_mad[Bam.id]
	ofh = open(out,'w')
	ofh.write('\t'.join(('#chr','start','end','type','size','id','coverage','coverage_GCcorrected','discordant_ratio','split_ratio','snv_coverage','heterozygous_allele_ratio','snvs','het_snvs'))+'\n')
	for call in sv:
		(c,s,e,cl) = call
		if sv_gc.get(call)==None or master_sv.get(call)==None: continue
		snv_depth,snvs,had,hets=float('nan'),0,float('nan'),0
		if snv_feats.get(call)!=None: snv_depth,snvs=snv_feats[call]
		if het_feats.get(call)!=None: had,hets=het_feats[call]
		chrom_key = match_chrom_key(call,Bam,gen)
		if Pre.chr_cov.get((Bam.id,chrom_key))==None: continue
		GC_content = sv_gc[call]
		sv_span,flank_span,windows = master_sv[call]
		size = int(e)-int(s)
		ci = Pre.insert_size[Bam.id]+(5*Pre.insert_mad[Bam.id])
		discordant,split,concordant = discordant_split_sv(flank_span,Itr,size,ci,windows,Bam.chr_flag,legacyM)
		if size > 1000:
			"""read count coverage estimation"""
			(Aln_count,bp_span) = count_reads(sv_span,Itr,ci,Bam.chr_flag,legacyM)
			cov=float('nan')
			if Pre.chr_cov[(Bam.id,chrom_key)] != 0: cov = normalize_coverage(Aln_count,bp_span,Pre.chr_cov[(Bam.id,chrom_key)],Pre.read_len[Bam.id])
			gc_norm_factor = 1.0
			if Conf.gc_dict.get(('RC',GC_content))==None:
				if Conf.gc_dict.get(('DOC',GC_content)) !=None: gc_norm_factor=Conf.gc_dict[('DOC',GC_content)]
			else: gc_norm_factor = Conf.gc_dict[('RC',GC_content)]
			gc_cov = cov * gc_norm_factor
		else:
			"""median depth of coverage estimation"""
			median_doc = depth_of_coverage(sv_span,Bam.fh,Bam.chr_flag,Pre.read_len[Bam.id])
			cov=float('nan')
			if Pre.chr_cov[(Bam.id,chrom_key)]!= 0: cov = median_doc/Pre.chr_cov[(Bam.id,chrom_key)]
			gc_norm_factor=1.0
			if Conf.gc_dict.get(('DOC',GC_content)) != None: gc_norm_factor=Conf.gc_dict[('DOC',GC_content)]
			gc_cov = cov * gc_norm_factor
		if float(concordant) == 0.0:
			discordant_ratio = str(round(float(discordant)/1.0,3))
			split_ratio = str(round(float(split)/1.0,3))
		else:
			discordant_ratio = str(round(float(discordant)/float(concordant),3))
			split_ratio = str(round(float(split)/float(concordant),3))
		ofh.write('\t'.join(map(str,(c,s,e,cl,size,Bam.id,round(float(cov),3),round(float(gc_cov),3),discordant_ratio,split_ratio,round(float(snv_depth),3),round(float(had),3),snvs,hets)))+'\n')
	Itr.close()
	ofh.close()
cdef append_list2dict(a):
	d={}
	for i in a:
		(c,s,e,cl,tag) = i
		if d.get(tag) == None: d[tag]=[(c,s,e)]
		else: d[tag].append((c,s,e))
	return d
cdef count_reads(sv_list,Itr,ci,chr_flag,legacyM):
	"""count reads within SV span for coverage estimation"""
	cdef unsigned int Aln_count
	cdef unsigned int bp_span
	Aln_count,bp_span=0,0
	for (c,s,e) in sv_list:
		c =  match_chrom_prefix(c,chr_flag)
		region= str('{}:{}-{}').format(c,int(s)+1,e)
		bp_span+=int(e)-int(s)
		"""count each Aln within the span"""
		for Aln in Itr.fetch(region=region):
			if (Aln.is_reverse == Aln.mate_is_reverse or Aln.is_proper_pair == False or Aln.is_qcfail == True or Aln.mapping_quality < 10
			   or ( (Aln.is_secondary and legacyM == True) or (Aln.is_supplementary and legacyM == False) )
			   or Aln.is_unmapped or Aln.mate_is_unmapped or Aln.is_duplicate or abs(Aln.tlen) >= ci or Aln.tid != Aln.rnext):
				continue
			Aln_count+=1
	return(Aln_count,bp_span)
cdef depth_of_coverage(sv_list,bamfh,chr_flag,Aln_length):
	"""return median depth of coverage for SV <= 1kb"""
	pos_doc={}
	for (c,s,e) in sv_list:
		c =  match_chrom_prefix(c,chr_flag)
		region= str('{}:{}-{}').format(c,int(s)+1,e)
		depth_result = pysam.depth("-a", "-Q" "40", "-r", region, "-l", str(Aln_length-10), bamfh)
		str_flag=0
		if isinstance(depth_result,str):
			depth_result = depth_result.split('\n')
			str_flag=1
		for x in depth_result:
			r = x.rstrip('\n').split('\t')
			if str_flag == 1:
				if len(r)!=3: continue
				pos_doc[int(r[1])]=int(r[2])
			else: pos_doc[int(r[1])]=int(r[2])
	if len(pos_doc) !=0:
		temp=[]
		for x in pos_doc: temp.append(float(pos_doc[x]))
		return np.median(temp)
	else: return 0
cdef discordant_split_sv(flank_list,Itr,size,ci,windows,chr_flag,legacyM):
	"""count discordant paired ends and split Alns"""
	cdef unsigned int discordant_count
	cdef unsigned int split_count
	cdef unsigned int concordant_count	
	discordant_count,split_count,concordant_count=0,0,0
	if flank_list==None: return (0,0,0)
	else:
		for (c,s,e) in flank_list:
			c =  match_chrom_prefix(c,chr_flag)
			region= str('{}:{}-{}').format(c,int(s)+1,e)
			for Aln in Itr.fetch(region=region,until_eof=True):
				if (Aln.is_qcfail or Aln.is_unmapped or Aln.mate_is_unmapped or Aln.tid != Aln.rnext or Aln.is_reverse == Aln.mate_is_reverse or Aln.is_duplicate): continue
				else:
					mate_pos = Aln.pnext
					if is_discordant(Aln,windows,ci,mate_pos) == True: discordant_count+=1
					if is_split(Aln,windows,c,legacyM)  == True: split_count+=1
					if (Aln.is_proper_pair and Aln.mapping_quality >= 10 and ( (not Aln.is_secondary and legacyM == True) or (not Aln.is_supplementary and legacyM == False) ) and abs(Aln.tlen) < ci):
						concordant_count+=1
					else: continue
		return(discordant_count,split_count,concordant_count)
cdef filter_sv(sv,Bam,gen,tmp_dir):
	#returns master_sv[sv_dict[int(i)]]=([sv_span],[flank_span,windows])
	cdef unsigned int tag=0
	annot_sv, sv_dict, start_flank, end_flank=[], {}, [],[]
	for r in sv:
		sv_dict[tag]=r
		(c,s,e,cl) = r
		chrom=format_chrom(c)
		c = match_chrom_prefix(c,Bam.chr_flag) # chromosome matches ones in the BAM file
		if c not in accepted_chrom(Bam.chr_flag,gen):
			sys.stderr.write('WARNING: Skipping {}. SV2 does not genotype SVs on contig {}\n'.format(r,c))
			continue
		if skip_sv(c,s,e,cl,Bam.refs[c])==True: continue
		annot_sv.append((chrom,s,e,cl,tag))
		start_flank.append((chrom,s,int(s)+1,cl,tag))
		end_flank.append((chrom,e,int(e)+1,cl,tag))
		tag+=1
	"""add 500bp to each start and end position"""
	_prefix='chr'
	if Bam.chr_flag==True: _prefix='' # if reference lengths are prefixed with chr, pass an empty string
	tmp_genome = tmp_dir+'/{}.genome'.format(Bam.id)
	tmpfh = open(tmp_genome,'w')
	for ref in Bam.refs: tmpfh.write('{}{}\t{}\n'.format(_prefix,ref,Bam.refs[ref]))
	tmpfh.close()
	start_flank500bp=BedTool(start_flank).slop(b=500,g=tmp_genome)
	end_flank500bp=BedTool(end_flank).slop(b=500,g=tmp_genome)
	start_flank_dict=dict((i[4],(i[1],i[2])) for i in start_flank500bp)
	end_flank_dict=dict((i[4],(i[1],i[2])) for i in end_flank500bp)
	os.remove(tmp_genome)
	"""mask regions """
	masked_sv = [(x[0],x[1],x[2],x[3],x[4]) for x in BedTool(annot_sv).subtract(BedTool('{}{}_excluded.bed.gz'.format(Config().resource_path(),gen)))]
	masked_start_flank = [(x[0],x[1],x[2],x[3],x[4]) for x in start_flank500bp.subtract(BedTool('{}{}_excluded.bed.gz'.format(Config().resource_path(),gen)))]
	masked_end_flank= [(x[0],x[1],x[2],x[3],x[4]) for x in end_flank500bp.subtract(BedTool('{}{}_excluded.bed.gz'.format(Config().resource_path(),gen)))]
	"""take the union of the two masked call sets"""
	if len(masked_sv) == 0:
		print 'FATAL ERROR: SVs completely overlap with excluded regions'
		sys.stderr.write('FATAL ERROR: SVs completely overlap with excluded regions\n')
		sys.exit(1)
	sv_tag = list(zip(*masked_sv)[4])
	passed_tag = list(set(sv_tag))
	masked_sv_dict=append_list2dict(masked_sv)
	masked_start_flank_dict=append_list2dict(masked_start_flank)
	masked_end_flank_dict=append_list2dict(masked_end_flank)
	master_sv={}
	for i in passed_tag:
		windows = (start_flank_dict[i],end_flank_dict[i])
		sv_span = masked_sv_dict[i]
		flank_span=None
		if masked_start_flank_dict.get(i) != None and masked_end_flank_dict.get(i) != None:
			flank_span = masked_start_flank_dict[i] + masked_end_flank_dict[i]
		master_sv[sv_dict[int(i)]]=(sv_span,flank_span,windows)
	return master_sv
cdef GC(nuc):
	cdef int A=0
	cdef int C=0
	cdef int G=0
	cdef int T=0
	cdef int N=0
	for x in nuc:
		A+=int(x[5])
		C+=int(x[6])
		G+=int(x[7])
		T+=int(x[8])
		N+=int(x[9])
	return ((G+C)*100)/(A+C+G+T+N)
cdef gc_sv(sv,master_sv,gen,Conf):
	sv_gc={}
	fasta = Conf.fasta(gen)
	for call in sv:
		if master_sv.get(call)==None: continue
		svs = master_sv[call][0]
		if Conf.fasta_chr_flag == False: svs = [tuple(map(lambda y: str.replace(str(y),'chr',''),x)) for x in svs]
		GC_content = int(5 * round(float(GC(BedTool(svs).nucleotide_content(fi=fasta)))/5))
		sv_gc[call]=GC_content
	return sv_gc
cdef is_discordant(Aln,windows,ci,matepos):
	"""returns True if Aln is discordant"""
	((s1,e1),(s2,e2)) = windows
	if abs(Aln.tlen) >= ci:
		"""insert size is greater than 5MAD from median"""
		if (int(s1) <= Aln.pos <= int(e1) and int(s2) <= int(matepos) <= int(e2)) or (int(s2) <= Aln.pos <= int(e2) and int(s1) <= int(matepos) <= int(e1)):
			return True
		else:
			return False
cdef is_split(Aln,windows,c,legacyM):
	"""returns True if Aln is split"""
	((s1,e1),(s2,e2)) = windows
	second_align=[]
	if legacyM == True:
		if Aln.is_secondary == True: second_align = Aln.get_tag("SA").split(',')
	else:
		if Aln.is_supplementary == True: second_align = Aln.get_tag("SA").split(',')
	if len(second_align) == 0: return False
	else:
		if second_align[0] != c: return False
		else:
			second_align[1] = int(second_align[1])
			if ( ((int(s1) <= Aln.pos <= int(e1)) and (int(s2) <= second_align[1]-1 <= int(e2)))
			 or ((int(s2) <= Aln.pos <= int(e2)) and (int(s1) <= second_align[1]-1 <= int(e1)))):
				"""supplementary/secondary alignment must be near the opposite breakpoint of the primary alignment"""
				return True
			else:
				return False
cdef normalize_coverage(double Aln_count, double span, double chr_cov, double Aln_length): return (((Aln_count/span)*Aln_length)/chr_cov)