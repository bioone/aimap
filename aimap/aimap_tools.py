# Copyright 2018-2023, Agriculture and Biology
# Author: sai wang
#
# This code is part of the aimap package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import os
import shutil
import gffutils

import pandas as pd

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from . import aimap_config
import bisect

# Generate single TRIM_GALORE command line
def construct_trim_galore_single_cmdline(outdir,length,fname,
                               trim_galore_exe=aimap_config.TRIM_GALORE_DEFAULT):
    """Returns a single trim_galore command.

    - trim_galore_exe - path to TRIM_GALORE executable
    """
    cmd="{0} -O {1} -q 30 --phred33 --trim-n --length {2} {3}"
    return cmd.format(trim_galore_exe, outdir, length, fname)


def construct_trim_galore_paired_cmdline(outdir,length,fname1,fname2,
                               trim_galore_exe=aimap_config.TRIM_GALORE_DEFAULT):
    """Returns a paired trim_galore command.

    - trim_galore_exe - path to TRIM_GALORE executable
    """
    cmd="{0} -O {1} -q 30 --phred33 --trim-n --length {2} --paired {3} {4}"
    return cmd.format(trim_galore_exe, outdir, length, fname1,fname2)


def construct_bwa_index_cmdline(fname, bwa_exe=aimap_config.BWA_DEFAULT):
    cmd = "{0} index {1}"
    return cmd.format(bwa_exe, fname)


def construct_bwa_mem_single_cmdline(genome,fqfile,outdir,outname,threads,bwa_exe=aimap_config.BWA_DEFAULT):
    cmd = "{0} mem {1} {2} -t {3} -o {4}/{5}.sam"
    return cmd.format(bwa_exe, genome,fqfile,threads,outdir,outname)


def construct_bwa_mem_paired_cmdline(genome,fqfile1,fqfile2,outdir,outname,threads,bwa_exe=aimap_config.BWA_DEFAULT):
    cmd = "{0} mem {1} {2} {3} -t {4} -o {5}/{6}.sam"
    return cmd.format(bwa_exe, genome,fqfile1,fqfile2,threads,outdir,outname)

def construct_samtools_sort_cmdline(insamfile,outdir,outname,threads,samtools_exe=aimap_config.SAMTOOLS_DEFAULT):
    cmd = "{0} sort -@ {1} -o {2}/{3}-sorted.bam {4}"
    return cmd.format(samtools_exe,threads,outdir,outname, insamfile)
    
def construct_bcftools_mpileup_cmdline(sorted_bamfile,genome,outdir,outname,threads,bcftools_exe=aimap_config.BCFTOOLS_DEFAULT):
    cmd = "{0} mpileup --threads {1} -d 50000 -f {2} {3} -o {4}/{5}.bcf"
    return cmd.format(bcftools_exe,threads,genome, sorted_bamfile,outdir,outname)
    
def construct_bcftools_call_cmdline(inbcffile,outdir,outname,threads,bcftools_exe=aimap_config.BCFTOOLS_DEFAULT):
    cmd = "{0} call -m -v -o {1}/{2}.vcf -O v --threads {3} {4}"
    return cmd.format(bcftools_exe,outdir,outname,threads, inbcffile)

def construct_bcftools_view_cmdline(invcffile,outdir,outname,threads,bcftools_exe=aimap_config.BCFTOOLS_DEFAULT):
    cmd = "{0} view --types snps -o {1}/{2}-filt.vcf -O v --threads {3} {4}"
    return cmd.format(bcftools_exe, outdir,outname,threads,invcffile)
    
def construct_bcftools_query_cmdline(invcffile,outdir,outname,bcftools_exe=aimap_config.BCFTOOLS_DEFAULT):
    cmd = "{0} query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%DP4\n' {1} >{2}/{3}.tab"
    return cmd.format(bcftools_exe, invcffile,outdir,outname)


def convert_pileup_to_table(infile,outdir,outname):
    outfile = open("%s/%s_modification.txt" % (outdir, outname),"w")
    with open('%s'% infile, 'r') as pileup:
        data=pileup.readlines()
        data=filter(None, data)
        for line in data:
            linsplit=line.strip().split("\t")
            r_al=int(linsplit[-1].split(",")[2])
            l_al=int(linsplit[-1].split(",")[3])
            t_al=r_al+l_al
            r_ref=int(linsplit[-1].split(",")[0])
            l_ref=int(linsplit[-1].split(",")[1])
            h_t=r_ref+l_ref+t_al
            if (linsplit[2]=="A" or linsplit[2]=="a") and (linsplit[3]=="G" or linsplit[3]=="g"):
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+linsplit[3]+"\t"+linsplit[4]+"\t"+str(h_t)+"\t"+str(t_al/float(h_t))+"\t"+str(t_al)+"\t"+str(r_al)+"\t"+str(l_al)+"\n")
            if (linsplit[2]=="T" or linsplit[2]=="t") and ((linsplit[3]=="C" or linsplit[3]=="c")):
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+linsplit[3]+"\t"+linsplit[4]+"\t"+str(h_t)+"\t"+str(t_al/float(h_t))+"\t"+str(t_al)+"\t"+str(r_al)+"\t"+str(l_al)+"\n")
            if (linsplit[2]=="C" or linsplit[2]=="c") and (linsplit[3]=="T" or linsplit[3]=="t"):
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+linsplit[3]+"\t"+linsplit[4]+"\t"+str(h_t)+"\t"+str(t_al/float(h_t))+"\t"+str(t_al)+"\t"+str(r_al)+"\t"+str(l_al)+"\n")
            if (linsplit[2]=="G" or linsplit[2]=="g") and ((linsplit[3]=="A" or linsplit[3]=="a")):
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+linsplit[3]+"\t"+linsplit[4]+"\t"+str(h_t)+"\t"+str(t_al/float(h_t))+"\t"+str(t_al)+"\t"+str(r_al)+"\t"+str(l_al)+"\n") 
    outfile.close()

def get_founction_Prokaryote(infile,outdir,outname,genomefile,anno_file):
    result_file=open("%s/%s_full_result.txt"%(outdir,outname),"w")
    result_file.write("Accession"+"\t"+"Position"+"\t"+"Old_base"+"\t"+"New_base"+"\t"+"Raw read depth"+"\t"+"Coverage"+"\t"+"Edit_level"+
                  "\t"+"snp_coverage"+"\t"+"snp_f_coverage"+"\t"+"snp_r_coverage"+"\t"+"Gene_biotype"+"\t"+"Gene_name"+"\t"+"Gene_strand"+"\t"+"Product"+
                  "\t"+"Amino acid_change"+"\n")
    db = gffutils.create_db("%s"% anno_file, ':memory:', force=True, keep_order=True,merge_strategy='merge', id_spec="ID",sort_attribute_values=True)
    start=[db[i.id].start for i in db.all_features(featuretype='gene')]
    end=[db[i.id].end for i in db.all_features(featuretype='gene')]
    id_list=[i.id for i in db.all_features(featuretype='gene')]
    with open('%s'% infile, 'r') as pileup:
        data=pileup.readlines()
        data=filter(None, data)
        for line in data:
            linsplit=line.strip().split("\t")
            accession=linsplit[0]
            position=int(linsplit[1])
            oldbase=linsplit[2]
            newbase=linsplit[3]
            idx = bisect.bisect_left(end, position)
            if idx < len(end) and start[idx] <= position <= end[idx]:
                id=id_list[idx]
                if db[id].seqid == accession and db[id].attributes['gene_biotype'][0]=="protein_coding":
                    if db[id].strand== "+":
                        if (oldbase=="A" or oldbase=="a") and (newbase=="G" or oldbase=="g"):
                            change_loc=position-db[id].start
                            old_seq=Seq(db[id].sequence("%s"% genomefile, use_strand=False))
                            new_seq = MutableSeq(old_seq)
                            new_seq[change_loc] = newbase
                            new_seq=Seq(new_seq)
                            CDS_id=[h.id for h in db.children(id,featuretype="CDS")][0]
                            product=db[CDS_id].attributes['product'][0]
                            old_pro=old_seq.translate()
                            new_pro=new_seq.translate()
                            if old_pro==new_pro:
                                result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"
                                                              +str(linsplit[4])+"\t"+str(linsplit[5])+"\t"+str(linsplit[6])+"\t"+str(linsplit[7])
                                                              +"\t"+str(linsplit[8])+"\t"+str(linsplit[9])+"\t"+"CDS"+"\t"+db[id].attributes['Name'][0]
                                                              +"\t"+"+"+"\t"+str(product)+"\t"+"no change"+"\n")
                            else:
                                for n in range(len(old_pro)):
                                    if old_pro[n]==new_pro[n]:
                                        pass
                                    else:
                                        result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"
                                                              +str(linsplit[4])+"\t"+str(linsplit[5])+"\t"+str(linsplit[6])+"\t"+str(linsplit[7])
                                                              +"\t"+str(linsplit[8])+"\t"+str(linsplit[9])+"\t"+"CDS"+"\t"+db[id].attributes['Name'][0]
                                                              +"\t"+"+"+"\t"+str(product)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\n")
                    if db[id].strand== "-" :
                        if (oldbase=="T" or oldbase=="t") and (newbase=="C" or newbase=="c"):
                            change_loc=position-db[id].start
                            old_seq=Seq(db[id].sequence("%s"% genomefile, use_strand=False))
                            new_seq = MutableSeq(old_seq)
                            new_seq[change_loc] = newbase
                            new_seq=Seq(new_seq)
                            CDS_id=[h.id for h in db.children(id,featuretype="CDS")][0]
                            product=db[CDS_id].attributes['product'][0]
                            old_seq=old_seq.reverse_complement()
                            new_seq=new_seq.reverse_complement()
                            old_pro=old_seq.translate()
                            new_pro=new_seq.translate()
                            if old_pro==new_pro:
                                result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"
                                                              +str(linsplit[4])+"\t"+str(linsplit[5])+"\t"+str(linsplit[6])+"\t"+str(linsplit[7])
                                                              +"\t"+str(linsplit[8])+"\t"+str(linsplit[9])+"\t"+"CDS"+"\t"+db[id].attributes['Name'][0]
                                                              +"\t"+"+"+"\t"+str(product)+"\t"+"no change"+"\n")                            
                            else:
                                for n in range(len(old_pro)):
                                    if old_pro[n]==new_pro[n]:
                                        pass
                                    else:
                                        result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"
                                                              +str(linsplit[4])+"\t"+str(linsplit[5])+"\t"+str(linsplit[6])+"\t"+str(linsplit[7])
                                                              +"\t"+str(linsplit[8])+"\t"+str(linsplit[9])+"\t"+"CDS"+"\t"+db[id].attributes['Name'][0]
                                                              +"\t"+"-"+"\t"+str(product)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\n")
                if db[id].seqid == accession and db[id].attributes['gene_biotype'][0]!="protein_coding":
                    result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"+str(linsplit[4])+"\t"+str(linsplit[5])+"\t"+str(linsplit[6])+"\t"+str(linsplit[7])
                                                              +"\t"+str(linsplit[8])+"\t"+str(linsplit[9])+"\t"+db[id].attributes['gene_biotype'][0]+"\t"+db[id].attributes['Name'][0]+"\t"+" "+"\t"+" "+"\t"+" "+"\n")
            else:
                if ((oldbase=="T" or oldbase=="t") and (newbase=="C" or newbase=="c")) or ((oldbase=="A" or oldbase=="a") and (newbase=="G" or oldbase=="g")):
                    result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"+str(linsplit[4])+"\t"+str(linsplit[5])+"\t"+str(linsplit[6])+"\t"+str(linsplit[7])
                                                              +"\t"+str(linsplit[8])+"\t"+str(linsplit[9])+"\t"+"Intergenic region"+"\t"+" "+"\t"+" "+"\t"+" "+"\t"+" "+"\n")
    result_file.close()
