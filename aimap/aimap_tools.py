# Copyright 2018-2019, Agriculture and Biology
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
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from . import aimap_config


# Generate single TRIM_GALORE command line
def construct_trim_galore_single_cmdline(outdir,length,fname,
                               trim_galore_exe=aimap_config.TRIM_GALORE_DEFAULT):
    """Returns a single trim_galore command.

    - trim_galore_exe - path to TRIM_GALORE executable
    """
    cmd="{0} -O {1} -q 30 --phred33 --length {2} {3}"
    return cmd.format(trim_galore_exe, outdir, length, fname)


def construct_trim_galore_paired_cmdline(outdir,length,fname1,fname2,
                               trim_galore_exe=aimap_config.TRIM_GALORE_DEFAULT):
    """Returns a paired trim_galore command.

    - trim_galore_exe - path to TRIM_GALORE executable
    """
    cmd="{0} -O {1} -q 30 --phred33 --length {2} --paired {3} {4}"
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


def construct_samtools_view_cmdline(insamfile,outdir,outname,threads,samtools_exe=aimap_config.SAMTOOLS_DEFAULT):
    cmd = "{0} view -bS {1} -@ {2} -o {3}/{4}.bam"
    return cmd.format(samtools_exe, insamfile,threads,outdir,outname)


def construct_samtools_sort_cmdline(inbamfile,outdir,outname,threads,samtools_exe=aimap_config.SAMTOOLS_DEFAULT):
    cmd = "{0} sort -@ {1} -o {2}/{3}-sorted.bam {4}"
    return cmd.format(samtools_exe,threads,outdir,outname, inbamfile)


def construct_samtools_mpileup_cmdline(sorted_bamfile,genome,outdir,outname,threads,samtools_exe=aimap_config.SAMTOOLS_DEFAULT):
    cmd = "{0} mpileup -f {1} {2} -o {3}/{4}.pileup"
    return cmd.format(samtools_exe,genome, sorted_bamfile,outdir,outname)


def convert_pileup_to_table(infile,outdir,outname):
    outfile = open("%s/%s_all_position_info.txt" % (outdir, outname),"w")
    with open('%s'% infile, 'r') as pileup:
        data=pileup.readlines()
        data=filter(None, data)
        for line in data:
            linsplit=line.split()
            if linsplit[4].count('+') >0 or linsplit[4].count('-')>0:
                continue
            if linsplit[2]=="A" or linsplit[2]=="a":
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+str(linsplit[3])+"\t"+str(int(linsplit[3])-linsplit[4].count('C')-linsplit[4].count('c')-linsplit[4].count('G')-linsplit[4].count('g')-linsplit[4].count('T')-linsplit[4].count('t'))+"\t"+str(linsplit[4].count('C')+linsplit[4].count('c'))+"\t"+str(linsplit[4].count('G')+linsplit[4].count('g'))+"\t"+str(linsplit[4].count('T')+linsplit[4].count('t'))+"\n")
            if linsplit[2]=="C" or linsplit[2]=="c":
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+str(linsplit[3])+"\t"+str(linsplit[4].count('A')+linsplit[4].count('a'))+"\t"+str(int(linsplit[3])-linsplit[4].count('A')-linsplit[4].count('a')-linsplit[4].count('G')-linsplit[4].count('g')-linsplit[4].count('T')-linsplit[4].count('t'))+"\t"+str(linsplit[4].count('G')+linsplit[4].count('g'))+"\t"+str(linsplit[4].count('T')+linsplit[4].count('t'))+"\n")
            if linsplit[2]=="G" or linsplit[2]=="g":
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+str(linsplit[3])+"\t"+str(linsplit[4].count('A')+linsplit[4].count('a'))+"\t"+str(linsplit[4].count('C')+linsplit[4].count('c'))+"\t"+str(int(linsplit[3])-linsplit[4].count('A')-linsplit[4].count('a')-linsplit[4].count('C')-linsplit[4].count('c')-linsplit[4].count('T')-linsplit[4].count('t'))+"\t"+str(linsplit[4].count('T')+linsplit[4].count('t'))+"\n")
            if linsplit[2]=="T" or linsplit[2]=="t":
                outfile.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+str(linsplit[3])+"\t"+str(linsplit[4].count('A')+linsplit[4].count('a'))+"\t"+str(linsplit[4].count('C')+linsplit[4].count('c'))+"\t"+str(linsplit[4].count('G')+linsplit[4].count('g'))+"\t"+str(int(linsplit[3])-linsplit[4].count('A')-linsplit[4].count('a')-linsplit[4].count('C')-linsplit[4].count('c')-linsplit[4].count('G')-linsplit[4].count('g'))+"\n")
    outfile.close

def get_a_i_map(infile,outdir,outname,coverage,editing_level):
    f = open("%s/%s_a_i.txt"%(outdir,outname),"w")
    with open('%s'% infile, 'r') as pileup:
        data=pileup.readlines()
        data=filter(None, data)
        for line in data:
            linsplit=line.split("\t")
            if int(linsplit[3]) <int(coverage):
                continue
            if linsplit[2]=="A" or linsplit[2]=="a":
                if int(linsplit[6])/float(linsplit[3])>float(editing_level):
                    f.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+"G"+"\t"+str(int(linsplit[6])/float(linsplit[3]))+"\t"+linsplit[3]+"\t"+linsplit[4]+"\t"+linsplit[5]+"\t"+linsplit[6]+"\t"+linsplit[7])
            if linsplit[2]=="T" or linsplit[2]=="t":
                if int(linsplit[5])/float(linsplit[3])>float(editing_level):
                    f.write(linsplit[0]+"\t"+linsplit[1]+"\t"+linsplit[2]+"\t"+"C"+"\t"+str(int(linsplit[5])/float(linsplit[3]))+"\t"+linsplit[3]+"\t"+linsplit[4]+"\t"+linsplit[5]+"\t"+linsplit[6]+"\t"+linsplit[7])
    f.close


def get_founction_Prokaryote(infile,outdir,outname,genomefile,anno_file):
    result_file=open("%s/%s_result.txt"%(outdir,outname),"w")
    result_file.write("Accession"+"\t"+"Position"+"\t"+"Old_base"+"\t"+"New_base"+"\t"+"Coverage"+"\t"+"Edit_level"+"\t"+"Gene_biotype"+"\t"+"Gene_name"+"\t"+"Product"+"\t"+"Amino acid_change"+"\n")
    db = gffutils.create_db("%s"% anno_file, ':memory:', force=True, keep_order=True,merge_strategy='merge', id_spec="ID",sort_attribute_values=True)
    with open('%s'% infile, 'r') as pileup:
        data=pileup.readlines()
        data=filter(None, data)
        for line in data:
            linsplit=line.split("\t")
            accession=linsplit[0]
            position=int(linsplit[1])
            oldbase=linsplit[2]
            newbase=linsplit[3]
            for i in db.all_features(featuretype='gene'):
                if position>=db[i.id].start and position<=db[i.id].end:
                    if db[i.id].seqid == accession and db[i.id].attributes['gene_biotype'][0]=="protein_coding":
                        change_loc=position-db[i.id].start
                        old_seq=Seq(i.sequence("%s"% genomefile, use_strand=False),IUPAC.unambiguous_dna)
                        new_seq = old_seq.tomutable()
                        new_seq[change_loc] = newbase
                        new_seq=new_seq.toseq()
                        CDS_id=[h.id for h in db.children(i.id,featuretype="CDS")][0]
                        product=db[CDS_id].attributes['product'][0]
                        if db[i.id].strand== "+":
                            old_pro=old_seq.translate()
                            new_pro=new_seq.translate()
                            for n in range(len(old_pro)):
                                if old_pro[n]==new_pro[n]:
                                    pass
                                else:
                                    result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"+str(linsplit[5])+"\t"+str(linsplit[4])+"\t"+"CDS"+"\t"+db[i.id].attributes['Name'][0]+"\t"+str(product)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\n")
                        if db[i.id].strand== "-":
                            old_seq=old_seq.reverse_complement()
                            new_seq=new_seq.reverse_complement()
                            old_pro=old_seq.translate()
                            new_pro=new_seq.translate()
                            for n in range(len(old_pro)):
                                if old_pro[n]==new_pro[n]:
                                    pass
                                else:
                                    result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"+str(linsplit[5])+"\t"+str(linsplit[4])+"\t"+"CDS"+"\t"+db[i.id].attributes['Name'][0]+"\t"+str(product)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\n")
                    if db[i.id].seqid == accession and db[i.id].attributes['gene_biotype'][0]!="protein_coding":
                        result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"+str(linsplit[5])+"\t"+str(linsplit[4])+"\t"+db[i.id].attributes['gene_biotype'][0]+"\t"+db[i.id].attributes['Name'][0]+"\t"+"-"+"\t"+"-"+"\n")
            start=[db[i.id].start for i in db.all_features(featuretype='gene')]
            end=[db[i.id].end for i in db.all_features(featuretype='gene')]
            ranges=zip(start,end)
            if any(lower <= position <= upper for (lower, upper) in ranges):
                continue
            else: 
                result_file.write(str(accession)+"\t"+str(position)+"\t"+str(oldbase)+"\t"+str(newbase)+"\t"+str(linsplit[5])+"\t"+str(linsplit[4])+"\t"+"Intergenic region"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\n")
    result_file.close
