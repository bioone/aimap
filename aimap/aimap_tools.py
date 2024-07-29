# Copyright 2018-2024, Agriculture and Biology
# Author: sai wang
#
# This code is part of the aimap package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import os

import pandas as pd

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from . import aimap_config
import bisect

# Generate single FASTP command line
def construct_fastp_single_cmdline(outdir, outname,length, fname,threads,
                               fastp_exe=aimap_config.FASTP_DEFAULT):
    """Returns a single fastp command.

    - fastp_exe - path to FASTP executable
    """
    cmd="{0} -o {1}/{2}.fq -q 30 -l {3} -i {4} -w {5}"
    return cmd.format(fastp_exe, outdir, outname,length, fname,threads)


def construct_fastp_paired_cmdline(fname1,fname2, outname1, outname2,length,threads,
                               fastp_exe=aimap_config.FASTP_DEFAULT):
    """Returns a paired fastp command.

    - fastp_exe - path to FASTP executable
    """
    cmd="{0} -i {1} -I {2} -o {3} -O {4} -q 30 -l {5} -w {6}"
    return cmd.format(fastp_exe,fname1,fname2, outname1, outname2,length,threads)


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
    cmd = "{0} mpileup --threads {1} -d 50000 -L 50000 -f {2} {3} -o {4}/{5}.bcf -Q 0 -q 0 -A"
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

def get_anno_file(anno_file,outdir):
    file_id="genome"
    output1=open(f"{outdir}/{file_id}_gene_anno.tsv","w")
    output2=open(f"{outdir}/{file_id}_cds_anno.tsv","w")
    temp_file=open(anno_file,"r")
    for line in temp_file.readlines():
        if line[0]=="#":
            pass
        elif line.split("\t")[2]=="gene" or line.split("\t")[2]=="pseudogene":
            chr=line.split("\t")[0]
            start=line.split("\t")[3]
            end=line.split("\t")[4]
            strand=line.split("\t")[6]
            temp_dict=dict()
            attr_list=line.strip().split("\t")[8].split(";")
            for attr in attr_list:
                temp_dict[attr.split("=")[0]]=attr.split("=")[1]
            locus=temp_dict["locus_tag"]
            if "gene" in temp_dict:
                gene_name=temp_dict["gene"]
            else:
                gene_name="-"
            if "gene_biotype" in temp_dict:
                gene_bio_type=temp_dict["gene_biotype"]
            else:
                gene_bio_type="-"
            
            if "ID" in temp_dict:
                gene_ID=temp_dict["ID"]
            output1.write(f"{chr}\t{gene_ID}\t{locus}\t{gene_name}\t{gene_bio_type}\t{start}\t{end}\t{strand}\n")
        elif line.split("\t")[2]=="CDS":
            chr=line.split("\t")[0]
            start=line.split("\t")[3]
            end=line.split("\t")[4]
            strand=line.split("\t")[6]
            temp_dict=dict()
            attr_list=line.strip().split("\t")[8].split(";")
            for attr in attr_list:
                temp_dict[attr.split("=")[0]]=attr.split("=")[1]
            id="cds_"+str(start)+"_"+str(end)
            if "Parent" in temp_dict:
                locus_name=temp_dict["Parent"]
            output2.write(f"{chr}\t{locus_name}\t{id}\t{start}\t{end}\t{strand}\n")
            
    output1.close()
    output2.close()
    df1=pd.read_csv(f"{outdir}/{file_id}_gene_anno.tsv",sep="\t",names=("chr","gene_ID","locus_tag","gene_name","gene_biotype","gene_start","gene_end","gene_strand"))
    df2=pd.read_csv(f"{outdir}/{file_id}_cds_anno.tsv",sep="\t",names=("chr","gene_ID","cds_name","cds_start","cds_end","cds_strand"))
    last_df=df1.merge(df2,how="left",on=["chr","gene_ID"])
    last_df.to_csv(f"{outdir}/{file_id}_new_anno.tsv",sep="\t",index=False)


base_dict={"T":"A","A":"T","C":"G","G":"C"}

def extract_cds_from_genome(acc,seq_record,start,end,strand):
    genome_len=len(str(seq_record.seq))
    if strand=="+":
        if end<=genome_len:
            cds_seq=str(seq_record[start-1:end].seq)
        elif end>genome_len:
            cds1=str(seq_record[start-1:].seq)
            cds2=str(seq_record[0:end-genome_len].seq)
            cds_seq=cds1+cds2
    elif strand=="-":
        if end<=genome_len:
            cds_seq=str(seq_record[start-1:end].reverse_complement().seq)
        elif end>genome_len:
            cds1=str(seq_record[start-1:].reverse_complement().seq)
            cds2=str(seq_record[0:end-genome_len].reverse_complement().seq)
            cds_seq=cds2+cds1
    return Seq(cds_seq)


def get_founction_Prokaryote(infile,outdir,outname,genomefile,anno_file):
    anno_df=pd.read_csv(anno_file,sep="\t",header=0)
    orchid_dict = SeqIO.to_dict(SeqIO.parse(genomefile, "fasta"))
    result_file=open("%s/%s_full_result.txt"%(outdir,outname),"w")
    result_file.write("Accession"+"\t"+"Position"+"\t"+"Old_base"+"\t"+"New_base"+"\t"+"Raw read depth"+"\t"+"Coverage"+"\t"+"Edit_level"+
                  "\t"+"snp_coverage"+"\t"+"snp_f_coverage"+"\t"+"snp_r_coverage"+"\t"+"Gene_biotype"+"\t"+"locus_tag"+"\t"+"Gene_name"+"\t"+"Gene_strand"+"\t"+"gene_start"+
                  "\t"+"gene_end"+"\t"+"Amino acid_change"+"\t"+"codon_num"+"\n")

    with open('%s'% infile, 'r') as pileup:
        data=pileup.readlines()
        data=filter(None, data)
        for line in data:
            linsplit=line.strip().split("\t")
            accession=linsplit[0]
            position=int(linsplit[1])
            oldbase=linsplit[2]
            newbase=linsplit[3]
            temp_df=anno_df[(anno_df["chr"]==accession) & (anno_df["gene_start"] <= int(position)) & (anno_df["gene_end"]>= int(position))]
            if len(list(temp_df["locus_tag"]))==0:
                result_file.write(line.strip()+"\t"+"Intergenic region"+"\t"+" "+"\t"+" "
                                                              +"\t"+" "+"\t"+" "+"\t"+" "+"\t"+"Intergenic region"+"\t\n")
                
            elif len(list(temp_df["locus_tag"]))==1:
                gene_start=int(list(temp_df["gene_start"])[0])
                gene_end=int(list(temp_df["gene_end"])[0])
                gene_strand=list(temp_df["gene_strand"])[0]
                gene_locus=list(temp_df["locus_tag"])[0]
                gene_name=list(temp_df["gene_name"])[0]
                gene_biotype=list(temp_df["gene_biotype"])[0]
                if (gene_strand=="+") and (gene_biotype in ["-","protein_coding"]):
                    gene_start=int(list(temp_df["cds_start"])[0])
                    gene_end=int(list(temp_df["cds_end"])[0])
                    change_loc=position-gene_start
                    #old_seq=orchid_dict[accession].seq[(gene_start-1):gene_end]
                    old_seq=extract_cds_from_genome(accession,orchid_dict[accession],gene_start,gene_end,gene_strand)
                    new_seq = MutableSeq(old_seq)
                    new_seq[change_loc] = newbase
                    new_seq=Seq(new_seq)
                    old_pro=old_seq.translate()
                    new_pro=new_seq.translate()
                    codon_num=(change_loc+1)%3
                    if codon_num==0:
                        last_codon=3
                    else:
                        last_codon=codon_num
                    if old_pro==new_pro:
                        result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                              +"\t"+"+"+"\t"+str(gene_start)+"\t"+str(gene_end)+"\t"+"no change"+"\t"+str(last_codon)+"\n")
                        
                    else:
                        for n in range(len(old_pro)):
                            if old_pro[n]==new_pro[n]:
                                pass
                            else:
                                result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                              +"\t"+"+"+"\t"+str(gene_start)+"\t"+str(gene_end)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\t"+str(last_codon)+"\n")
                                continue
                elif (gene_strand=="-") and (gene_biotype in ["-","protein_coding"]):
                    gene_start=int(list(temp_df["cds_start"])[0])
                    gene_end=int(list(temp_df["cds_end"])[0])
                    change_loc=gene_end-position
                    #old_seq=orchid_dict[accession].seq[(gene_start-1):gene_end]
                    old_seq=extract_cds_from_genome(accession,orchid_dict[accession],gene_start,gene_end,gene_strand)
                    new_seq = MutableSeq(old_seq)
                    new_seq[change_loc] = base_dict[newbase]
                    new_seq=Seq(new_seq)
                    old_pro=old_seq.translate()
                    new_pro=new_seq.translate()
                    codon_num=(gene_end-position+1)%3
                    if codon_num==0:
                        last_codon=3
                    else:
                        last_codon=codon_num
                    if old_pro==new_pro:
                        result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                              +"\t"+"-"+"\t"+str(gene_start)+"\t"+str(gene_end)+"\t"+"no change"+"\t"+str(last_codon)+"\n")
                    else:
                        for n in range(len(old_pro)):
                            if old_pro[n]==new_pro[n]:
                                pass
                            else:
                                result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                              +"\t"+"-"+"\t"+str(gene_start)+"\t"+str(gene_end)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\t"+str(last_codon)+"\n")
                                continue
                else:
                    result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                        +"\t"+""+"\t"+str(gene_start)+"\t"+str(gene_end)+"\t"+""+"\t"+""+"\n")
            elif len(list(temp_df["locus_tag"]))>1:
                last_df=anno_df[(anno_df["chr"]==accession) & (anno_df["cds_start"] <= int(position)) & (anno_df["cds_end"]>= int(position))]
                cds_locus=list(last_df["cds_name"])
                if len(cds_locus)==1:
                    cds_name=cds_locus[0]
                    cds_start=int(list(last_df["cds_start"])[0])
                    cds_end=int(list(last_df["cds_end"])[0])
                    cds_strand=list(last_df["cds_strand"])[0]
                    gene_biotype=list(last_df["gene_biotype"])[0]
                    gene_locus=list(last_df["locus_tag"])[0]
                    gene_name=list(last_df["gene_name"])[0]
                    cds_len=cds_end-cds_start
                    cds_filt=(cds_len+1)%3
                    if cds_filt==0:
                        if (cds_strand=="+") and (gene_biotype in ["-","protein_coding"]):
                            change_loc=position-cds_start
                            #old_seq=orchid_dict[accession].seq[(cds_start-1):cds_end]
                            old_seq=extract_cds_from_genome(accession,orchid_dict[accession],cds_start,cds_end,cds_strand)
                            new_seq = MutableSeq(old_seq)
                            new_seq[change_loc] = newbase
                            new_seq=Seq(new_seq)
                            old_pro=old_seq.translate()
                            new_pro=new_seq.translate()
                            codon_num=(change_loc+1)%3
                            if codon_num==0:
                                last_codon=3
                            else:
                                last_codon=codon_num
                            if old_pro==new_pro:
                                result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                                  +"\t"+"+"+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+"no change"+"\t"+str(last_codon)+"\n")
                            else:
                                for n in range(len(old_pro)):
                                    if old_pro[n]==new_pro[n]:
                                        pass
                                    else:
                                        result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                                  +"\t"+"+"+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\t"+str(last_codon)+"\n")
                                        continue
                        elif (cds_strand=="-") and (gene_biotype in ["-","protein_coding"]):
                            change_loc=cds_end-position
                            #old_seq=orchid_dict[accession].seq[(cds_start-1):cds_end]
                            old_seq=extract_cds_from_genome(accession,orchid_dict[accession],cds_start,cds_end,cds_strand)
                            new_seq = MutableSeq(old_seq)
                            new_seq[change_loc] = base_dict[newbase]
                            new_seq=Seq(new_seq)
                            old_pro=old_seq.translate()
                            new_pro=new_seq.translate()
                            codon_num=(cds_end-position+1)%3
                            if codon_num==0:
                                last_codon=3
                            else:
                                last_codon=codon_num
                            if old_pro==new_pro:
                                result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                                  +"\t"+"-"+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+"no change"+"\t"+str(last_codon)+"\n")
                            else:
                                for n in range(len(old_pro)):
                                    if old_pro[n]==new_pro[n]:
                                        pass
                                    else:
                                        result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                                                                  +"\t"+"-"+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+old_pro[n]+str(n+1)+new_pro[n]+"\t"+str(last_codon)+"\n")
                        else:
                            result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                            +"\t"+""+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+"cds_filt_left"+"\t"+""+"\n")
                    else:
                        result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                            +"\t"+""+"\t"+str(cds_start)+"\t"+str(cds_end)+"\t"+"cds_len_not_3"+"\t"+""+"\n")
                else:
                    gene_start=int(list(temp_df["gene_start"])[0])
                    gene_end=int(list(temp_df["gene_end"])[0])
                    gene_strand=list(temp_df["gene_strand"])[0]
                    gene_biotype=list(temp_df["gene_biotype"])[0]
                    gene_locus=list(temp_df["locus_tag"])[0]
                    gene_name=list(temp_df["gene_name"])[0]
                    result_file.write(line.strip()+"\t"+gene_biotype+"\t"+gene_locus+"\t"+gene_name
                        +"\t"+""+"\t"+str(gene_start)+"\t"+str(gene_end)+"\t"+"no_determined_cds"+"\t"+""+"\n")
    result_file.close()

def get_A_I_result(infile,outdir,outname):
    infile=open(infile,"r")
    all_lines=infile.readlines()
    outfile=open(f"{outdir}/{outname}_A_I_result.tsv","w")
    outfile.write("Accession	Position	Old_base	New_base	Raw_read_depth	Coverage	Edit_level	snp_coverage	snp_f_coverage	snp_r_coverage	Gene_biotype	locus_tag	Gene_name	Gene_strand	gene_start	gene_end	Amino acid_change	codon_num	SRA\n")
    for line in all_lines:
        all_line_info=line.split("\t")
        if len(all_line_info)>19:
            pass
        if all_line_info[10]=="Intergenic region" and all_line_info[2]=="A" and all_line_info[3]=="G":
            outfile.write(line)
        elif all_line_info[10]=="Intergenic region" and all_line_info[2]=="T" and all_line_info[3]=="C":
            outfile.write(line)
        elif all_line_info[10]=="protein_coding" and all_line_info[2]=="T" and all_line_info[3]=="C" and all_line_info[13]=="-":
            outfile.write(line)
        elif all_line_info[10]=="protein_coding" and all_line_info[2]=="A" and all_line_info[3]=="G" and all_line_info[13]=="+":
            outfile.write(line)
        elif all_line_info[10]=="protein_coding" and all_line_info[2]=="T" and all_line_info[3]=="C" and all_line_info[13]=="" and all_line_info[16]=="no_determined_cds":
            outfile.write(line) 
        elif all_line_info[10]=="protein_coding" and all_line_info[2]=="A" and all_line_info[3]=="G" and all_line_info[13]=="" and all_line_info[16]=="no_determined_cds":
            outfile.write(line)
    outfile.close()
    
