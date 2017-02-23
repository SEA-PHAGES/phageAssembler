#!/usr/bin/env python2.7

# By Dan Russell, November 2016

"""
Requires python2.7+ for the argparse module.

Arguments:
	fastq      		# Required. The file of reads in fastq format.
        genome_name		# Optional. Specify the name of this genome. Default will be the input fastq file before any special characters.
        num_reads		# Optional. Number of reads to try assembling.  Default: 80000.
        adapter_list		# Optional. Specify a file of adapters to use.  Default: Adapters.txt

"""

PATH_TO_NEWBLER = "~/454/bin/"
PATH_TO_ACEUTIL = "~/phageAssembler/AceUtil"
DEFAULT_NUM_READS = 80000
DEFAULT_READ_LENGTH = 140
DEFAULT_ADAPTER_LIST = "~/phageAssembler/Adapters/Adapters.fasta"
DEFAULT_BLAST_DATABASE = "~/phageAssembler/BLASTdbs/Actino_DB"

#from datetime import datetime
import subprocess
import argparse
import sys
import os
import re
from shutil import copy, move
from Bio import SeqIO
from Bio.SeqUtils import GC

# Make parser for options listed above
parser = argparse.ArgumentParser(description='A script to assemble phage genomes.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('fastq', help="A file of reads in fastq format.")
parser.add_argument('-g','--genome_name', help="The name of the genome you are assembling.  Default will be name of the fastq before any special characters.")
parser.add_argument('-n','--num_reads', help="The number of reads to try assembling.", type=int, default=DEFAULT_NUM_READS)
parser.add_argument('-l','--avg_read_length', help="Average read length used to calculate coverage.", type=int, default=DEFAULT_READ_LENGTH)
parser.add_argument('-a','--adapter_list', help="A fasta-formatted file of adapter sequences to be trimmed. Default file contains TruSeq Illumina adapters.", default=DEFAULT_ADAPTER_LIST)
parser.add_argument('-c','--reads_percent_cutoff', help="Contigs with more than this percentage of assembled reads will be blasted and AceUtiled.", type=int, default=2)

# Parse command-line options and save values
args = parser.parse_args()

fastq = args.fastq
if not os.path.isfile(fastq):
    sys.exit("ERROR: Couldn't find the file %s" % fastq)

if args.genome_name:
    genome_name = args.genome_name
else:
    genome_name = re.split('\W+|_', os.path.basename(args.fastq))[0]

num_reads = args.num_reads
avg_read_length = args.avg_read_length
reads_percent_cutoff = args.reads_percent_cutoff

if args.adapter_list:
    adapter_list = args.adapter_list
else:
    adapter_list = DEFAULT_ADAPTER_LIST


#Log file stuff
log_file_name = '%s_phageAssembler.log' % genome_name
log_file = open(log_file_name,'w')
def printlog(message):
    print message
    log_file.write(message + "\n")

#Input
printlog("\n***INPUT***")
printlog("\tGenome: %s" % genome_name)
printlog("\tWill assemble %s reads from the file %s." % (str(num_reads), fastq))
printlog("\tWill trim reads using adapters found in %s." % adapter_list)

def wc(filename):
    from subprocess import check_output
    return int(check_output(["wc", "-l", filename]).split()[0])

total_fastq_reads = wc(fastq)

def subsample_fastq_file(filename, number_of_reads, head_tail="head", new_file_if_all_reads=False):
    total_reads = wc(filename)/4
    if total_reads < number_of_reads:
        printlog("\tFewer reads in file than total requested.")
        if new_file_if_all_reads: 
            new_filename = '%s_All_%sk_Reads.fastq' % (genome_name, str(total_reads/1000))
            copy(filename,new_filename)
            printlog("\tCreated new file, %s, with all reads for assembly." % new_filename)
            return new_filename
        else:
            printlog("\tWill use entire original file, %s, for assembly." % filename)
            return filename
    new_filename = '%s_%skReads.fastq' % (genome_name, str(num_reads/1000))
    subprocess.call('head -n %s %s > %s' % (str(num_reads*4),filename,new_filename), shell=True)
    printlog("\tCreated new file %s with %s reads." % (new_filename, str(num_reads)))
    return new_filename

#Downsample Fastq
printlog("\n***CREATING FASTQ FOR ASSEMBLY***")
assembly_fastq = subsample_fastq_file(fastq, num_reads)
printlog("\tReads file for Newbler: %s" % assembly_fastq)

def run_newbler(projectname,adaptersfile,fastqfile):
    newbler_command = "%srunAssembly -o %s -vt %s -consed %s >> %s" % (PATH_TO_NEWBLER, projectname, adaptersfile, fastqfile, log_file_name)
#    printlog("\tWill run: %s" % newbler_command)
    subprocess.call(newbler_command, shell=True)

printlog("\n***ASSEMBLING WITH NEWBLER***")
log_file.close()
run_newbler(genome_name,adapter_list,assembly_fastq)
log_file = open(log_file_name,'a')
printlog("\tNewbler assembly complete.")

#Tidy Up Directories
printlog("\n***CLEANING UP DIRS***")
cwd = os.getcwd()
project_dir = cwd + "/%s" % genome_name
consed_dir = project_dir + "/consed" 
subprocess.call(["rm","%s/sff_dir" % consed_dir])
printlog("\tRemoved fake sff_dir.")
subprocess.call(["mkdir", "%s/solexa_dir" % consed_dir])
printlog("\tCreated new solexa_dir in consed folder to hold reads file.")
move(assembly_fastq,"%s/solexa_dir/" % consed_dir)
printlog("\tMoved %s to new solexa_dir." % assembly_fastq)
os.chdir("./%s/" % genome_name)

def parse_metrics(metricsfile):
    f = open(metricsfile,'r')
    metrics = {}
    large_done = False
    for line in f:
        if 'numberOfContigs' in line:
            if large_done:
                metrics['all_contigs'] = int(line.strip().split()[-1].split(';')[-2])
            else:
                metrics['large_contigs'] = int(line.strip().split()[-1].split(';')[-2])
                large_done = True
        if 'largestContigSize' in line:
            metrics['largest_contig'] = int(line.strip().split()[-1].split(';')[-2])
        if 'numberAssembled' in line:
            metrics['aligned_reads'] = int(line.split()[2].split(',')[0])
    f.close()
    return metrics

def parse_contigs(contigsfile,totalreads):
    f = open(contigsfile,'r')
    contigs = []
    for line in f:
        if line[0] == ">":
            contig = line[1:].rstrip().split()
            contig[1] = int(contig[1].split('=')[1])
            contig[2] = int(contig[2].split('=')[1])
            contig.append(round(100*float(contig[2])/totalreads,1))
            contig.append(round(avg_read_length*float(contig[2])/contig[1],1))
            contig.append(int(round(avg_read_length*float(total_fastq_reads/4)*contig[3]/(100*contig[1]))))
            contigs.append(contig)
    f.close()
    return contigs

#Show Assembly Results
printlog("\n***ASSEMBLY RESULTS***")
assembly_metrics = parse_metrics("%s/454NewblerMetrics.txt" % project_dir )
contig_list = parse_contigs("%s/454AllContigs.fna" % project_dir,assembly_metrics['aligned_reads'])
total_contigs = len(contig_list)
printlog("\t%s total contigs." % str(total_contigs))
contigs_to_blast = []
if total_contigs > 20:
    printlog("\tThe first 20 are (Name, length in bp, # reads, ~coverage, % reads):")
    i=0
    while i<20:
        out_string = "\t\t%s\t%s\t%s\t%s-fold\t(%s%% of assembled reads)" % (contig_list[i][0], contig_list[i][1], contig_list[i][2], contig_list[i][4], contig_list[i][3])
        if contig_list[i][3] > reads_percent_cutoff:
            out_string += "*"
            contigs_to_blast.append(contig_list[i][0])
        printlog(out_string)
        i += 1
    printlog("\t* These contigs have > %s%% of the assembled reads and will thus be blasted." % reads_percent_cutoff)
else:
    printlog("\tThey are (Name, length in bp, # reads, ~coverage, % reads):")
    for contig in contig_list:
        out_string = "\t\t%s\t%s\t%s\t%s-fold\t(%s%% of assembled reads)" % (contig[0], str(contig[1]), str(contig[2]), str(contig[4]), str(contig[3]))
        if contig[3] > reads_percent_cutoff:
            out_string += "*"
            contigs_to_blast.append(contig[0])
        printlog(out_string)
    printlog("\t* These contigs have > %s%% of the assembled reads and will thus be blasted." % reads_percent_cutoff)

def blast_contigs(seqfile,database,outfile="blast_results.txt"):
    blast_command = "blastn -db %s -query %s -out %s" % (database,seqfile,outfile)
    subprocess.call(blast_command,shell=True)
    return outfile

def biopy_blast(queryseq,database,outfile="blast_output.xml"):
    blast_command = "blastn -db %s -query %s -out %s -outfmt 5" % (database,queryseq,outfile)
    subprocess.call(blast_command,shell=True)
    result_handle = open(outfile,'r')
    from Bio.Blast import NCBIXML
    return NCBIXML.read(result_handle)

def display_blast_results(record, reblast=False):
    potential_clusters = []
    base_one = ()
    printlog("\n\n\tQuery: %s" % record.query)
    if len(record.alignments) > 10:
        printlog("\t\tFirst 10 hits:")
    else:
        printlog("\t\tAll hits:")
    if len(record.alignments) == 0:
        printlog("\t\t\tNo hits found.")
    else:
        for alignment in record.alignments[:10]:
            title_split=alignment.title.split()
            outline='\t\t\t' + ' '.join([title_split[1],title_split[2],title_split[3]]) + ' (%s bp)' % str(alignment.length)
            if title_split[-2] == 'Cluster':
                outline += ", Cluster %s" % title_split[-1]
                potential_clusters.append(title_split[-1])
            printlog(outline)
    printlog("\t\tBest hit details:")
    try:
        best = record.alignments[0]
        title_split=best.title.split()
        best_title = ' '.join([title_split[1],title_split[2],title_split[3]]) + ' (%s bp)' % str(best.length)
        printlog("\t\t\tSubject: %s" % best_title)
        i = 1
        for hsp in best.hsps[:10]:
            printlog("\t\t\t\tHit #%s" % str(i))
            printlog("\t\t\t\tLength: %s bp Score: %s E-value: %s" % (hsp.align_length,hsp.score,hsp.expect))
            printlog("\t\t\t\tIdentities: %s/%s Gaps: %s/%s." % (hsp.identities,hsp.align_length,hsp.gaps,hsp.align_length))
            if hsp.frame[1] == -1:
                strand = "Minus"
            else:
                strand = "Plus"
            printlog("\t\t\t\tStrand: Plus/%s" % strand)
            printlog("\t\t\t\tQuery: %s" % hsp.query_start)
            printlog("\t\t\t\tSbjct: %s" % hsp.sbjct_start)
            if hsp.sbjct_start == 1:
                base_one = (best_title,hsp.query_start,hsp.sbjct_start,reblast,record.query,)
            i += 1
    except:
        pass
    return (potential_clusters, base_one,)

def parse_blast(blastresults):
    best_hits = []
    with open(blastresults,'r') as in_file:
        for line in in_file:
            if (len(best_hits) <= 10) and " phage " in line:
                data = line.split()
                best_hits.append((data[2],data[-2],data[-1],))
    return best_hits

def rc_reblast(seq_file):
#    rc_out = open('%s_rc.fasta' % blast_result.query.split()[0],'w')
    out_file = '%s_rc.fasta' % seq_file.split('.')[0]
    rc_out = open(out_file,'w')
    seq = SeqIO.read(open(seq_file,'r'),'fasta')
    rc_out.write(">%s_reverse_complement\n" % seq.name)
    seq = seq.reverse_complement()
    rc_out.write(str(seq.seq))
    rc_out.close()
    return biopy_blast(out_file, DEFAULT_BLAST_DATABASE, outfile='%s_blast.xml' % out_file.split('.')[0])
    
    
 
    

#BLAST
printlog("\n***BLAST***")
printlog("\tRunning local blast of %s contig(s) against %s database..." % (str(len(contigs_to_blast)),DEFAULT_BLAST_DATABASE))
#all_contigs_file = cwd + "/%s/454AllContigs.fna" % genome_name
#blast_output = blast_contigs(all_contigs_file, DEFAULT_BLAST_DATABASE)
#print "\tBLAST complete."
#print "\tParsing BLAST results..."
#blast_results = parse_blast(blast_output)
#print "\tBest matches (Name, Score, E-value):"
#for result in blast_results:
#    print "\t\t%s\t%s\t%s" % (result[0],result[1],result[2])
all_contig_objects=SeqIO.parse(open('%s/454AllContigs.fna' % project_dir,'r'),'fasta')
blasted_contigs = []
for contig in all_contig_objects:
    if contig.id in contigs_to_blast:
        SeqIO.write(contig, '%s.fasta' % contig.id, 'fasta')
        blasted_contigs.append(biopy_blast('%s.fasta' % contig.id, DEFAULT_BLAST_DATABASE, outfile='%s_blast.xml' % contig.id))

reblasted_contigs = []
cluster_guesses = []
base_ones = []
for result in blasted_contigs:
    cg = display_blast_results(result)
    cluster_guesses.append((result.query.split()[0],cg[0],))
    if cg[1]:
        base_ones.append(cg[1])
    try:
        if result.alignments[0].hsps[0].frame[1] == -1:
            reblasted_contigs.append(rc_reblast('%s.fasta' % result.query.split()[0]))
    except:
        pass

if reblasted_contigs:
    printlog("\n\tRe-blasting %s contig(s) in reverse orientation." % str(len(reblasted_contigs)))
    for result in reblasted_contigs:
        cg = display_blast_results(result, reblast=True)
        if cg[1]:
            base_ones.append(cg[1])

def run_AceUtil(acefile,contig=None):
    try:
        outfile = acefile.rsplit('.',1)[0] + "." + str(int(acefile.rsplit('.',1)[1])+1)
    except:
        outfile = acefile + ".aceUtil"
    AceUtil_command = "java -jar %s/AceUtil.jar %s %s " % (PATH_TO_ACEUTIL, acefile, outfile)
    if contig:
        AceUtil_command += contig
    AceUtil_command += " >> %s" % (log_file_name)
    subprocess.call(AceUtil_command, shell=True)
    return outfile

#AceUtil
os.chdir("%s" % cwd)
printlog("\n***ACE UTIL***")
#Temp code until AceUtil fixed
#if len(contig_list) < 2:
#    printlog("\tSkipping AceUtil because there's only one contig.\n\tThis will be changed when AceUtil is fixed.")
#else:
#    printlog("\tRunning AceUtil...")
#    log_file.close()
#    ace_out = run_AceUtil('%s/edit_dir/454Contigs.ace.1' % consed_dir)
#    log_file = open(log_file_name,'a')
aceutil_infile = "%s/edit_dir/454Contigs.ace.1" % consed_dir
for contig in contigs_to_blast:
    printlog("\tRunning AceUtil on %s..." % contig)
    log_file.close()
    aceutil_outfile = run_AceUtil(aceutil_infile,contig=contig)
    aceutil_infile = aceutil_outfile
    log_file = open(log_file_name,'a')
printlog("\tAceUtil analysis complete.")

#Report
printlog("\n***REPORT***")

printlog("\tCluster Guess")
if cluster_guesses:
    for contig in cluster_guesses:
        gs = ', '.join(contig[1][:5])
        printlog("\t\t%s\tCluster of top hits: %s" % (contig[0], gs))
    if len(set(contig[1][:5])) == 1:
        printlog("\t\tProbable cluster: %s" %  contig[1][0])
    else:
        printlog("\t\tUnable to make single cluster guess from blast results.")
else:
    printlog("\t\tUnable to determine a likely cluster.")

printlog("\tBase One Guess")
if base_ones:
    for base_one in base_ones:
        out = "\t\tIn the blast hit to %s, query position %s matches subject position %s." % (base_one[0], str(base_one[1]), str(base_one[2]))
        out2 = "\t\tLikely Base 1 position: %s in %s" % (base_one[1], base_one[4]) 
        if base_one[3]:
            out += "  (After contig was reverse-complemented.)"
        printlog(out)
        printlog(out2)
else:
    printlog("\t\tUnable to find Base 1.")

printlog("\tGC Info")
i=0
all_contig_objects=SeqIO.parse(open('%s/454AllContigs.fna' % project_dir,'r'),'fasta')
for contig in all_contig_objects:
    if i==10:
        break
    printlog("\t\t%s\t%s %%" % (contig.id, round(GC(contig.seq),1)))
    i += 1

printlog("\tCoverage Info")
for contig in contig_list:
    printlog("\t\t%s\t%s (assembled)\t%s (estimated for entire fastq)" % (contig[0],contig[4],contig[5])) 

log_file.close()

