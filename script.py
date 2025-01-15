# requires pandas, pyBigWig, Biopython
import pandas
import pyBigWig
from Bio.Seq import Seq

# define file paths for inputs and output
genome_file = "/Users/santiagochang/Desktop/JSB/WBcel235/genome.fna"
gtf_file = "/Users/santiagochang/Desktop/JSB/WBcel235/genomic.gtf"
phylop_file = "/Users/santiagochang/Desktop/JSB/WBcel235/phylop.bw"
output_path = "/Users/santiagochang/Desktop/JSB/WBcel235/output.tsv"

# input: fasta file
# output: chromosome ID to nucleotide sequence and chromosome ID to chromosome label dictionaries
def construct_chromosome_dictionaries(fasta_file):
    # initialize dictionaries and current chromosome variable as empty
    chrom_id_to_nuc_seq = {}
    chrom_id_to_chrom_label = {}
    curr_chrom = ''
    # loop through lines in fasta file
    for line in open(fasta_file, 'r'):
        # if line is a header line
        if line[0] == '>':
            # extract chromosome ID, set current chromosome to chromosome ID
            curr_chrom = (line[1:].split())[0]
            # add ID-list (empty) and ID-label pairs to nucleotide and label dictionaries respectively
            chrom_id_to_nuc_seq[curr_chrom] = []
            chrom_id_to_chrom_label[curr_chrom] = line.split()[-1]
            # specific to C elegans fafsa file, mitochondrial gene header ends in 'genome'
            if line.split()[-1] == 'genome':
                chrom_id_to_chrom_label[curr_chrom] = 'M'
        # if line is a sequence line
        else:
            # append line to current chromosome's list
            chrom_id_to_nuc_seq[curr_chrom].append(line.strip())
    # for each chromosome, concatenate strings in chromosome's list of nucleotide lines, map chromosome ID to new string
    for chrom_id in chrom_id_to_nuc_seq:
        chrom_id_to_nuc_seq[chrom_id] = ''.join(chrom_id_to_nuc_seq[chrom_id])
    # return chromosome ID to nucleotide sequence and chromosome ID to chromosome label dictionaries
    return chrom_id_to_nuc_seq, chrom_id_to_chrom_label

# input: chromosome ID, chromosome ID to nucleotide dictionary, start index, stop index, strand direction
# output: string of nucleotide sequence
def get_nucleotide_sequence(chromosome_id, nucleotide_dict, start, stop, strand):
    # declare sequence (Biopython's Seq) as string defined by chromosome and indices (string: 0-indexed; gtf: 1-indexed)
    seq = Seq(nucleotide_dict[chromosome_id][int(start)-1:int(stop)])
    # if reverse strand, set sequence to its own reverse complement
    if strand == '-':
        seq = seq.reverse_complement()
    # return sequence as string type
    return str(seq)

# input: phyloP data, chromosome label, start index, stop index, strand direction
# output: list of phyloP scores
def get_phylop_scores(phylop_data, chromosome_label, start, stop, strand):
    # extract scores defined by chromosome and indices with pyBigWig function (phyloP: 0-indexed; gtf: 1-indexed)
    scores = phylop_data.values(chromosome_label, start-1, stop)
    # if reverse strand, reverse the list of scores
    if strand == '-':
        scores.reverse()
    # return scores as list type
    return scores

# input: list and separator string
# output: string of concatenated list elements separated by separator
def list_to_string(l, separator):
    string = ''
    if len(l) != 0:
        string = separator.join(list(map(str, l)))
    return string

# input: none
# output: empty variables corresponding to column variable type in processed data frame
def reset():
    return '', '', [], [], [], [], [], [], '', '', [], '', '', '', False

# input: gtf data frame, nucleotide sequence and chromosome label dictionaries, phyloP data
# output: processed data frame with information on IDs, nucleotide sequences, indices, phyloP scores
def construct_data_frame(gtf_data_frame, nucleotide_dictionary, chromosome_dictionary, phylop_data):
    # initialize empty processed rows list, progress tracker, relevant variables
    processed_rows = []
    progress = 0
    protein_coding_gene, isoforms, isoform_index, min_start = False, [], -1, 0
    (chrom_id, trans_id,
     cds, cdna, cds_starts, cds_ends, cdna_starts, cdna_ends,
     chrom_label, strand,
     phylop_scores, phylop_range,
     cds_start_in_mrna, cds_length,
     mrna_transcript) = reset()
    # loop through rows in gtf data frame
    for index in range(len(gtf_data_frame)):
        # print progress in 1% intervals
        if int(index/(len(gtf_data_frame)-1)*100) == progress:
            print(f'Construction is {progress}% complete.')
            progress+=1
        # initialize row to respective row in gtf based on index
        row = gtf_data_frame.iloc[index]
        # if row's feature is a protein coding gene, reset the isoform and minimum start information
        if row['feature'] == 'gene' and 'protein_coding' in row['attribute']:
            protein_coding_gene = True
            isoforms = []
            isoform_index = -1
            min_start = 0
        # if row's feature is a mRNA transcript, cDNA (and potentially CDS) strand data described in following lines
        elif row['feature'] == 'transcript' and 'gbkey "mRNA"' in row['attribute'] and protein_coding_gene:
            # extract chromosome ID, transcript ID, chromosome label, strand direction from row and assign to variables
            chrom_id = row['seqID']
            # if you prefer gene_id over transcript_id, then replace ".split('; ')[1]" with ".split(';')[0]"
            trans_id = row['attribute'].split('; ')[1].split()[1].strip('"')
            chrom_label = chromosome_dictionary[chrom_id]
            strand = row['strand']
            mrna_transcript = True
        # if row's feature is an exon coming from a protein coding gene and mRNA transcript
        elif row['feature'] == 'exon' and protein_coding_gene and mrna_transcript:
            # extract start index, stop index from row and append to cDNA start and cDNA end indices lists
            start = int(row['start'])
            end = int(row['end'])
            cdna_starts.append(start)
            cdna_ends.append(end)
            # append nucleotide sequence defined by chromosome, start index, stop index, strand direction to cDNA list
            cdna.append(get_nucleotide_sequence(chrom_id, nucleotide_dictionary, start, end, strand).strip())
            # get phyloP scores defined by chromosome, start index, stop index and add to phyloP list
            phylop_scores.extend(get_phylop_scores(phylop_data, 'chr' + chrom_label, start, end, strand))
        # if row's feature is a transcript ID or stop codon pertaining to CDS
        elif protein_coding_gene and mrna_transcript:
            # extract start index, stop index from row and append to CDS start and CDS end indices lists
            start = int(row['start'])
            end = int(row['end'])
            cds_starts.append(start)
            cds_ends.append(end)
            # append nucleotide sequence defined by chromosome, start index, stop index, strand direction to CDS list
            cds.append(get_nucleotide_sequence(chrom_id, nucleotide_dictionary, start, end, strand).strip())
            # if you prefer gene_id over transcript_id, then add the following code block
            """
            if len(cds) == 1 and trans_id[-1] != row['feature'][-1]:
                trans_id += (row['feature'][-1])
            """
        # if last row or next row's feature is transcript and information recorded, cDNA and CDS information complete
        if ((index+1 == len(gtf_data_frame) or gtf_data_frame.iloc[index+1]['feature'] == 'transcript')
                and protein_coding_gene and mrna_transcript):
            # if CDS information exists
            if len(cds) != 0:
                # test for minimum start and reassign representative isoform index as necessary
                temp = int(cds_starts[0])
                if strand == '-':
                    temp *= -1
                if temp < min_start:
                    min_start = temp
                    isoform_index = len(isoforms)
                # convert all lists to strings formatted as per example
                cds = list_to_string(cds, '')
                cdna = list_to_string(cdna, '')
                cds_starts = list_to_string(cds_starts, ',')
                cds_ends = list_to_string(cds_ends, ',')
                cdna_starts = list_to_string(cdna_starts, ',')
                cdna_ends = list_to_string(cdna_ends, ',')
                phylop_scores = (list_to_string(phylop_scores, ',').replace('[', '').replace(']', '')
                                                                   .replace(' ', '').replace('nan', 'NA'))
                # calculate proportion of available phyloP scores for bases in strand
                phylop_range = 1 - (phylop_scores.split(',').count('NA'))/(len(phylop_scores.split(',')))
                # if start index for CDS in mRNA not defined, find CDS start index in mRNA and length
                if cds_start_in_mrna == '':
                    for starting_position in range(0, len(cdna)):
                        if cdna[starting_position:len(cds)+starting_position] == cds:
                            cds_start_in_mrna = starting_position
                            cds_length = len(cds)
                            break
                # append ordered values to isoforms list, with order defined below
                isoforms.append([trans_id,
                                 cds,
                                 cdna,
                                 cds_starts,
                                 cds_ends,
                                 cdna_starts,
                                 cdna_ends,
                                 chrom_label,
                                 strand,
                                 phylop_scores,
                                 phylop_range,
                                 cds_start_in_mrna,
                                 cds_length])
            # reset variables to prepare for next transcript
            (chrom_id, trans_id,
             cds, cdna, cds_starts, cds_ends, cdna_starts, cdna_ends,
             chrom_label, strand,
             phylop_scores, phylop_range,
             cds_start_in_mrna, cds_length,
             mrna_transcript) = reset()
        # if last row or next row's feature is gene, append representative isoform to processed rows and reset values
        if ((index + 1 == len(gtf_data_frame) or gtf_data_frame.iloc[index + 1]['feature'] == 'gene')
                and protein_coding_gene):
            if len(isoforms) != 0:
                processed_rows.append(isoforms[isoform_index])
            protein_coding_gene, isoforms, isoform_index, min_start = False, [], -1, 0
    # return processed data frame generated from processed rows list
    return pandas.DataFrame(processed_rows, columns = ['transcript_id',
                                                       'cds',
                                                       'cdna',
                                                       'cds_starts',
                                                       'cds_ends',
                                                       'cdna_starts',
                                                       'cdna_ends',
                                                       'chromosome',
                                                       'strand',
                                                       'phyloP_scores',
                                                       'phyloP_range',
                                                       'cds_start_in_mRNA',
                                                       'cds_length'])

# from genome file (fasta), generate nucleotide sequence and chromosome label dictionaries
nuc_dict, chrom_dict = construct_chromosome_dictionaries(genome_file)

# from gtf file, generate gtf data frame, omitting rows with feature defined as 'gene' or 'start_codon'
gtf_df = pandas.read_csv(gtf_file, sep='\t', comment = '#', header = None)
gtf_df.columns = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
gtf_df = gtf_df[gtf_df['feature'] != 'start_codon']

# open phyloP file (BigWig)
phylop = pyBigWig.open(phylop_file)

# construct processed data frame and write to .tsv file
processed_df = construct_data_frame(gtf_df, nuc_dict, chrom_dict, phylop)
processed_df.to_csv(output_path, sep = '\t', index = False)

# close phyloP file (BigWig)
phylop.close()