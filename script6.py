# requires pandas
import pandas

# define file paths
uORFs_file = "/Users/santiagochang/Desktop/JSB/WBcel235/output_uORF.tsv"
ASequences_file = "/Users/santiagochang/Desktop/JSB/WBcel235/updated_ASequences.tsv"
output_file = "/Users/santiagochang/Desktop/JSB/WBcel235/ASequences_with_score.tsv"

# define nucleotide dictionary
nucleotide_dictionary = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

# input: uORFs data frame with each row representing a CDS
# output: matrix with proportions of nucleotides by each position -4 through -1
def proportions_by_position(uORFs_df):
    # define proportions matrix s.t. proportions[pos][n] corresponds to nucleotide x at position i-5
    proportions = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    # only consider CDS with minimum 4 preceding nucleotides
    filtered_df = uORFs_df[(uORFs_df['cds_start_in_mRNA'] >= 4)]
    total = 0
    # iterate by row of data frame
    for row in filtered_df.itertuples(index=False):
        # increment total
        total += 1
        # get 4 preceding nucleotides of CDS
        cdna = row[2].upper()
        cds_start = int(row[11])
        kozak = cdna[cds_start-4:cds_start]
        # increment matrix values based on position and nucleotide
        for pos in range(len(kozak)):
            proportions[pos][nucleotide_dictionary[kozak[pos]]] += 1
    # divide matrix elements by total to get proportion
    for pos in range(len(proportions)):
        for n in range(len(proportions[pos])):
            proportions[pos][n] /= total
    # output proportions
    print(f'Pos -4: A = {proportions[0][0]}; T = {proportions[0][1]}; G = {proportions[0][2]}; C = {proportions[0][3]}')
    print(f'Pos -3: A = {proportions[1][0]}; T = {proportions[1][1]}; G = {proportions[1][2]}; C = {proportions[1][3]}')
    print(f'Pos -2: A = {proportions[2][0]}; T = {proportions[2][1]}; G = {proportions[2][2]}; C = {proportions[2][3]}')
    print(f'Pos -1: A = {proportions[3][0]}; T = {proportions[3][1]}; G = {proportions[3][2]}; C = {proportions[3][3]}')
    # return proportions matrix
    print(total)
    return proportions

# input: ASequences data frame, proportions matrix
# output: data frame containing only uORFs with 4 preceding nucleotides, with kozak score and phyloP score (for ATG)
def calculate_kozak_score(ASequences_df, prop_mat):
    # calculate max score possible
    max_score = 0
    for i in range(len(prop_mat)):
        max_score += max(prop_mat[i])
    # filter out rows with uORFs that lack 4 preceding nucleotides (kozak sequence)
    df = ASequences_df[ASequences_df['kozak'].str.len() == 7]
    df = df.reset_index(drop=True)
    # iterate by row of data frame
    for row in range(len(df)):
        # get kozak sequence
        sequence = df.loc[row, 'kozak'][0:4]
        score = 0
        # calculate score
        for i in range(len(sequence)):
            score += prop_mat[i][nucleotide_dictionary[sequence[i]]]
        score /= max_score
        # append score to data frame
        df.loc[row, 'kozak_score'] = score
        # get mean phyloP score of ATG nucleotides
        phylop = df.loc[row, 'phyloP_scores'].split(',')[0:3]
        phylop = [float(s) for s in phylop if s != 'NA']
        # append score to data frame
        df.loc[row, 'phyloP_ATG'] = sum(phylop)/3
        # get mean phyloP score of end codon nucleotides
        phylop = df.loc[row, 'phyloP_scores'].split(',')[-3:]
        phylop = [float(s) for s in phylop if s != 'NA']
        # append score to data frame
        df.loc[row, 'phyloP_stop'] = sum(phylop)/3
    return df

# read in files
uORFs = pandas.read_csv(uORFs_file, sep='\t')
ASequences = pandas.read_csv(ASequences_file, sep='\t')
# get proportions
prop = proportions_by_position(uORFs)
# calculate kozak score and ATG mean phyloP score
updated_df = calculate_kozak_score(ASequences, prop)
# output to file
updated_df.to_csv(output_file, sep = '\t', index = False)