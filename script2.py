# requires pandas
import pandas

# define file paths for input data frame and output data frame
input_file = "/Users/santiagochang/Desktop/JSB/WBcel235/filtered_output.tsv"
output_file = "/Users/santiagochang/Desktop/JSB/WBcel235/phylop_by_position.tsv"

# input: filtered data frame generated by script.py
# output: data frame with phyloP score, position, and strand
def organize_scores_by_position(data_frame):
    # initialize array to hold rows and int to track progress
    data = []
    progress = 0
    # iterate through each row in data frame
    for index, row in data_frame.iterrows():
        # print progress in 1% intervals
        if int(index / (len(data_frame) - 1) * 100) == progress:
            print(f'Construction of data frame is {progress}% complete.')
            progress += 1
        # extract phyloP scores
        scores = row.iloc[9].split(',')
        # append row to data array with phyloP score, position, and strand
        for i in range(int(row.iloc[11]), int(row.iloc[11])+int(row.iloc[12]), 3):
            if scores[i] != 'NA':
                data.append([float(scores[i]), 1, row.iloc[8]])
            if scores[i+1] != 'NA':
                data.append([float(scores[i+1]), 2, row.iloc[8]])
            if scores[i+2] != 'NA':
                data.append([float(scores[i+2]), 3, row.iloc[8]])
    # convert to data frame
    phylop_scores_by_position = pandas.DataFrame(data, columns = ['score', 'position', 'strand'])
    # return data frame
    return phylop_scores_by_position

# read in input file
input_df = pandas.read_csv(input_file, sep='\t', comment = '#', header = 0)
# construct phyloP score data frame
output_df = organize_scores_by_position(input_df)
# output to file
output_df.to_csv(output_file, sep = '\t', index = False)