import sys, numpy, pandas

input_file = 'EC_Alexia_EC72_LIPIDOMIC_TARGETED_POS_B2_1.txt'

#
# 1. read and format info
#
quantification ={}
with open(input_file, 'rb') as f:
    for line in f:

        #
        # fix format issues
        #
        decoded_line = line.decode('unicode_escape')
        v = decoded_line.split('\t')
        info = [element.replace('\r\n', '') for element in v]

        if 'Compound' in info[0]:
            compound = v[0].replace('\r\n', '')
            container = {}

        #
        # recover info
        #
        if len(info) >= 3:
            if 'Lipidomics' in info[2]:
                sample_name = info[2] + '_' + info[-1]
                intensity_string = info[-3]
                if len(intensity_string) != 0:
                    intensity = float(intensity_string)
                else:
                    intensity = numpy.nan
                container[sample_name] = intensity

        #
        # save
        #
        if len(info) == 1:
            if len(container.keys()) != 0:

                quantification[compound] = container


#
# 2. convert to df and store
#
df = pandas.DataFrame(quantification)
df.sort_index(axis='index', inplace=True)
df.to_csv('quantification.tsv', sep='\t')
df.to_excel('quantification.xlsx')
