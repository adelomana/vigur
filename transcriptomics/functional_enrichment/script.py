###
# This script reads json files from http://amigo.geneontology.org/amigo and format data into tables for publication
###

import os, json

import pyensembl
release = pyensembl.EnsemblRelease()

def formatter(enrichment, g):

    print(enrichment)
    print()

    term = enrichment['term']['label']

    level = 0
    if term != 'UNCLASSIFIED':
        level = enrichment['term']['level']

    background_rank = enrichment['number_in_reference']
    found_rank = enrichment['input_list']['number_in_list']
    expected_rank = enrichment['input_list']['expected']
    fold_enrichment = enrichment['input_list']['fold_enrichment']
    sign = enrichment['input_list']['plus_minus']
    pvalue = enrichment['input_list']['pValue']
    mapped_ids = enrichment['input_list']['mapped_id_list']['mapped_id']

    gene_names_string = ', '.join(mapped_ids)

    if level == 0:
        g.write('\n')

    g.write('{}\t'.format(level))
    g.write('{}\t'.format(term))
    g.write('{}\t'.format(background_rank))
    g.write('{}\t'.format(found_rank))
    g.write('{}\t'.format(expected_rank))
    g.write('{}\t'.format(fold_enrichment))
    g.write('{}\t'.format(sign))
    g.write('{}\t'.format(pvalue))
    g.write('{}\t'.format(gene_names_string))

    g.write('\n')

    return None

# 0. user-defined variables
json_dir = '/home/adrian/projects/vigur/results/transcriptomics/annotation/json_dir/'
table_dir = '/home/adrian/projects/vigur/results/transcriptomics/annotation/table_dir/'

# 1. read the data
json_files = os.listdir(json_dir)

# 2. create tables
for json_file in json_files:

    print('working with {}'.format(json_file))

    output_file = table_dir + json_file.replace('.json', '') + '.tsv'
    input_file = json_dir + json_file

    g = open(output_file,'w')
    g.write('Level\tTerm\tBackground rank\tFound rank\tExpected rank\tFold enrichment\tSign\tP value\tGene names\n')

    with open(input_file, 'r') as f:
        data = json.load(f)
        for category in data['overrepresentation']['group']:
            if category != '':
                if isinstance(category['result'], list):
                    for working_dictionary in category['result']:
                        formatter(working_dictionary, g)
                elif isinstance(category['result'], dict):
                    working_dictionary = category['result']
                    formatter(working_dictionary, g)
                else:
                    raise ValueError('unexpected element detected')
