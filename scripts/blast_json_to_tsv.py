import pandas as pd
import argparse
import json

parser = argparse.ArgumentParser(
    description="""Converts a BLAST JSON file into a TSV""")
parser.add_argument('-i', '--input', help="input BLAST JSON", required=True)
parser.add_argument('-o', '--output', help="output BLAST TSV", required=True)
parser.add_argument('-b', '--blast_db', help="BLAST DB type", default="nt")
parser.add_argument('-s', '--sample', help="sample", default="")
parser.add_argument('-k', '--kraken_db', help="Kraken DB type", default="")
parser.add_argument('-d', '--kraken_species_id', help="Kraken species ID", default="")
args = parser.parse_args()


blast_hits_list = []
json_file = open(args.input, "r")
blastn_dict = json.load(json_file)
json_file.close()
for blast_result in blastn_dict['BlastOutput2']:
    query_name = blast_result['report']['results']['search']['query_title']
    query_length = blast_result['report']['results']['search']['query_len']
    hits = blast_result['report']['results']['search']['hits']
    if len(hits) > 0:
        for hit_count in range(len(hits)):
            hit =  hits[hit_count]
            if args.blast_db == "plasmid":
                new_addition = [hit_count + 1, query_name, query_length, 
                                hit['description'][0]['title'], hit['description'][0]['accession'], 
                                hit['hsps'][0]['align_len'], hit['len'], (hit['hsps'][0]['align_len']/query_length)*100,
                                hit['hsps'][0]['hit_from'], hit['hsps'][0]['hit_to']]
            elif args.blast_db == "viral":
                new_addition = [hit_count + 1, query_name, query_length,
                                hit['description'][0]['title'], hit['description'][0]['accession'], 
                                hit['hsps'][0]['align_len'], hit['len'], (hit['hsps'][0]['align_len']/query_length)*100,
                                hit['hsps'][0]['hit_from'], hit['hsps'][0]['hit_to']]
            else:
                new_addition = [hit_count + 1, query_name, query_length,
                                hit['description'][0]['sciname'], hit['description'][0]['sciname'].split(" ")[0],
                                hit['description'][0]['taxid'], hit['description'][0]['accession'], 
                                hit['hsps'][0]['align_len'], hit['len'], (hit['hsps'][0]['align_len']/query_length)*100,
                                hit['hsps'][0]['hit_from'], hit['hsps'][0]['hit_to']]
            if '' not in new_addition:
                blast_hits_list.append(new_addition)


if len(blast_hits_list) > 0:
    print("{} blast hits found".format(len(blast_hits_list)))
    if args.blast_db == "plasmid":
        blast_results_df = pd.DataFrame(blast_hits_list, 
                                    columns=['hit_rank', 'query', 'query_length', 'blast_plasmid', 'accession', 
                                    'top_alignment_length', 'reference_length', 'percent_aligned', 'hit_from', 'hit_to'])
    elif args.blast_db == "viral":
        blast_results_df = pd.DataFrame(blast_hits_list, 
                                    columns=['hit_rank', 'query', 'query_length', 'blast_virus', 'accession', 
                                    'top_alignment_length', 'reference_length', 'percent_aligned', 'hit_from', 'hit_to'])
    else:
        blast_results_df = pd.DataFrame(blast_hits_list, 
                                        columns=['hit_rank', 'query', 'query_length', 'blast_species', 'blast_genus',
                                        'hit_taxid', 'accession', 'top_alignment_length', 'reference_length', 'percent_aligned',
                                        'hit_from', 'hit_to'])
    blast_results_df = blast_results_df.assign(blast_db=args.blast_db)
    if len(args.sample) > 0:
        blast_results_df = blast_results_df.assign(sample=args.sample)
    if len(args.kraken_db) > 0:
        blast_results_df = blast_results_df.assign(kraken_db=args.kraken_db)
    if len(args.kraken_species_id) > 0:
        blast_results_df = blast_results_df.assign(kraken_species_id=args.kraken_species_id)
    blast_results_df.to_csv(args.output, header = True, index = False, sep="\t")
else:
    print("No blast hits found in BLAST JSON")