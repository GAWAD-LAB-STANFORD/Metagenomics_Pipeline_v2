import pandas as pd
import argparse


parser = argparse.ArgumentParser(
    description="""Merges a BLAST TSV with a Kraken report TSV""")
parser.add_argument('-b', '--blast', help="input BLAST file", required=True)
parser.add_argument('-k', '--kraken', help="input Kraken file", required=True)
args = parser.parse_args()


blast_results_df = pd.read_csv(args.blast, sep="\t")
sample = blast_results_df['sample'].unique().tolist()[0]
columns = ['kraken_percent_fragments_covered', 'kraken_fragments_covered', 'kraken_fragments_assigned', 'rank_code', 
           'kraken_species_id', 'kraken_species', 'kraken_domain', 'kraken_kingdom', 'kraken_phylum', 'kraken_class', 
           'kraken_order', 'kraken_family', 'kraken_genus']


# Determine taxonomic information for each species
kraken_df = pd.read_csv(args.kraken, sep="\t")
kraken_df.columns = ['percent_fragments_covered', 'fragments_covered', 'fragments_assigned', 'rank_code', 'taxid', 'sciname']
kraken_species_list = []
domain, kingdom, phylum, k_class, order, family, genus = '', '', '', '', '', '', ''
for index, row in kraken_df.iterrows():
    if row['rank_code'] == 'U' or row['rank_code'] == 'R':
        domain, kingdom, phylum, k_class, order, family, genus = '', '', '', '', '', '', ''
    elif row['rank_code'] == 'D':
        domain, kingdom, phylum, k_class, order, family, genus = row['sciname'], '', '', '', '', '', ''
    elif row['rank_code'] == 'K':
        kingdom, phylum, k_class, order, family, genus = row['sciname'], '', '', '', '', ''
    elif row['rank_code'] == 'P':
        phylum, k_class, order, family, genus = row['sciname'], '', '', '', ''
    elif row['rank_code'] == 'C':
        k_class, order, family, genus = row['sciname'], '', '', ''
    elif row['rank_code'] == 'O':
        order, family, genus = row['sciname'], '', ''
    elif row['rank_code'] == 'F':
        family, genus = row['sciname'], ''
    elif row['rank_code'] == 'G':
        genus = row['sciname']
    elif row['rank_code'] == 'S' or row['rank_code'] == 'S1':
        new_species = row.tolist() + [domain.strip(), kingdom.strip(), phylum.strip(), k_class.strip(), order.strip(), family.strip(), genus.strip()]
        new_species[6] = new_species[6].strip()
        kraken_species_list.append(new_species)


kraken_species_df = pd.DataFrame(kraken_species_list, columns=columns)
kraken_species_df = kraken_species_df.assign(sample=sample)


# Ensure columns that will be compared are of the same data type
kraken_species_df['sample'] = kraken_species_df['sample'].astype(str)
blast_results_df['sample'] = blast_results_df['sample'].astype(str)
kraken_species_df['kraken_species_id'] = kraken_species_df['kraken_species_id'].astype(int)
blast_results_df['kraken_species_id'] = blast_results_df['kraken_species_id'].astype(int)


# Merge all data into single pandas dataframe
merged_df = pd.merge(blast_results_df, kraken_species_df, how='left', on=['sample', 'kraken_species_id'])
merged_df.to_csv(args.blast, header = True, index = False, sep="\t")