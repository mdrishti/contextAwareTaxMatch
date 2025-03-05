import pandas as pd
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('wd_sparql_file', type=str, help="Enter the file which contains wikidata mappings to other databases")
parser.add_argument('verbatim_file', type=str, help="Enter the file which contains verbatim interactions for GloBI")
args = parser.parse_args()
wd_sparql_file = args.wd_sparql_file
verbatim_file = args.verbatim_file



# 1-read and transform the CSV file
wd_sparql_df = pd.read_csv(wd_sparql_file, sep=",", dtype=str)

# 1-rename columns (Skipping the first row)
prefixes = {
    1: "EOL:", 2: "GBIF:", 3: "NCBI:", 4: "OTT:", 5: "ITIS:",
    6: "IRMNG:", 7: "COL:", 8: "NBN:", 9: "WORMS:", 10: "BOLD:",
    11: "PLAZI:", 12: "APNI:",  13: "msw3:", 14: "INAT_TAXON:", 15: "EPPO:"
}

for col, prefix in prefixes.items():
    if col < len(wd_sparql_df.columns):
        wd_sparql_df.iloc[0:, col] = prefix + wd_sparql_df.iloc[0:, col]

# 1-replace URL patterns
wd_sparql_df.replace({"http://www.wikidata.org/entity/": "Wikidata:", '"': ''}, regex=True, inplace=True)
#wd_sparql_df.to_csv(output_file, sep="\t", index=False, header=False)

# 2- process the gzipped verbatim interactions file
verbatim_globi_df = pd.read_csv(verbatim_file, usecols=['sourceTaxonId','sourceTaxonName'], sep="\t", dtype=str) #read source
verbatim_globi_df.columns = ["TaxonId", "TaxonName"]
df1 = pd.read_csv(verbatim_file, usecols=['targetTaxonId','targetTaxonName'], sep="\t", dtype=str) #read target
df1.columns = ["TaxonId", "TaxonName"]

verbatim_globi_df = pd.concat([verbatim_globi_df, df1], axis=0) #concat source and target

# 2- replace placeholders and sort -u
verbatim_globi_df.replace({
    "https://www.wikidata.org/wiki/": "Wikidata:",
    "https://www.wikidata.org/entity/": "Wikidata:",
    "urn:lsid:marinespecies.org:taxname": "WORMS",
    "urn:lsid:irmng.org:taxname": "IRMNG",
    "http://www.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=BOLD": "BOLD",
    "https://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=": "ITIS:",
    "https://www.inaturalist.org/taxa/": "INAT_TAXON:",
    "https://www.gbif.org/species/": "GBIF:",
    "https://species.nbnatlas.org/species/": "NBN:",
    "https://gd.eppo.int/taxon/": "EPPO:",
    "^tsn": "ITIS",
    "GBIF: +": "GBIF:",
    "gbif: +": "GBIF:",
    "gbif:": "GBIF:",
}, regex=True, inplace=True)
verbatim_globi_df = verbatim_globi_df.drop_duplicates()
#verbatim_df.to_csv(out_verbatim_file, sep="\t", index=False, header=False)


# 3-merge wd_sparql_temp1 with verbatim interactions
id_map = {row[i]: row[len(row) - 1] for _, row in wd_sparql_df.iterrows() for i in range(len(row) - 1) if pd.notna(row[i])}


final_output = "wd_globi_source_target_names.txt"
with open(final_output, 'w') as out:
    for _, row in verbatim_globi_df.iterrows():
        key, value = row[0], row[1]
        if key != "":
            if key in id_map.keys():
                mapped_id, mapped_value = key, id_map[key]
                print(mapped_id," ",mapped_value)
                match_status = "NAME-MATCH-YES" if mapped_value.lower() == value.lower() else "NAME-MATCH-NO"
                out.write(f"{key}\t{value}\t{mapped_id}\t{mapped_value}\t{match_status}\n")
            else:
                out.write(f"{key}\t{value}\t\tID-NOT-FOUND\n")
        else:
            out.write(f"{key}\t{value}\t\tID-NOT-PRESENT\n")

