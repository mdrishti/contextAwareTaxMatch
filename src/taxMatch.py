import pandas as pd
import gzip
import sys
import argparse
from itertools import zip_longest
import time 

sys.path.append('./functions')  # Add the 'src' directory to the sys.path
import dataProcessing as dpx

predefined_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# function to extract ranks into separate columns
def extract_ranks(structure, values):
    # split structure and values into lists
    rank_list = [r.strip() for r in structure.split("|")]
    value_list = [v.strip() for v in values.split("|")]
    rank_dict = dict(zip_longest(rank_list, value_list, fillvalue=""))
    return {rank: rank_dict.get(rank, "") for rank in predefined_ranks}

# function to split taxon paths and rank names
def safe_extract_ranks(row):
    if pd.notna(row["TaxonPathName"]) and pd.notna(row["TaxonRankName"]):
        return extract_ranks(row["TaxonRankName"], row["TaxonPathName"])
    else:
        return pd.Series({rank: pd.NA for rank in predefined_ranks})

parser = argparse.ArgumentParser()
parser.add_argument('wd_sparql_file', type=str, help="Enter the file which contains wikidata mappings to other databases")
parser.add_argument('verbatim_file', type=str, help="Enter the file which contains verbatim interactions for GloBI")
parser.add_argument('wd_lineage_aligned_file', type=str, help="Enter the wd lineage file")
parser.add_argument('output_file', type=str, help="Enter the file to store the taxonomy matches")
parser.add_argument('wdLineageRepeats_file', type=str, help="Enter the file having taxonomic names which are repeated in the wd lineage file")
args = parser.parse_args()
wd_sparql_file = args.wd_sparql_file
verbatim_file = args.verbatim_file
outputFileX = args.output_file
wd_lineage_file = args.wd_lineage_aligned_file
wd_repeats_file = args.wdLineageRepeats_file


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

# 1-map of ids
cols_to_map = wd_sparql_df.columns[:-1]  # all columns except the last one
id_map = (
    wd_sparql_df.melt(id_vars=wd_sparql_df.columns[-1], value_vars=cols_to_map, value_name="key")
    .dropna(subset=["key"])
    .set_index("key")[wd_sparql_df.columns[-1]]
    .to_dict()
)
cols_to_map_WD = wd_sparql_df.columns[1:-1]  # all columns except the first and last one
id_map_WD = (
    wd_sparql_df.melt(id_vars=wd_sparql_df.columns[0], value_vars=cols_to_map, value_name="key")
    .dropna(subset=["key"])
    .set_index("key")[wd_sparql_df.columns[0]]
    .to_dict()
)


# 2- process the gzipped verbatim interactions file
verbatim_globi_df = pd.read_csv(verbatim_file, usecols=['sourceTaxonId','sourceTaxonName','sourceTaxonPathNames','sourceTaxonPathRankNames'], sep="\t", dtype=str) #read source
verbatim_globi_df.columns = ["TaxonId", "TaxonName","TaxonPathName","TaxonRankName"]
df1 = pd.read_csv(verbatim_file, usecols=['targetTaxonId','targetTaxonName','targetTaxonPathNames','targetTaxonPathRankNames'], sep="\t", dtype=str) #read target
df1.columns = ["TaxonId", "TaxonName","TaxonPathName","TaxonRankName"]
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
verbatim_globi_df_backup= verbatim_globi_df.copy()

# 2 - add ranks to the df
expanded_verbatim_globi_df = verbatim_globi_df.apply(safe_extract_ranks, axis=1, result_type="expand")
verbatim_globi_df = pd.concat([verbatim_globi_df, expanded_verbatim_globi_df], axis=1)
verbatim_globi_df_backup1= verbatim_globi_df.copy()
#verbatrim_globi_df = pd.concat([verbatrim_globi_df, expanded_verbatim_globi_df], axis=1)
#verbatim_df.to_csv(out_verbatim_file, sep="\t", index=False, header=False)

# 3- first layer of extractions using wd_sparql_df and verbatim interactions. Assign NAME-MATCH-YES/NO, ID-NOT-FOUND, ID-NOT-PRESENT
verbatim_globi_df = dpx.initialTaxMatchDfY(verbatim_globi_df, id_map, id_map_WD)
verbatim_globi_df_backup2= verbatim_globi_df.copy()


# 4- second layer of extraction #########################
# 4- process wd lineage file
wd_lineage_df = pd.read_csv(wd_lineage_file, sep=",", dtype=str)
wd_lineage_df["WdID"] = wd_lineage_df["WdID"].str.replace("http://www.wikidata.org/entity/", "Wikidata:", regex=False)
wd_name_to_id_set = set(wd_lineage_df["WdName"])

# 4- process repeated taxonomic names in the lineage
wd_repeats_lineage = pd.read_csv(wd_repeats_file)
wd_repeats_lineage = wd_repeats_lineage.fillna("")
wd_repeats_lineage["WdID"] = wd_repeats_lineage["WdID"].str.replace("http://www.wikidata.org/entity/", "Wikidata:", regex=False)
wd_lineage_set = set(wd_repeats_lineage["WdName"])
matched_df = pd.DataFrame(columns=['TaxonId', 'TaxonName', 'TaxonPathName', 'TaxonRankName', 'kingdom','phylum', 'class', 'order', 'family', 'genus', 'species', 'Mapped_ID','Mapped_Value', 'Match_Status','Status'])

#4- make a set and dict each of the lineage and the repeats BEFORE the next two functions of get_best_wikidata_id and process_row
wd_lineage_dict = (
    wd_repeats_lineage.set_index(["WdName", "family", "class", "order", "phylum", "kingdom"])["WdID"]
    .groupby(level=[0, 1, 2, 3, 4, 5])
    .apply(list)
    .to_dict()
)
mask = ~wd_lineage_df["WdName"].isin(wd_lineage_set) #~ depicts not in the wd_lineage_set
wd_name_to_id = wd_lineage_df.loc[mask].set_index("WdName")[
    ["WdID", "family", "class", "order", "phylum", "kingdom"]
].apply(lambda x: tuple(x), axis=1).to_dict()


#4-function to score and get the best wikidata-id according to rank matching in the case of wd repeats
def get_best_wikidata_id(taxon_name, family, tax_class, order, phylum, kingdom):
    possible_keys = [k for k in wd_lineage_dict.keys() if k[0] == taxon_name and pd.notna(k[0]) and k[0] != ""]
    best_match = None
    best_score = 0
    for key in possible_keys:
        score=0
        ranks = [family, tax_class, order, phylum, kingdom]
        score = sum([
            1 if (pd.notna(key[1]) and key[1] != "" and key[1] == family) else 0,
            1 if (pd.notna(key[2]) and key[2] != "" and key[2] == tax_class) else 0,
            1 if (pd.notna(key[3]) and key[3] != "" and key[3] == order) else 0,
            1 if (pd.notna(key[4]) and key[4] != "" and key[4] == phylum) else 0,
            1 if (pd.notna(key[5]) and key[5] != "" and key[5] == kingdom) else 0
        ])
        if score > best_score:
            best_match = key
            best_score = score
    return (best_match, wd_lineage_dict[best_match]) if best_match else (None,None)


#4-function to process the rows for "ID-NOT-FOUND" cases by checking through three cases -repeats elif direct match elif still not found
def process_row(row):
    taxon_name = row["TaxonName"]
    #if (taxon_name,) in wd_lineage_dict:
    possible_keys = [k for k in wd_lineage_dict.keys() if k[0] == taxon_name]
    if possible_keys:
        #best_wd_ids = wd_name_to_id[possible_keys[0]]  # Or apply ranking logic
        family = row.get("family", "")
        family = family if pd.notna(family) else ""
        tax_class = row.get("class", "")
        tax_class = tax_class if pd.notna(tax_class) else ""
        order = row.get("order", "")
        order = order if pd.notna(order) else ""
        phylum = row.get("phylum", "")
        phylum = phylum if pd.notna(phylum) else ""
        kingdom = row.get("kingdom", "")
        kingdom = kingdom if pd.notna(kingdom) else ""
        tempVar, tempVarX = get_best_wikidata_id(taxon_name, family, tax_class, order, phylum, kingdom) # return both key and value - key (WdName, ranks...) is first, value (WdID) is later
        if tempVar:
            best_wd_id = tempVarX
            row["family"] = tempVar[1]
            row["class"] = tempVar[2]
            row["order"] = tempVar[3]
            row["phylum"] = tempVar[4]
            row["kingdom"] = tempVar[5]
        else:
            best_wd_id = None
            #best_wd_id = row["TaxonId"]
        status = "ID-MATCHED-BY-NAME-DUPL-duplicate" if tempVar else "ID-MATCHED-BY-NAME-DUPL-mismatch"
    elif taxon_name in wd_name_to_id_set:
        tempVar = wd_name_to_id.get(taxon_name, (None,)) #retrieve full row-Value followed by index-based alignment
        best_wd_id = tempVar[0]
        row["family"] = tempVar[1]
        row["class"] = tempVar[2]
        row["order"] = tempVar[3]
        row["phylum"] = tempVar[4]
        row["kingdom"] = tempVar[5]
        #best_wd_id = wd_name_to_id[(taxon_name,)]
        status = "ID-MATCHED-BY-NAME-direct"
    else:
        best_wd_id = None
        status = row["Match_Status"]
        #best_wd_id = row["TaxonId"]
    row["Mapped_ID_WD"] = best_wd_id
    row["Match_Status"] = status
    return row

#4- Match as many ID-NOT-FOUND as possible
start=time.time() #sanity check for time
#matched_df_part1 = verbatim_globi_df.loc[verbatim_globi_df["Match_Status"] == "ID-NOT-FOUND"].iloc[1:100].apply(process_row, axis=1) # snaity check for few rows
#matched_df = verbatim_globi_df.loc[verbatim_globi_df["Match_Status"] == "ID-NOT-FOUND"].apply(process_row, axis=1) # separate df 
mask = verbatim_globi_df["Match_Status"] == "ID-NOT-FOUND" 
verbatim_globi_df.loc[mask] = verbatim_globi_df.loc[mask].apply(process_row, axis=1).apply(pd.Series) #same df
verbatim_globi_df_backup3= verbatim_globi_df.copy()
mask = verbatim_globi_df["Match_Status"] == "ID-NOT-PRESENT" 
verbatim_globi_df.loc[mask] = verbatim_globi_df.loc[mask].apply(process_row, axis=1).apply(pd.Series) #same df
verbatim_globi_df_backup4= verbatim_globi_df.copy()
end=time.time() - start
print(end)
verbatim_globi_df.to_csv(outputFileX, index=False)



# code that could be useful 
'''wd_name_to_id = (
    wd_lineage_df.set_index(["WdName", "family", "class", "order", "phylum", "kingdom"])["WdID"]
    .groupby(level=[0, 1, 2, 3, 4, 5])
    .apply(list)
    .to_dict()
)
wd_name_to_id = wd_lineage_df[wd_lineage_df["WdName"] not in wd_lineage_set].set_index("WdName")["WdID"].to_dict()
wd_name_to_id = (
    wd_lineage_df[~wd_lineage_df["WdName"].isin(wd_lineage_set)] #~ depicts not inthe wd_lineage_set
    .set_index("WdName")["WdID"]
    .to_dict()
)'''



