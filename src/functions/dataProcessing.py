import pandas as pd
import gzip

def initialTaxMatchDf(verbatim_globi_df, id_map):
    records = []
    for _, row in verbatim_globi_df.iterrows():
        key, value = row.iloc[0], row.iloc[1]
        #key, value = row[0], row[1]
        if pd.notna(key):
            if key in id_map:
                mapped_id, mapped_value = key, id_map[key]
                match_status = "NAME-MATCH-YES" if mapped_value.lower() == value.lower() else "NAME-MATCH-NO"
                records.append({"Key": key, "Value": value, "Mapped_ID": mapped_id, 
                                "Mapped_Value": mapped_value, "Match_Status": match_status})
            else:
                records.append({"Key": key, "Value": value, "Mapped_ID": None, 
                                "Mapped_Value": None, "Match_Status": "ID-NOT-FOUND"})
        else:
            records.append({"Key": key, "Value": value, "Mapped_ID": None, 
                            "Mapped_Value": None, "Match_Status": "ID-NOT-PRESENT"})
    return pd.DataFrame(records)  # convert list to DataFrame


def initialTaxMatchDfX(verbatim_globi_df, id_map):
    # Initialize columns
    verbatim_globi_df["Mapped_ID"] = None
    verbatim_globi_df["Mapped_Value"] = None
    verbatim_globi_df["Match_Status"] = None
    for idx, row in verbatim_globi_df.iterrows():
        key, value = row.iloc[0], row.iloc[1]
        if pd.notna(key):
            if key in id_map:
                mapped_id, mapped_value = key, id_map[key]
                match_status = "NAME-MATCH-YES" if mapped_value.lower() == value.lower() else "NAME-MATCH-NO"
            else:
                mapped_id, mapped_value, match_status = None, None, "ID-NOT-FOUND"
        else:
            mapped_id, mapped_value, match_status = None, None, "ID-NOT-PRESENT"
        verbatim_globi_df.at[idx, "Mapped_ID"] = mapped_id
        verbatim_globi_df.at[idx, "Mapped_Value"] = mapped_value
        verbatim_globi_df.at[idx, "Match_Status"] = match_status
    return verbatim_globi_df


def initialTaxMatchDfY(verbatim_globi_df, id_map):
    # map TaxonID to corresponding TaxonName from id_map (returns NaN if not found)
    verbatim_globi_df["Mapped_Value"] = verbatim_globi_df.iloc[:, 0].map(id_map)
    # Assign Mapped_ID (same as TaxonID if found in id_map, else None)
    verbatim_globi_df["Mapped_ID"] = verbatim_globi_df.iloc[:, 0].where(
        verbatim_globi_df["Mapped_Value"].notna()
    )
    # compare TaxonName with mapped value (case insensitive)
    verbatim_globi_df["Match_Status"] = (
        verbatim_globi_df["Mapped_Value"].str.lower()
        == verbatim_globi_df.iloc[:, 1].str.lower()
    ).map({True: "NAME-MATCH-YES", False: "NAME-MATCH-NO"})
    # handle cases where TaxonID is missing in id_map
    verbatim_globi_df["Match_Status"] = verbatim_globi_df["Match_Status"].where(
        verbatim_globi_df["Mapped_Value"].notna(), "ID-NOT-FOUND"
    )
    # Handle cases where TaxonID is NA or empty
    verbatim_globi_df["Match_Status"] = verbatim_globi_df["Match_Status"].where(
        verbatim_globi_df.iloc[:, 0].notna() & (verbatim_globi_df.iloc[:, 0] != ""),
        "ID-NOT-PRESENT",
    )
    return verbatim_globi_df

