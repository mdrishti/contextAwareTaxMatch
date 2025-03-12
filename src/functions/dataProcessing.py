import pandas as pd
import gzip


def initialTaxMatch(final_output, verbatim_globi_df, id_map):
    with open(final_output, 'w') as out:
        for _, row in verbatim_globi_df.iterrows():
            key, value = row[0], row[1]
            if pd.notna(key):
                if key in id_map.keys():
                    mapped_id, mapped_value = key, id_map[key]
                    print(mapped_id," ",mapped_value)
                    match_status = "NAME-MATCH-YES" if mapped_value.lower() == value.lower() else "NAME-MATCH-NO"
                    out.write(f"{key}\t{value}\t{mapped_id}\t{mapped_value}\t{match_status}\n")
                else:
                    out.write(f"{key}\t{value}\t\tID-NOT-FOUND\n")
            else:
                out.write(f"{key}\t{value}\t\tID-NOT-PRESENT\n")
    return

def initialTaxMatchDf(verbatim_globi_df, id_map):
    records = []
    for _, row in verbatim_globi_df.iterrows():
        key, value = row[0], row[1]
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
        key, value = row[0], row[1]
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
