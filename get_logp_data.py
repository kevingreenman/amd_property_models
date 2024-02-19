import json

import pandas as pd
import requests

from model_retraining import (
    create_new_model_directory,
    create_new_train_data_file,
    submit_new_chemprop_job,
)
from utils import (
    exclude_all_but_most_recent_date,
    get_date,
    get_mongodb_connection_info,
    run_mongodb_query,
)


def get_logp_from_mongodb(content):
    """Get the UV-Vis spectra using MongoDB IDs from a previous query."""

    http_address, mongo_port, username, password = get_mongodb_connection_info()

    # Extract data from special query
    ids = [
        content[1]["search_results"][x]["_id"]
        for x in range(0, len(content[1]["search_results"]))
    ]

    # Query by Object IDs found from special query above
    data_entries = {}
    for _id in ids:
        incoming_request = {
            "auth": {"port": mongo_port, "user": username, "password": password},
            "request": {
                "request_type": "query_data",
                "search_type": "object_id",
                "collection": "data",
                "object_id": _id,
            },
        }
        resp = requests.post(
            http_address,
            json=incoming_request,
            headers={"content-type": "application/authjson"},
        )
        id_content = json.loads(resp.content.decode("utf_8"))
        entry = id_content[1]["search_results"][0]["data_entry"]
        data_entries[_id] = entry

    df = pd.DataFrame()
    df["substance_name"] = [data_entries[_id]["data_name"].split("_")[0] for _id in ids]
    df["smiles"] = [data_entries[_id]["compound_smiles"] for _id in ids]
    df["solvent"] = [
        data_entries[_id]["solvent"] for _id in ids
    ]
    df["experimentdate"] = [data_entries[_id]["origin_date"] for _id in ids]
    df["origin1"] = [data_entries[_id]["data_origin"][0] for _id in ids]
    df["origin2"] = [data_entries[_id]["data_origin"][1] for _id in ids]
    # df['campaign_name'] = [data_entries[_id]['campaign_name'] for _id in ids]
    df["logP"] = [
        data_entries[_id]["logp_value"] for _id in ids
    ]
    df["mongodbid"] = ids

    df.dropna(subset=["smiles", "logP"], inplace=True)

    # Only return database entries that were measured on the AMD platform
    df = df[
        (df["origin1"] == "measured") & (df["origin2"] == "amd_platform_HPLC")
    ].drop(columns=["origin1", "origin2"])

    # Convert date from string to datetime.time object
    df["experimentdate"] = df["experimentdate"].apply(lambda x: get_date(x))

    return df


def main():
    content = run_mongodb_query(
        search_type="special_search", collection="data", search_term=["exp_logp_hplc"]
    )
    df = get_logp_from_mongodb(
        content
    )
    df = exclude_all_but_most_recent_date(
        df
    )
    create_new_train_data_file(
        df,
        old_results_dir_name=None,
        old_results_file_name=None,
        new_results_dir_name="v3",
        new_results_file_name="v3.csv",
        property="logp",
    )
    create_new_model_directory(
        new_model_dir_name="logp_model_v4",
        new_results_dir_name="v3",
        new_results_file_name="v3.csv",
        property="logp",
    )
    submit_new_chemprop_job(
        new_model_dir_name="logp_model_v4",
        property="logp",
    )


if __name__ == "__main__":
    main()
