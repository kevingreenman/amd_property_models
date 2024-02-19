import json

import numpy as np
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


def get_peak_from_peaks_list(peaks_list):
    """Get the peak from a list of peaks.

    Args:
        peaks_list (list): List of peaks (with one element).

    Returns:
        float: Peak wavelength.
    """
    if len(peaks_list) == 0:
        peak = np.nan
    else:
        peak = peaks_list[0]
    return peak


def get_solvent(solvent_list):
    """Get a single solvent from a list of solvents (handle inconsistencies in solvent format in MongoDB).

    Args:
        solvent_list (list): List of solvents.

    Returns:
        str: Solvent.
    """
    if isinstance(solvent_list[-1], list):
        return solvent_list[-1][2]
    elif isinstance(solvent_list[-1], str):
        return solvent_list[-1]
    else:
        raise ValueError("solvent_list must be a list of strings or lists")


def get_uvvis_peaks_from_spectra(row):
    """Get the peaks from the UV-Vis spectra.

    Args:
        row (pd.Series): Row from the DataFrame.

    Returns:
        list: List of peaks.
    """
    peakwavs_max = row["wavelengths"][np.argmax(row["intensities"])]
    return peakwavs_max


def get_uvvis_spectra_from_mongodb(content, peaks_method="direct"):
    """Get the UV-Vis spectra using MongoDB IDs from a previous query.

    Args:
        content: Content from the MongoDB query.
        peaks_method (str, optional): Method for getting the peaks. Use "direct" if ["exp_abs_stick"] was used in the MongoDB query and "spectra" if ["exp_abs_spec"] was used.

    Returns:
        pd.DataFrame: DataFrame with the UV-Vis spectra.
    """

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
    df["substance_name"] = [
        data_entries[_id]["data_name"].split("_")[0] for _id in ids
    ]
    df["smiles"] = [data_entries[_id]["compound_smiles"] for _id in ids]
    df["solvent"] = [
        get_solvent(data_entries[_id]["solvent"]) for _id in ids
    ]
    df["wavelengths"] = [data_entries[_id]["xdata"] for _id in ids]
    df["intensities"] = [
        data_entries[_id]["ydata"][0] for _id in ids
    ]  # Using [0] drops uncertainties
    df["experimentdate"] = [data_entries[_id]["origin_date"] for _id in ids]
    df["origin1"] = [data_entries[_id]["data_origin"][0] for _id in ids]
    df["origin2"] = [data_entries[_id]["data_origin"][1] for _id in ids]
    # df['campaign_name'] = [data_entries[_id]['campaign_name'] for _id in ids]
    if peaks_method == "direct":
        df["peakwavs_max"] = df["wavelengths"].apply(get_peak_from_peaks_list)
    elif peaks_method == "spectra":
        df["peakwavs_max"] = df.apply(
            get_uvvis_peaks_from_spectra, axis=1
        )
    else:
        raise ValueError("peaks_method must be 'direct' or 'spectra'")
    df["mongodbid"] = ids

    df.dropna(subset=["smiles", "solvent", "peakwavs_max"], inplace=True)

    # Only return database entries that were measured on the AMD platform
    df = df[(df["origin1"] == "measured") & (df["origin2"] == "amd_platform")].drop(
        columns=["origin1", "origin2"]
    )

    # Convert date from string to datetime.time object
    df["experimentdate"] = df["experimentdate"].apply(lambda x: get_date(x))

    return df


def main():
    content = run_mongodb_query(
        search_type="special_search", collection="data", search_term=["exp_abs_stick"]
    )
    df = get_uvvis_spectra_from_mongodb(
        content,
        peaks_method="direct",
    )
    df = exclude_all_but_most_recent_date(
        df
    )
    create_new_train_data_file(
        df,
        old_results_dir_name="v2c",
        old_results_file_name="smiles_target_train_v1_v2abc.csv",
        new_results_dir_name="v3",
        new_results_file_name="smiles_target_train_v1_v2abc_v3.csv",
        property="uvvis",
    )
    create_new_model_directory(
        new_model_dir_name="20221020_new_uvvis_measurements_close_loop_combined_v3",
        new_results_dir_name="v3",
        new_results_file_name="smiles_target_train_v1_v2abc_v3.csv",
        property="uvvis",
    )
    submit_new_chemprop_job(
        new_model_dir_name="20221020_new_uvvis_measurements_close_loop_combined_v3",
        property="uvvis",
    )


if __name__ == "__main__":
    main()