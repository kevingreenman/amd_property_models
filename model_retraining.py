import os
from pathlib import Path
from string import Template
import subprocess

import pandas as pd


def create_new_train_data_file(
    df,
    old_results_dir_name,
    old_results_file_name,
    new_results_dir_name,
    new_results_file_name,
    property="uvvis",
):
    """Create a new training data file."""

    datasets_dir = "datasets"

    new_data_dir = Path(f"{datasets_dir}/{property}/{new_results_dir_name}")
    new_data_dir.mkdir(parents=True, exist_ok=True)

    if property == "uvvis":
        cols = ["smiles", "solvent", "peakwavs_max"]
    elif property == "logp":
        cols = ["smiles", "logP"]
    elif property == "photodeg":
        cols = ["smiles", "photodeg_kinetic_rate"]
    else:
        raise ValueError("Invalid property.")

    new_df = df[cols].copy()

    if property == "uvvis":
        old_df = pd.read_csv(
            f"{datasets_dir}/{property}/{old_results_dir_name}/{old_results_file_name}"
        )

        combined_df = pd.concat([old_df, new_df], axis=0)

        new_results_dir = Path(f"{datasets_dir}/{property}/{new_results_dir_name}")
        new_results_dir.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(
            f"{datasets_dir}/{property}/{new_results_dir_name}/{new_results_file_name}",
            index=False,
        )
    elif property == "logp":
        new_df.to_csv(f"{datasets_dir}/{property}/{new_results_dir_name}/{new_results_file_name}", index=False)
    elif property == "photodeg":
        new_df.to_csv(f"{datasets_dir}/{property}/{new_results_dir_name}/{new_results_file_name}", index=False)

    return


def create_new_model_directory(
    new_model_dir_name: str,
    new_results_dir_name: str,
    new_results_file_name: str,
    property="uvvis",
):
    models_dir = f"models/{property}/production"
    new_model_dir = Path(f"{models_dir}/{new_model_dir_name}")
    new_model_dir.mkdir(parents=True, exist_ok=True)

    if property == "uvvis":
        replace_dict = {
            "TRAIN_DATA_FILE_FROM_PYTHON": f"{property}/{new_results_dir_name}/{new_results_file_name}"
        }
    elif property == "logp":
        replace_dict = {
            "TRAIN_DATA_FILES_LIST_FROM_PYTHON": "v0/smiles_target_all.csv v1/v1.csv" + f" {new_results_dir_name}/{new_results_file_name}"
        }
    elif property == "photodeg":
        replace_dict = {
            "TRAIN_DATA_FILE_FROM_PYTHON": f"{property}/{new_results_dir_name}/{new_results_file_name}",
        }

    with open(f"models/{property}/production/TEMPLATE_run_chemprop.sh", "r") as f:
        src = Template(f.read())
        result = src.substitute(replace_dict)

    with open(f"{models_dir}/{new_model_dir_name}/run_chemprop.sh", "w") as f:
        f.write(result)

    return


def submit_new_chemprop_job(
    new_model_dir_name: str,
    output_file_name: str = "slurm_output.txt",
    property="uvvis",
):
    os.chdir(f"models/{property}/production/{new_model_dir_name}")
    with open(output_file_name, 'w') as f:
        proc = subprocess.Popen(
            ["sbatch", "run_chemprop.sh"], stdout=f, stderr=f
        )
        proc.wait()
    os.chdir("../../..")
    return
