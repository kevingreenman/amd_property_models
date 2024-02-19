#!/bin/bash
#SBATCH -J uvvis_chemprop
#SBATCH -o uvvis_chemprop-%j.out
#SBATCH -t 1-00:00:00
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem=30gb
#SBATCH --gres=gpu:volta:1

echo "Date              = $$(date)"
echo "Hostname          = $$(hostname -s)"
echo "Working Directory = $$(pwd)"
echo ""
cat $$0
echo ""

source /etc/profile
module load anaconda/2022a
source activate chemprop

CHEMPROP_DIR='/home/gridsan/kgreenman/chemprop'
DATASETS_DIR='/home/gridsan/kgreenman/amd_property_models/datasets'
CONFIG_FILE_DIR='/home/gridsan/kgreenman/amd_property_models/models/uvvis'

# Train and test an ensemble of Chemprop models
python $$CHEMPROP_DIR/train.py --data_path $$DATASETS_DIR/$TRAIN_DATA_FILE_FROM_PYTHON --separate_val_path $$DATASETS_DIR/uvvis/v0/smiles_target_val.csv --separate_test_path $$DATASETS_DIR/uvvis/v0/smiles_target_test.csv --dataset_type regression --save_dir $$(pwd) --metric rmse --epochs 200 --gpu 0 --ensemble_size 5 --config_path $$CONFIG_FILE_DIR/sigopt_chemprop_lambda_max_best_hyperparams_small.json --number_of_molecules 2

python $$CHEMPROP_DIR/predict.py --test_path $$DATASETS_DIR/uvvis/v0/smiles_target_test.csv --checkpoint_dir $$(pwd) --preds_path preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance


# Compress directory for new trained model, upload to Dropbox, and get shareable link
CURRENT_DIR=$$(basename $$(pwd))
cd ..
tar -czf $$CURRENT_DIR.tar.gz $$CURRENT_DIR/
