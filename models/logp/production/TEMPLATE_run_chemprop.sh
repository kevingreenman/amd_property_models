#!/bin/bash
#SBATCH -J logp_chemprop
#SBATCH -o logp_chemprop-%j.out
#SBATCH -t 3-00:00:00
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

# Train
datasets_dir='/home/gridsan/kgreenman/amd_property_models/datasets/logp'
logp_dir='/home/gridsan/kgreenman/amd_property_models/models/logp'
# location of the modified chemprop directory to include relevant features. Do not use your normal chemprop directory from github.
# Comes from https://github.com/cjmcgill/chemprop.git branch "logp".
chemprop_dir=$${logp_dir}/chemprop
# location to save the splits and results.
model_dir=$$(pwd)
# data to be used for training the model. List of files, e.g. "v0.csv v1.csv v2.csv"
train_data=$TRAIN_DATA_FILES_LIST_FROM_PYTHON

depth=3
hidden_size=1100
ffn_num_layers=3
ffn_hidden_size=1100
dropout=0.05
activation=PReLU
final_lr=0.00013402665639290512
init_lr=5.515076512972011e-07
max_lr=0.0014333625202976184
warmup_epochs=25
epochs=100
num_folds=5
ensemble_size=3
aggregation=norm

python $$logp_dir/split.py \
--train_only_data $$train_data \
--train_eval_data None \
--eval_only_data None \
--num_splits 1 \
--splits_dir $$model_dir
--data_dir $$datasets_dir

split_dir=$$model_dir/split_0

python $$logp_dir/make_features.py \
$$split_dir/trainval.csv \
$$split_dir/trainval_features.csv

python $$chemprop_dir/train.py \
--data_path $$split_dir/trainval.csv \
--features_path $$split_dir/trainval_features.csv \
--dataset_type regression \
--save_dir $$split_dir \
--split_sizes 0.89 0.11 0 \
--depth $$depth \
--hidden_size $$hidden_size \
--ffn_num_layers $$ffn_num_layers \
--ffn_hidden_size $$ffn_hidden_size \
--dropout $$dropout \
--activation $$activation \
--final_lr $$final_lr \
--init_lr $$init_lr \
--max_lr $$max_lr \
--warmup_epochs $$warmup_epochs \
--epochs $$epochs \
--num_folds $$num_folds \
--ensemble_size $$ensemble_size \
--aggregation $$aggregation


# Predict
# csv file with the smiles you want to make predictions for
smiles_path=$${datasets_dir}/v2/v2.csv
# path to save the predictions
preds_path='test_preds.csv'

python $$logp_dir/make_features.py \
$$smiles_path \
$$model_dir/test_features.csv

python $$chemprop_dir/predict.py \
--test_path $$smiles_path \
--features_path $$model_dir/test_features.csv \
--preds_path $$preds_path \
--checkpoint_dir $$model_dir \
--ensemble_variance
