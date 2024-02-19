#!/bin/bash
#SBATCH -J photodeg_chemprop
#SBATCH -o photodeg_chemprop-%j.out
#SBATCH -t 5-00:00:00
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem=30gb
#SBATCH --gres=gpu:volta:1

echo "Date              = $$$$(date)"
echo "Hostname          = $$$$(hostname -s)"
echo "Working Directory = $$$$(pwd)"
echo ""
cat $$$$0
echo ""

source /etc/profile
module load anaconda/2022a
source activate chemprop

# if any data set is None, then enter "None". Can support space separated lists as an argument.
# train_data=$${@}
train_data="data/feb2022_expt.csv data/jun2022_expt.csv data/sep2021_expt.csv data/sep2022_expt.csv data/jul2022_expt.csv"
# location to save the splits and results
splits_dir="splits"
# number of splits to use in crossvalidation
num_splits=1
num_folds=25

qm_chemprop_dir=chemprop-atom-bond
chemprop_dir=chemprop
photodeg_model_dir=.
scratch=scratch_photodeg

hidden_size=1650
ffn_hidden_size=2400
depth=2
ffn_num_layers=3
dropout=0
batch_size=5
epochs=100

python $$photodeg_model_dir/split.py \
--train_data $$train_data \
--num_splits $$num_splits \
--splits_dir $$splits_dir

for j in $$(eval echo "{1..$$num_splits}")
do
    split=$$j

    python $$photodeg_model_dir/index_smiles.py \
    $$splits_dir/split_$$split/trainval.csv \
    $$scratch/index_input.csv

    python $$qm_chemprop_dir/predict.py \
    --test_path $$scratch/index_input.csv \
    --checkpoint_path $$qm_chemprop_dir/trained_model/QM_137k.pt \
    --preds_path $$scratch/qm_features.pkl

    python $$photodeg_model_dir/split_features.py \
    $$scratch/qm_features.pkl \
    $$splits_dir/split_$$split/trainval_atom_features.pkl \
    $$splits_dir/split_$$split/trainval_bond_features.pkl

    python $$photodeg_model_dir/index_smiles.py \
    $$splits_dir/split_$$split/eval.csv \
    $$scratch/index_input.csv

    python $$qm_chemprop_dir/predict.py \
    --test_path $$scratch/index_input.csv \
    --checkpoint_path $$qm_chemprop_dir/trained_model/QM_137k.pt \
    --preds_path $$scratch/qm_features.pkl

    python $$photodeg_model_dir/split_features.py \
    $$scratch/qm_features.pkl \
    $$splits_dir/split_$$split/eval_atom_features.pkl \
    $$splits_dir/split_$$split/eval_bond_features.pkl

    sub_dir=iter$$i

    python $$chemprop_dir/train.py \
    --data_path $$splits_dir/split_$$split/trainval.csv \
    --split_sizes 0.89 0.11 0 \
    --save_dir $$splits_dir/split_$$split \
    --dataset_type regression \
    --batch_size $$batch_size --epochs $$epochs \
    --depth $$depth --ffn_num_layers $$ffn_num_layers \
    --hidden_size $$hidden_size --ffn_hidden_size $$ffn_hidden_size \
    --dropout $$dropout \
    --num_folds $$num_folds \
    --atom_descriptors feature \
    --atom_descriptors_path $$splits_dir/split_$$split/trainval_atom_features.pkl \
    --bond_features_path $$splits_dir/split_$$split/trainval_bond_features.pkl \
    --adding_h

    python $$chemprop_dir/predict.py \
    --checkpoint_dir $$splits_dir/split_$$split \
    --test_path $$splits_dir/split_$$split/eval.csv \
    --preds_path $$splits_dir/split_$$split/eval_preds.csv \
    --atom_descriptors feature \
    --atom_descriptors_path $$splits_dir/split_$$split/eval_atom_features.pkl \
    --bond_features_path $$splits_dir/split_$$split/eval_bond_features.pkl \
    --ensemble_variance
done

python $$photodeg_model_dir/gather_eval.py \
--splits_dir $$splits_dir \
--num_splits $$num_splits \
--train_data $$train_data
