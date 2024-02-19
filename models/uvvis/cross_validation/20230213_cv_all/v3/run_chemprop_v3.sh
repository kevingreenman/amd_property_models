#!/bin/bash
#SBATCH -J uvvis_chemprop_v3
#SBATCH -o uvvis_chemprop_v3-%j.out
#SBATCH -t 5-00:00:00
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem=30gb
#SBATCH --gres=gpu:volta:1

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
cat $0
echo ""

source /etc/profile
module load anaconda/2022a
source activate chemprop

CHEMPROP_DIR='/home/gridsan/kgreenman/chemprop'
DATASETS_DIR='/home/gridsan/kgreenman/amd_property_models/datasets'
CONFIG_FILE_DIR='/home/gridsan/kgreenman/amd_property_models/models/uvvis'

python ../split.py --train_only_data $DATASETS_DIR/uvvis/v0/smiles_target_train.csv $DATASETS_DIR/uvvis/v0/smiles_target_val.csv --train_eval_data $DATASETS_DIR/uvvis/v1/v1.csv $DATASETS_DIR/uvvis/v2c/v2abc.csv $DATASETS_DIR/uvvis/v3/v3.csv --eval_only_data None --splits_dir splits_v3 --num_splits 5

for i in `seq 0 4`
do

# Train with cross-validation
python $CHEMPROP_DIR/train.py \
--data_path splits_v3/split_$i/trainval.csv \
--separate_test_path splits_v3/split_$i/eval.csv \
--dataset_type regression \
--save_dir $(pwd)/$i \
--metric rmse \
--extra_metrics mae \
--epochs 200 \
--gpu 0 \
--ensemble_size 5 \
--config_path $CONFIG_FILE_DIR/sigopt_chemprop_lambda_max_best_hyperparams_small.json \
--number_of_molecules 2 \

# Predict
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v0/smiles_target_test.csv --checkpoint_dir $(pwd)/$i --preds_path test_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v1/v1.csv --checkpoint_dir $(pwd)/$i --preds_path v1_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v2a/v2a.csv --checkpoint_dir $(pwd)/$i --preds_path v2a_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v2b/v2b.csv --checkpoint_dir $(pwd)/$i --preds_path v2b_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v2c/v2c.csv --checkpoint_dir $(pwd)/$i --preds_path v2c_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v2c/v2abc.csv --checkpoint_dir $(pwd)/$i --preds_path v2_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance
python $CHEMPROP_DIR/predict.py --test_path $DATASETS_DIR/uvvis/v3/v3.csv --checkpoint_dir $(pwd)/$i --preds_path v3_preds.csv --gpu 0 --number_of_molecules 2 --ensemble_variance

done

# Compress directory for new trained model, upload to Dropbox, and get shareable link
CURRENT_DIR=$(basename $(pwd))
cd ..
tar -czf $CURRENT_DIR.tar.gz $CURRENT_DIR/

