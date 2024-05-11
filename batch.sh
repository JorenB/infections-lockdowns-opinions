#!/bin/bash
#SBATCH --job-name=mp
#SBATCH --partition=hef
#SBATCH --time=04-00:00:00
#SBATCH --array=0-9
#SBATCH --output=log/measure_%A_%a.out
#SBATCH --error=err/error_%A_%a.err
#SBATCH --mem=20G

#ids=("1")
ids=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#ids=("16" "17" "18" "19" "20")
#ids=("11" "12" "13" "14" "15" "16" "17" "18" "19" "20")
#ids=("21" "22" "23" "24" "25" "26" "27" "28" "29" "30")
#ids=("31" "32" "33" "34" "35" "36" "37" "38" "39" "40")
#ids=("41" "42" "43" "44" "45" "46" "47" "48" "49" "50")
#ids=("51" "52" "53" "54" "55" "56" "57" "58" "59" "60")
#ids=("61" "62" "63" "64" "65" "66" "67" "68" "69" "70")
#ids=("71" "72" "73" "74" "75")
# ids=("1" "2")
echo "batch job " $SLURM_ARRAY_TASK_ID
srun -u matlab -nodisplay -r "run('main(${ids[SLURM_ARRAY_TASK_ID]})'); exit;"
