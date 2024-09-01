#!/bin/bash
#SBATCH --job-name="ldsr_fpQTLs"
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

source ~/.bashrc
source activate ldsc

LDSC="/mnt/isilon/sfgi/programs/ldsc"


# hg38
cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
BFILE=$cdir/plink_files_hg38/1000G.EUR.hg38
FREQ=$cdir/plink_files_hg38/1000G.EUR.hg38.
BASELINE=$cdir/baselineLD_v2.2_hg38/baselineLD.
WEIGHTS=$cdir/weights_hg38/weights.hm3_noMHC.
PRINT_SNPS=$cdir/list.txt

# echo "plink_files_hg38/1000G.EUR.hg38.1.bim:"
# wc -l $BFILE.1.bim
# echo ""

# echo "baselineLD_v2.2_hg38/baselineLD.1.annot.gz:"
# gunzip < "${BASELINE}1.annot.gz" | wc -l 
# echo ""

# echo "weights_hg38/weights.hm3_noMHC.1.l2.ldscore.gz:"
# gunzip < "${WEIGHTS}1.l2.ldscore.gz" | wc -l 
# echo ""

# echo "plink_files_hg38/1000G.EUR.hg38.1.frq:"
# wc -l "${FREQ}1.frq"
# echo ""

# Run to reset results
# for n in annotations/* ; do
#     arrIN=(${n//\// })
#     NAME=${arrIN[-1]}
#     
#     rm annotations/$NAME/out/*
# done

# for n in annotations/* ; do
#     arrIN=(${n//\// })
#     NAME=${arrIN[-1]}
    
for NAME in PRINT_no_gaussian PRINT_no_gaussian_FDR10; do

    BED_IN="annotations/$NAME/$NAME.bed"
    LD_SCORE="annotations/$NAME/ld_score"

    mkdir -p annotations/$NAME/ld_score
    mkdir -p annotations/$NAME/out

    echo "Running LDSR for functional annotation '$NAME'"

    # Create annotation files
    if [ -f "$LD_SCORE/$NAME.22.l2.ldscore.gz" ]; then
        echo "Skipping annotation creation for $NAME because $LD_SCORE/$NAME.22.l2.ldscore.gz already exists..." 
    else
        # Bed to Annot
        for j in {1..22}; do
        python $LDSC/make_annot.py \
            --bed-file $BED_IN \
            --bimfile $BFILE.$j.bim \
            --annot-file $LD_SCORE/$NAME.$j.annot.gz
        done

        #Annot to LD files
        for j in {1..22}; do
        python $LDSC/ldsc.py \
            --l2 \
            --bfile $BFILE.$j \
            --ld-wind-cm 1 \
            --annot $LD_SCORE/$NAME.$j.annot.gz \
            --thin-annot \
            --out $LD_SCORE/$NAME.$j \
            --print-snps $PRINT_SNPS
        done
    fi

    # Run LDSR, liver curated sumstats
    SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/munged_sumstats"
    for sumstat in "$SUMSTATS"/*.sumstats.gz; do
        echo "Running regression on $sumstat..."

        # */NAFLD.sumstats.gz --> NAFLD
        arrIN=(${sumstat//\// })
        arrIN=(${arrIN[-1]//.sumstats./ })
        trait=${arrIN[0]}

        OUT="annotations/$NAME/out/$trait.$NAME"

        if [ -f "$OUT.results" ]; then
            echo "Skipping because $OUT.results already exists..."
            continue 
        fi

        # #run ldscore regression
        # #h2 outpt from munge sumstats
        python $LDSC/ldsc.py \
            --h2 $sumstat \
            --ref-ld-chr $LD_SCORE/$NAME.,$BASELINE \
            --w-ld-chr $WEIGHTS \
            --overlap-annot \
            --frqfile-chr $FREQ \
            --print-coefficients \
            --out $OUT
    done

    # Alkes Price group sumstats
    SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/alkesgroup_sumstats"
    for sumstat in "$SUMSTATS"/*.sumstats; do
        echo "Running regression on $sumstat..."

        # */NAFLD.sumstats --> NAFLD
        arrIN=(${sumstat//\// })
        arrIN=(${arrIN[-1]//.sumstats/ })
        trait=${arrIN[0]}

        OUT="annotations/$NAME/out/$trait.$NAME"

        if [ -f "$OUT.results" ]; then
            echo "Skipping because $OUT.results already exists..."
            continue 
        fi

        # #run ldscore regression
        # #h2 outpt from munge sumstats
        python $LDSC/ldsc.py \
            --h2 $sumstat \
            --ref-ld-chr $LD_SCORE/$NAME.,$BASELINE \
            --w-ld-chr $WEIGHTS \
            --overlap-annot \
            --frqfile-chr $FREQ \
            --print-coefficients \
            --out $OUT
    done
done
