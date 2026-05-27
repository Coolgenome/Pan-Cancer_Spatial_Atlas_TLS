#!/bin/bash
#We used iSTAR to detect TLSs in the spatial transcriptomics samples. 
#The implementation was adapted and extended from the original iSTAR repository (https://github.com/daviddaiweizhang/istar).
#For details regarding the gene imputation step, please refer to the repository.

prefix=$1  # e.g. data/demo/
device="cuda"  # "cuda" or "cpu"
pixel_size=0.5
echo $pixel_size > ${prefix}pixel-size.txt
python rescale.py ${prefix} --image --mask --locs --radius
python preprocess.py ${prefix} --image --mask

n_genes=1000
gene_txt1="/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_signature/genes.txt"
gene_txt2="/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_signature/tls_genes_all.txt"
python ${pipeline}select_genes.py --n-top=${n_genes} "${prefix}cnts.tsv" "${prefix}gene-names.txt"
python ${pipeline}Add_genes.py ${prefix} ${gene_txt1}
python ${pipeline}Add_genes.py ${prefix} ${gene_txt2}
python ${pipeline}plot_spots.py ${prefix}

python ${pipeline}extract_features.py ${prefix} --device=${device}
python ${pipeline}impute.py ${prefix} --epochs=400 --device=${device}
python ${pipeline}plot_imputed.py ${prefix}

tls_folder="/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/tools/TLS_signature/"
struct_name_file=${tls_folder}structures.txt
taxonomy_file=${tls_folder}structures.yml

while read struct; do
    gene_name_file=${tls_folder}${struct}/gene-names.txt
    python ${pipeline}marker_score.py ${prefix} $gene_name_file ${prefix}markers/phenotype/raw/${struct}
done < $struct_name_file

python ${pipeline}phenotype.py ${prefix}

cp ${prefix}markers/phenotype/raw/lymphoid.pickle ${prefix}embeddings-tls.pickle

python ${pipeline}cluster_TLS.py --n-clusters=5 --filter-size=32 --min-cluster-size=20 ${prefix}embeddings-tls.pickle ${prefix}tls_cluster/
python ${pipeline}combine_TLS_cluster.py ${prefix} 5 4 3

python3 ${pipeline}enlarge.py ${prefix}

mask="tls_cluster/nclusters005/masks_combined/3.png"
python3 ${pipeline}tls_segment.py ${prefix} ${mask} TLS_segmentation/

anno="/rsrch6/home/genomic_med/lwang22_lab/kevin/kevin/4.TLS/2_Gastric_G35_v2/tls_annotation.csv"
python3 ${pipeline}tls_highlight.py ${prefix} ${mask} TLS_segmentation/ ${anno}
python3 ${pipeline}tls_segment_TOP.py ${prefix} ${mask} TLS_segmentation/ ${anno}
python3 ${pipeline}tls_segment_TOP_spatial_plot.py ${prefix} ${mask} TLS_segmentation/ ${anno}
python3 ${pipeline}tls_segment_TOP_six_genes.py ${prefix} ${mask} TLS_segmentation/ ${anno}
python3 ${pipeline}tls_segment_TOP_all_genes.py ${prefix} ${mask} TLS_segmentation/ ${anno}

python3 ${pipeline}save_pickle.py ${prefix} ${mask} TLS_segmentation/ ${anno}
python3 ${pipeline}pickle_to_spot_ST.py ${prefix} TLS_segmentation/combined_labeled_array.pickle TLS_segmentation/ ${anno}
