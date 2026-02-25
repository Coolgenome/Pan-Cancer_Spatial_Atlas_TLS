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
gene_txt="./data/genes_add.txt"
#need clean genes_add
python select_genes.py --n-top=${n_genes} "${prefix}cnts.tsv" "${prefix}gene-names.txt"
python Add_genes.py ${prefix} ${gene_txt}
python plot_spots.py ${prefix}

python extract_features.py ${prefix} --device=${device}
python impute.py ${prefix} --epochs=400 --device=${device}
python plot_imputed.py ${prefix}


tls_folder="./data/TLS_signature/"
struct_name_file=${tls_folder}structures.txt
taxonomy_file=${tls_folder}structures.yml

while read struct; do
    gene_name_file=${tls_folder}${struct}/gene-names.txt
    python marker_score.py ${prefix} $gene_name_file ${prefix}markers/phenotype/raw/${struct}
done < $struct_name_file

python phenotype.py ${prefix}

cp ${prefix}markers/phenotype/raw/lymphoid.pickle ${prefix}embeddings-tls.pickle

python cluster_TLS.py --n-clusters=5 --filter-size=32 --min-cluster-size=20 ${prefix}embeddings-tls.pickle ${prefix}tls_quantify_filtersize_32/ncluster_5/size_20/
python combine_TLS_cluster.py ${prefix} 5 32 4 3

python preprocess_tls.py ${prefix} 5 32
python preprocess_tls_combined.py ${prefix} 5 32

python3 enlarge.py ${prefix}

mask="tls_quantify_filtersize_32/ncluster_5/size_20/nclusters005/masks_combined/3.png"

python3 tls_segment.py ${prefix} ${mask} TLS_segmentation/
