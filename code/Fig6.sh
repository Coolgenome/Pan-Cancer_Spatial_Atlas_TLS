#!/bin/bash
#We use HookNet-TLS to detect TLSs in the H&E samples.
#Please refer to the original paper for more details: https://www.nature.com/articles/s43856-023-00421-7

python3 -m hooknettls \
    hooknettls.default.image_path=${image_path} \
    hooknettls.default.mask_path=${mask_path}

patch1_path = os.path.join(./data/, "TLS1.jpeg")
patch1_step1_out_dir = "./result/TLS1_step1/"
patch1_step2_out_dir = "./result/TLS1_step2/"

python3 yolo_step1.py ${patch1_path} ${patch1_step1_out_dir}
python3 yolo_step2.py ${patch1_path} ${patch1_step2_out_dir}