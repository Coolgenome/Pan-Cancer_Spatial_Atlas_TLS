import os
from ultralytics import YOLO
import sys

model_path = os.path.join('./data/TLS_classification_HE_YOLO_step2.pt')
model=YOLO(model_path)

image_path = sys.argv[1]
out_dir = sys.argv[2]
model(source=image_path,project=out_dir,save_txt=True,imgsz=640)
