from ultralytics import YOLO

if __name__ == '__main__':
    # Load the 'Small' segmentation model
    model = YOLO('yolov8s-seg.pt') 

    # Start training
    model.train(
        data='E:/EucVision_YOLO/data.yaml',
        epochs=150,           
        imgsz=640,
        batch=8,              # 8 is the limit for 6GB VRAM at 640px
        workers=4,            # Balances the CPU data-loading threads
        save_period=10,       # Save a checkpoint every 10 epochs
        overlap_mask=False,   # Set to False to help keep individual tree boundaries distinct
        project='EucVision_Thesis',
        name='Eucalyptus_YOLO_v1',
        device=0
    )
    
    
# To start again
    
import os
# This line tells Windows to ignore the OpenMP library conflict
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

from ultralytics import YOLO

if __name__ == '__main__':
    # 1. Point YOLO to the LAST saved checkpoint
    model = YOLO('C:/Users/jakev/Documents/Python Project/runs/segment/EucVision_Thesis/Eucalyptus_YOLO_v13/weights/last.pt') 

    # 2. Tell it to resume
    model.train(resume=True)