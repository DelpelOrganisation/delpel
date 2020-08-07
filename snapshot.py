import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from myshow import myshow, myshow3d
import os
import json
from fpdf import FPDF




def create_snapshop(fileName_CT, fileName_segmentation, output_dir, patient, study, struct_name):
    
    image_CT = sitk.ReadImage(fileName_CT)
    image_seg = sitk.ReadImage(fileName_segmentation)
    
    ## resample segmentation image with CT image
    identity = sitk.Transform(3, sitk.sitkIdentity)
    image_seg_resampled = sitk.Resample(image_seg,image_CT.GetSize(),identity,sitk.sitkNearestNeighbor,image_CT.GetOrigin(),image_CT.GetSpacing(),image_CT.GetDirection())
    
    ## get center of the segmentation (in physical space, then in image indices)
    statistics_label_filter = sitk.LabelShapeStatisticsImageFilter()
    statistics_label_filter.Execute(image_seg_resampled)
    centroid = statistics_label_filter.GetCentroid(1)
#    print("Centroid : ", centroid)
    centroid_indices = list(centroid)
    for i in range(0,3):
        centroid_indices[i] = int((centroid[i]-image_CT.GetOrigin()[i])/image_CT.GetSpacing()[i])
    
    ## cast CT image to char
    image_CT_256 = sitk.Cast(sitk.IntensityWindowing(image_CT,windowMinimum=20-200, windowMaximum=20+200), sitk.sitkUInt8)
        
    overlay_img = sitk.LabelOverlay(image_CT_256, image_seg_resampled, opacity = 0.15)
    
#    seg = sitk.Cast(image_seg_resampled, sitk.sitkLabelUInt8)
#    
#    overlay_img_contour = sitk.LabelMapContourOverlay(seg, image_CT_256,opacity = 0.5, 
#                                                 contourThickness=[1,1,1]
#                                              )
#    
    img_xslices = overlay_img[centroid_indices[0], :, :]
    img_yslices = overlay_img[:, centroid_indices[1], :]
    img_zslices = overlay_img[:, :, centroid_indices[2]]
    

    myshow(img_xslices, title=patient + " " + struct_name, invertX=True, zoom=2.5, save=output_dir+patient+"_"+study+"_"+struct_name+"1.png", figsize=(7,7))
    myshow(img_yslices, title=patient + " " + struct_name, invertX=True, zoom=2.5, save=output_dir+patient+"_"+study+"_"+struct_name+"2.png", figsize=(7,7))
    myshow(img_zslices, title=patient + " " + struct_name, invertY=True, save=output_dir+patient+"_"+study+"_"+struct_name+"3.png", figsize=(7,7))
    
    plt.close("all")



directory = "/home/tbaudier/david/delpel/patient/saphir/output/"

 
        
with open("/home/tbaudier/david/delpel/delpel.json") as f:
  data = json.load(f)

#print(json.dumps(data, indent=4))

patients = list(data.get("patients"))


pdf = FPDF()
pdf.add_page()
pdf.set_font('Arial', 'B', 16)

## coordonnées des images dans pdf
x_start = 10
x_delta = 67
y_start = 40
y_delta = 67
size_x = 65
size_y = 65


for patient in patients:
    print(patient)
   
    pdf.write(10, "Patient : " + patient +"\n")

    data_patient = data.get("patients")[patient]
    for study, study_data in data_patient.items():
        print(study)
        pdf.write(10, "Study : " + study)
        structure_data = study_data.get("Structures")
        found_PTV=False
       
        for struct, struct_reference in structure_data.items():
        
            if struct_reference=="Bladder_ext":
                print(struct)
                directory_study = directory+"/patient."+patient+"/input/"+study
                fileName_CT = directory_study+"/CT.nii"
                fileName_segmentation = directory_study+"/"+struct+".nii"
                create_snapshop(fileName_CT, fileName_segmentation, "snapshots/",patient, study, struct_reference)
                
                ## PDF
                x = x_start
                y = y_start
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"3.png", x,y,size_x,size_y)
                x = x+x_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"2.png", x,y,size_x,size_y)
                x = x+x_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"1.png", x,y,size_x,size_y)
                
            if struct_reference=="Rectum_ext":
                print(struct)
                directory_study = directory+"/patient."+patient+"/input/"+study
                fileName_CT = directory_study+"/CT.nii"
                fileName_segmentation = directory_study+"/"+struct+".nii"
                create_snapshop(fileName_CT, fileName_segmentation, "snapshots/",patient, study, struct_reference)
                
                ## PDF
                x = x_start
                y = y_start+ y_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"3.png", x,y,size_x,size_y)
                x = x+x_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"2.png", x,y,size_x,size_y)
                x = x+x_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"1.png", x,y,size_x,size_y)
                
                
            if (struct_reference=="CTV_prostate" or struct_reference=="Prostate") and found_PTV==False:
                found_PTV=True
                print(struct)
                directory_study = directory+"/patient."+patient+"/input/"+study
                fileName_CT = directory_study+"/CT.nii"
                fileName_segmentation = directory_study+"/"+struct+".nii"
                create_snapshop(fileName_CT, fileName_segmentation, "snapshots/",patient, study, struct_reference)
                
                ## PDF
                x = x_start
                y = y_start+ 2*y_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"3.png", x,y,size_x,size_y)
                x = x+x_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"2.png", x,y,size_x,size_y)
                x = x+x_delta
                pdf.image("snapshots/"+patient+"_"+study+"_"+struct_reference+"1.png", x,y,size_x,size_y)
        
        pdf.add_page()       
                
pdf.output('snapshots.pdf', 'F')


        
        
        
    
    
    

    
    
    


#directory = 'E:/data_creatis_old/patient.21272778/input/1.2.840.113854.131546242490825631409557113510065697332.1/'
#
#fileName_CT = "CT.nii"
#fileName_segmentation = "8_VessieExt.nii"
#
#output_dir = "snapshots/"
#if not os.path.exists(output_dir):
#    os.makedirs(output_dir)
#
#create_snapshop(directory + fileName_CT, directory + fileName_segmentation, output_dir, "bladder")
