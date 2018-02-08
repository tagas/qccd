# Radar(Spider) plots for QCCD vs. others

library(fmsb)
pdf("spider_plots.pdf", width = 11, height = 6)
par(mfrow = c(2,4),mai = c(1, 0.6, 0.1, 0.6),mar=c(1, 0.5, 1, 0.1), cex.main = 1.6)


col_names = c("SIM", "SIM \n -ln", "SIM \n -G", "AN \n -G", "AN \n -s",
              "HN", "HN \n -s", "MN \n -U", "MN \n -G", "Tueb")
p_col = rgb(0.529412, 0.807843, 0.980392,0.7)
pf_col = rgb(0.678431, 0.847059, 0.901961,0.7)


# QCCD
cc_acc_df = as.data.frame(matrix(c(0.49,0.77,0.76,1,1,1,1,1, 0.82, 0.66), ncol = 10))
colnames(cc_acc_df) = col_names
cc_acc_df = rbind(rep(1,10) , rep(0,10) , cc_acc_df) 


radarchart( cc_acc_df  , axistype=1 , 
            pcol = p_col , pfcol = pf_col, plwd=4, 
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,
            vlcex=1.5 , title = "QCCD")

# CAM
cam_acc_df = as.data.frame(matrix(c(0.57,0.87,0.81,1,1,0.78,0.47, 0.04, 0.03, 0.57), ncol = 10))
colnames(cam_acc_df) =  col_names
cam_acc_df = rbind(rep(1,10) , rep(0,10) , cam_acc_df) 

radarchart( cam_acc_df  , axistype=1 , 
            pcol=p_col , pfcol=pf_col, plwd=4,
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
            ,title = "biCAM")


# IGCI
igci_acc_df = as.data.frame(matrix(c(0.42,0.52,0.54,0.41,0.35,0.81,0.85, 0.02, 0.97, 0.67), ncol = 10))
colnames(igci_acc_df) = col_names
igci_acc_df = rbind(rep(1,10) , rep(0,10) , igci_acc_df) 

radarchart( igci_acc_df  , axistype=1 , 
            pcol=p_col , pfcol=pf_col, plwd=4,
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
            ,title = "IGCI")

# Slope
slope_acc_df = as.data.frame(matrix(c(0.45,0.47,0.48,0.18, 0.28, 0.62, 0.74, 0.07, 0.96, 0.74), ncol = 10))
colnames(slope_acc_df)= col_names
slope_acc_df = rbind(rep(1,10) , rep(0,10) , slope_acc_df) 

radarchart( slope_acc_df  , axistype=1 , 
            pcol=p_col , pfcol=pf_col, plwd=4,
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
            ,title = "Slope")
s

resit_acc_df = as.data.frame(matrix(c(0.78,0.87,0.77,1,1, 0.31,0.08, 0.09,0.04,0.53), ncol = 10))
colnames(resit_acc_df) = col_names
resit_acc_df = rbind(rep(1,10) , rep(0,10) , resit_acc_df) 

# RESIT
radarchart(resit_acc_df  , axistype=1 , 
           pcol=p_col , pfcol=pf_col, plwd=4,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
           ,title = "RESIT")

# LINGAM
lin_acc_df = as.data.frame(matrix(c(0.42,0.29,0.33,0.04,0.04,0,0.01,0,0.02,0.3), ncol = 10))
colnames(lin_acc_df) = col_names
lin_acc_df = rbind(rep(1,10), rep(0,10), lin_acc_df) 

radarchart(lin_acc_df  , axistype=1 , 
           pcol=p_col , pfcol=pf_col, plwd=4,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
           ,title = "LINGAM")


# EMD
emd_acc_df = as.data.frame(matrix(c(0.47,0.52,0.6, 0.62, 0.47, 0.73,0.84,0.62,0.95, 0.61), ncol = 10))
colnames(emd_acc_df) = col_names
emd_acc_df = rbind(rep(1,10), rep(0,10), emd_acc_df)

radarchart(emd_acc_df  , axistype=1 , 
           pcol=p_col , pfcol=pf_col, plwd=4,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
           ,title = "EMD")


# GR-AN
gran_acc_df = as.data.frame(matrix(c(0.47,0.41,0.41,0.05,0.06,0.39,0.49, 0.94,0.49, 0.40), ncol = 10))
colnames(gran_acc_df) = col_names
gran_acc_df = rbind(rep(1,10), rep(0,10), gran_acc_df) 

radarchart(gran_acc_df  , axistype=1 , 
           pcol=p_col , pfcol=pf_col, plwd=4,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,5), cglwd=0.8,vlcex=1.5
           ,title = "GR-AN")
dev.off()

