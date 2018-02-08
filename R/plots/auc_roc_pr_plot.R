# AUC-ROC, AUC-PR and plot for roc curves

rm(list = ls())
library(ROCR)
library(PRROC)
library(readr)
library(naglr)
source("./cd_utils.R")

QCCD_AN_s <- read_delim("./../results/QCCD_AN_s.tab",
                                "\t", escape_double = FALSE, trim_ws = TRUE)
pairs_gt_s <- read_csv("./../data/ANHNMN_pairs/AN-s/pairs_gt.txt", 
                       col_names = FALSE)
pred <- prediction(QCCD_AN_s$Eps, pairs_gt_s$X1)
perf <- performance(pred,"tpr","fpr")
auc_AN_s <- performance(pred, "auc")@y.values[[1]] 
pr.curve(scores.class0 = QCCD_AN_s$Eps, weights.class0 = pairs_gt_s$X1, curve = T);
plot(perf, col = "blue", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type="l", lty=2)


# double check ROC/PR
# ROC

ROCR::performance(ROCR::prediction(QCCD_AN_s$Eps, pairs_gt_s$X1), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_AN_s$Eps, weights.class0 = pairs_gt_s$X1)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_AN_s$Eps, pairs_gt_s$X1), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_AN_s$Eps, weights.class0 = pairs_gt_s$X1)


# need these for the ggplot graph
auc_an_s_x = unlist(perf@x.values)
auc_an_s_y = unlist(perf@y.values)


QCCD_AN <- read_delim("./../results/QCCD_AN.tab",
                               "\t", escape_double = FALSE, trim_ws = TRUE)
pairs_gt_AN <- read_csv("./../data/ANHNMN_pairs/AN/pairs_gt.txt", 
                         col_names = FALSE)
pred_AN <- prediction(QCCD_AN$Eps, pairs_gt_AN)
perf_AN <- performance(pred_AN,"tpr","fpr")
auc_AN <- performance(pred_AN, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_AN$Eps, weights.class0 = pairs_gt_AN$X1, curve = T);
par(new = TRUE)
plot(perf_AN, col = "steelblue", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type="l", lty=3)

auc_an_x = unlist(perf_AN@x.values)
auc_an_y = unlist(perf_AN@y.values)


# ROC
ROCR::performance(ROCR::prediction(QCCD_AN$Eps, pairs_gt_AN$X1), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_AN$Eps, weights.class0 = pairs_gt_AN$X1)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_AN$Eps, pairs_gt_AN$X1), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_AN$Eps, weights.class0 = pairs_gt_AN$X1)



QCCD_H <- read_delim("./../results/QCCD_HN.tab",
                             "\t", escape_double = FALSE, trim_ws = TRUE)
pairs_gt_h <- read_csv("./../data/ANHNMN_pairs/HN/pairs_gt.txt", 
                       col_names = FALSE)
pred_h <- prediction(QCCD_H$Eps, pairs_gt_h)
perf_h <- performance(pred_h,"tpr","fpr")
auc_h <- performance(pred_h, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_H$Eps, weights.class0 = pairs_gt_h$X1, curve = T);
par(new = TRUE)
plot(perf_h, col = "mediumblue", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type="l", lty=4)

auc_h_x = unlist(perf_h@x.values)
auc_h_y = unlist(perf_h@y.values)


# ROC
ROCR::performance(ROCR::prediction(QCCD_H$Eps, pairs_gt_h$X1), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_H$Eps, weights.class0 = pairs_gt_h$X1)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_H$Eps, pairs_gt_h$X1), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_H$Eps, weights.class0 = pairs_gt_h$X1)


QCCD_HNs <- read_delim("./../results/QCCD_HN_s.tab",
                              "\t", escape_double = FALSE, trim_ws = TRUE)
pairs_gt_HNs <- read_csv("./../data/ANHNMN_pairs/HN-s/pairs_gt.txt", 
                        col_names = FALSE)
pred_HN <- prediction(QCCD_HNs$Eps, pairs_gt_HNs)
perf_HN <- performance(pred_HN,"tpr","fpr")
auc_HN <- performance(pred_HN, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_HNs$Eps, weights.class0 = pairs_gt_HNs$X1, curve = T);
par(new = TRUE)
plot(perf_HN, col = "slateblue2", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type="l", lty=4)

auc_hn_x = unlist(perf_HN@x.values)
auc_hn_y = unlist(perf_HN@y.values)

# ROC
ROCR::performance(ROCR::prediction(QCCD_HNs$Eps, pairs_gt_HNs$X1), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_HNs$Eps, weights.class0 = pairs_gt_HNs$X1)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_HNs$Eps, pairs_gt_HNs$X1), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_HNs$Eps, weights.class0 = pairs_gt_HNs$X1)



QCCD_sim <- read_delim("./../results/QCCD_SIM.tab", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
path = "../data/SIM_pairs/"
ext = "SIM/" 
pairmeta <- read.table(paste0(path, ext,"pairmeta.txt"),quote = "\"")
file.names <- dir(paste0(path,ext), pattern = ".txt")

ce_gt <- data.frame(pairmeta[, 2], pairmeta[, 5]) # the ground truth
ce_gt = ce_gt[,1]
n_pairs = 100
ce_gt[ce_gt == 2] = 0

pairs_gt_sim <- ce_gt

pred_sim <- prediction(QCCD_sim$Eps, pairs_gt_sim)
perf_sim <- performance(pred_sim,"tpr","fpr")
auc_sim <- performance(pred_sim, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_sim$Eps, weights.class0 = pairs_gt_sim, curve = T);
par(new = TRUE)
plot(perf_sim, col = "lightsalmon4", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type = "l", lty=1)

auc_sim_x = unlist(perf_sim@x.values)
auc_sim_y = unlist(perf_sim@y.values)



# ROC
ROCR::performance(ROCR::prediction(QCCD_sim$Eps,pairs_gt_sim), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_sim$Eps, weights.class0 = pairs_gt_sim)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_sim$Eps, pairs_gt_sim), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_sim$Eps, weights.class0 = pairs_gt_sim)



QCCD_simln <- read_delim("./../results/QCCD_SIM-ln.tab", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
path = "../data/SIM_pairs/"
ext = "SIM-ln/"
pairmeta <- read.table(paste0(path, ext,"pairmeta.txt"),quote = "\"")
file.names <- dir(paste0(path,ext), pattern = ".txt")

ce_gt <- data.frame(pairmeta[, 2], pairmeta[, 5]) # the ground truth
ce_gt = ce_gt[,1]
n_pairs = 100
ce_gt[ce_gt == 2] = 0

pairs_gt_simln <- ce_gt
QCCD_simln$Eps[is.na(QCCD_simln$Eps)]= 0.5
pred_simln <- prediction(QCCD_simln$Eps, pairs_gt_simln)
perf_simln <- performance(pred_simln,"tpr","fpr")
auc_simln <- performance(pred_simln, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_simln$Eps, weights.class0 = pairs_gt_simln, curve = T);
par(new = TRUE)
plot(perf_simln, col = "orange", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type = "l")

auc_simln_x = unlist(perf_simln@x.values)
auc_simln_y = unlist(perf_simln@y.values)


# ROC
ROCR::performance(ROCR::prediction(QCCD_simln$Eps,pairs_gt_simln), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_simln$Eps, weights.class0 = pairs_gt_simln)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_simln$Eps, pairs_gt_simln), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_simln$Eps, weights.class0 = pairs_gt_simln)


QCCD_simg <- read_delim("./../results/QCCD_SIM-G.tab", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
path = "../data/SIM_pairs/"
ext = "SIM-G/" 
pairmeta <- read.table(paste0(path, ext,"pairmeta.txt"),quote = "\"")
file.names <- dir(paste0(path,ext), pattern = ".txt")

ce_gt <- data.frame(pairmeta[, 2], pairmeta[, 5]) # the ground truth
ce_gt = ce_gt[,1]
n_pairs = 100
ce_gt[ce_gt == 2] = 0

pairs_gt_simg <- ce_gt
QCCD_simg$Eps[is.na(QCCD_simg$Eps)]= 0.5
pred_simg <- prediction(QCCD_simg$Eps, pairs_gt_simg)
perf_simg <- performance(pred_simg,"tpr","fpr")
auc_simg <- performance(pred_simg, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_simg$Eps, weights.class0 = pairs_gt_simg, curve = T);
par(new = TRUE)
plot(perf_simg, col = "lightsalmon1", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type = "l")


auc_simg_x = unlist(perf_simg@x.values)
auc_simg_y = unlist(perf_simg@y.values)




# ROC
ROCR::performance(ROCR::prediction(QCCD_simg$Eps,pairs_gt_simg), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_simg$Eps, weights.class0 = pairs_gt_simg)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_simg$Eps, pairs_gt_simg), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_simg$Eps, weights.class0 = pairs_gt_simg)



QCCD_MN_s_un <- read_delim("./../results/QCCD_MN_u.tab",
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
pairs_gt_umn <- read_csv("./../data/ANHNMN_pairs/MN-U/pairs_gt.txt", 
                         col_names = FALSE)
pred <- prediction(QCCD_MN_s_un$Eps, pairs_gt_umn$X1)
perf <- performance(pred,"tpr","fpr")
auc_umn <- performance(pred, "auc")@y.values[[1]]
pr.curve(scores.class0 = QCCD_MN_s_un$Eps, weights.class0 = pairs_gt_umn$X1, curve = T);
par(new = TRUE)
plot(perf, col = "indianred1", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type="l", lty=2)

auc_unm_x = unlist(perf@x.values)
auc_unm_y = unlist(perf@y.values)


# ROC
ROCR::performance(ROCR::prediction(QCCD_MN_s_un$Eps,pairs_gt_umn$X1), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_MN_s_un$Eps, weights.class0 = pairs_gt_umn$X1)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_MN_s_un$Eps, pairs_gt_umn$X1), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_MN_s_un$Eps, weights.class0 = pairs_gt_umn$X1)



QCCD_MN_gmn <- read_delim("./../results/QCCD_MN_g.tab",
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
pairs_gt_gmn <- read_csv("./../data/ANHNMN_pairs/MN-G/pairs_gt.txt", 
                         col_names = FALSE)
pred <- prediction(QCCD_MN_gmn$Eps, pairs_gt_gmn)
perf <- performance(pred,"tpr","fpr")
auc_gmn <- performance(pred, "auc")@y.values[[1]]
par(new = TRUE)
pr.curve(scores.class0 = QCCD_MN_gmn$Eps, weights.class0 = pairs_gt_gmn$X1, curve = T);
plot(perf, col = "indianred3", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, type="l", lty=2)

auc_gnm_x = unlist(perf@x.values)
auc_gnm_y = unlist(perf@y.values)

# ROC
ROCR::performance(ROCR::prediction(QCCD_MN_gmn$Eps,pairs_gt_gmn$X1), "auc")@y.values[[1]]
PRROC::roc.curve(QCCD_MN_gmn$Eps, weights.class0 = pairs_gt_gmn$X1)

# PR
pr <- ROCR::performance(ROCR::prediction(QCCD_MN_gmn$Eps, pairs_gt_gmn$X1), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(QCCD_MN_gmn$Eps, weights.class0 = pairs_gt_gmn$X1)




ref.uv$V6_01[ref.uv$V6=="->"]=1
ref.uv$V6_01[ref.uv$V6=="<-"]=0
tuebcop <- read.table("../results/QCCD_Tueb.tab", as.is = TRUE, header = TRUE, sep = "\t")

pred <- prediction(tuebcop$Eps ,ref.uv$V6_01)
perf <- performance(pred,"tpr","fpr")
par(new = TRUE)
plot(perf,col = "red", lwd = 2, cex.lab = 1.5, cey.lab = 1.5, cex.axis = 4, cey.axis = 2, type = "l")
auc_tueb = performance(pred, "auc")@y.values[[1]]

pr.curve(scores.class0 = tuebcop$Eps , weights.class0 = ref.uv$V6_01, curve = T);
roc.curve(scores.class0 = tuebcop$Eps, weights.class0 = ref.uv$V6_01, curve = T);

auc_tueb_x = unlist(perf@x.values)
auc_tueb_y = unlist(perf@y.values)

# ROC
ROCR::performance(ROCR::prediction(tuebcop$Eps,ref.uv$V6_01), "auc")@y.values[[1]]
PRROC::roc.curve(tuebcop$Eps, weights.class0 = ref.uv$V6_01)

# PR
pr <- ROCR::performance(ROCR::prediction(tuebcop$Eps, ref.uv$V6_01), "prec", "rec")
integrate(splinefun(pr@x.values[[1]], pr@y.values[[1]]), 0, 1)
PRROC::pr.curve(tuebcop$Eps, weights.class0 = ref.uv$V6_01)



df1 <- data.frame(x = auc_an_x, y = auc_an_y, Type = as.factor("AN"))
df2 <- data.frame(x = auc_an_s_x, y = auc_an_s_y, Type = as.factor("AN-s"))
df3 <- data.frame(x = auc_h_x, y = auc_h_y, Type = as.factor("HN"))
df4 <- data.frame(x = auc_hn_x, y = auc_hn_y, Type = as.factor("HN-s"))
df5 <- data.frame(x = auc_sim_x, y = auc_sim_y, Type = as.factor("SIM"))
df6 <- data.frame(x = auc_simln_x, y = auc_simln_y, Type = as.factor("SIM-ln"))
df7 <- data.frame(x = auc_simg_x, y = auc_simg_y, Type = as.factor("SIM-G"))
df8 <- data.frame(x = auc_unm_x, y = auc_unm_y, Type = as.factor("MN-U"))
df9 <- data.frame(x = auc_gnm_x, y = auc_gnm_y, Type = as.factor("MN-G"))
df10 <- data.frame(x = auc_tueb_x, y = auc_tueb_y, Type = as.factor("Tueb"))

df.merged <- rbind(df1, df2, df3, df4,df5, df6, df7, df8, df9, df10)

ggplot(df.merged, aes(x, y, colour = Type, linetype = Type, shape = Type)) + geom_line() + geom_point() +
  geom_abline(intercept = 0, aes(colour = "rnd.class"),linetype="dotted")+
  xlab("False Positive Rate") + ylab("True Positive Rate") +
  theme_naglr() +
  theme(plot.margin=grid::unit(rep(1,4), "mm"),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x =element_text(size=8),axis.text.y =element_text(size=8),legend.title=element_blank())


