rm(list = ls())
source("./cd_utils.R")

  cnt = 0
  pair = rep("", length(uv))
  eps = rep(0, length(uv))
  cds = rep("--", length(uv))
  time = rep(0, length(uv))

  for(i in 1:length(uv)){
  t1 = Sys.time()
  t = readI(uv[i])[,1:2]
  res = QCCDwrap(t[,1], t[,2])
  t2 = Sys.time()
  elapsed = as.numeric(difftime(t2,t1), units="secs")
  pair[i] = uv[i]
  eps[i] = res$eps
  cds[i] = res$cd
  time[i] = elapsed
  }
corr = rep(0,length(uv))
corr[ref.uv$V6 == cds] = 1

print("Correct causal directions recovered:")
print(sum(corr))

print("Weighted accuracy:")
print(sum(QCCD_Tueb$Correct*meta$V6[uv])/sum(meta$V6[uv]))

resQCCD = data.frame(Correct = corr, Eps = eps, Cds = cds, T = time)
write.table(resQCCD, file = "../results/QCCD_Tueb.tab", row.names = F, quote = F, sep="\t")

