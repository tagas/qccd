rm(list = ls())

library(dplyr)
library(readr)
source("./cd_utils.R")

n_pairs = 100
n_size = 1000

set.seed(3)

path = "../data/ANHNMN_pairs/"
ext = "MN-U/"

# load ground truth for pairs
pairs_gt <- read_csv(paste0(path, ext, "pairs_gt.txt"), col_names = FALSE)
pairs_gt = pairs_gt$X1



runANM = function(method){
  
  method_to_run = switch(method,
                        "QCCD" = QCCD)  

  pair_idx = rep("", n_pairs)
  eps = rep(0, n_pairs)
  cds = rep("--", n_pairs)
  time = rep(0, n_pairs)
  eps = rep(0, n_pairs)
  Correct = rep(0, n_pairs)

      n_corr = 0
      results = matrix(0,nrow = n_pairs, ncol = 2)
      
      for(i in 1: 100){
        print(paste("Pair #",i))
        t1 = Sys.time()
        pair =  read.table(paste0(path, ext, "pair_",i,".txt"), as.is = TRUE, header = TRUE, sep = ",", row.names = 1)
        
        plot(pair)
        rm(res)
        res = method_to_run(pair)
        t2 = Sys.time()
        elapsed = as.numeric(difftime(t2, t1), units = "secs")
        pair_idx[i] = i
        eps[i] = res$epsilon
        cds[i] = res$cd
        time[i] = elapsed
        
        # dummy prints for debug
        print(res$cd)
        print(pairs_gt[i])
        print("-----")
        results[i,1] = res$cd
        
        if(!is.na(res$cd)){
          correct = (res$cd == pairs_gt[i] )
          n_corr = n_corr + correct
        }
      }
      ncorrect = sum(results[,1]  == pairs_gt )
      print(paste("correct =", n_corr, "mean =", ncorrect))

      corr = rep(0, n_pairs)
      corr[pairs_gt == cds] = 1
      resDF = data.frame(Correct = corr, Eps = eps, Cds = cds,T = time)
      write.table(resDF,file = paste0("../results/",method,"_",substr(ext,1,nchar(ext)-1),".tab"),row.names = F, quote = F, sep = "\t")
  acc = ncorrect/100
  acc
}

methods <- list("QCCD") 
resultsAll <- lapply(methods, function (m) runANM(m))
