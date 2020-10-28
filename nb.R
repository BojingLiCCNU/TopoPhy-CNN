library(SpiecEasi)
# getwd()
#path="E:\\augdata"
#setwd(path)
setwd("F:/augdata")
source('mysynfun.R')
source('myrmvnegbin.R')
source('myget_comm_params.R')
source('myfitdistr.R')
source('myrmvnorm.R')
source('myfixInf.R')
source('myrmvpois.R')


# rawdata = as.matrix(read.table("D:\\augdata\\leandata.txt"))
# 
# 
# multiple=2
# target=multiple*nrow(rawdata)
# 
# syndata=mysynfun(rawdata,mar=2,distr='negbin',  n=target)
# #syndata=synth_comm_from_counts(rawdata,mar=2,distr='negbin',  n=target)
# 
# 
# 
# savefile <- paste("D:\\augdata\\nblean",  round(multiple), ".txt", sep = "");
# 
# write.table(syndata, savefile, sep="\t", row.names = F, col.names = F)







total_train_sample = as.matrix(read.table("D:\\augdata\\rawdata\\t2d_raw_data.txt"))
total_train_label = as.matrix(read.table("D:\\augdata\\rawdata\\t2d_raw_label.txt"))

class0 <- total_train_sample[total_train_label[,1] == 1, ];
class1 <- total_train_sample[total_train_label[,1] == 0, ];


multiple <- c(1)
for (i in 1:length(multiple)){

  n_class0 <- round(nrow(class0)*multiple[i])
  n_class1 <- round(nrow(class1)*multiple[i])

  X0 <- mysynfun(class0, mar=2, distr='pois',  n=n_class0 )
  X1 <- mysynfun(class1, mar=2, distr='pois',  n=n_class1 )
  X <- rbind(X0, X1);
  label0 <- matrix(0, n_class0, 1);
  label1 <- matrix(1, n_class1, 1);
  labels <- rbind(label0, label1);

  sample_file <- paste("D:\\augdata\\results\\pois_t2d_data_",  multiple[i], ".txt", sep = "");
  label_file <- paste("D:\\augdata\\results\\pois_t2d_label_",  multiple[i], ".txt", sep = "");
  write.table(X, sample_file, sep="\t", row.names = F, col.names = F)
  write.table(labels, label_file, sep="\t", row.names = F, col.names = F)

}

# synth_comm_from_counts

