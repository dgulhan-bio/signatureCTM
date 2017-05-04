signatures<-read.csv('nmfresult.dat',header = FALSE)
nsig<-dim(signatures)[[2]]

vocab_list <- read.table('vocabSorted.txt')
vocab_mtrx<-array(unlist(vocab_list), dim = c(nrow(vocab_list), ncol(vocab_list), length(vocab_list)))


for(i in 1:nsig){
  print(i)
  png(sprintf('barsnmf_num%d_totalnum%d.png',i,nsig), width=6000, height=1500, res=300)
  barplot(signatures[,i], col=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6), las=2, ylim=c(0, 0.25),names.arg = vocab_mtrx)
  dev.off()
}
