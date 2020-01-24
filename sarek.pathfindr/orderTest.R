library(data.table)
# tumor log ratio
cat("Testing temp order\n")
temp=fread(file = './Sarek/test/Ascat/P2233_101T.LogR' )[,-1]
colnames(temp)[3]='LogR'
temp$sample="P2233_101T" #samplename[s]
temp=temp[order(Chr,Position),c(4,1:3)]
temp$Ratio=2^temp$LogR
temp$smoothed=runmed(x = temp$Ratio,k = 99)