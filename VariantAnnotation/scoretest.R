library(VariantAnnotation)
vcf <- readVcf("mutect2test.vcf", "GRCh38")
names(info(vcf))
dim(info(vcf))
vcfInfo<-info(vcf)            # that is INFO data without header
vcfHeader<-info(header(vcf))  # that is INFO HEADER without data
vcfHeader@rownames
# naming the new column
newInfoCol<-DataFrame(Number=1,Type="Integer", Description="Bugger new data", row.names="RankScore")
newInfoCol@rownames
# merging the old info fields and the new column NAME - no data involved yet
vcfHeader <- rbind(vcfHeader,newInfoCol)
vcfHeader@rownames
vcfInfo$RankScore<-1:56
info(vcf)<-vcfInfo

info(header(vcf)) <- vcfHeader
rownames(info(header(vcf)))
writeVcf(vcf,"bugger.vcf")