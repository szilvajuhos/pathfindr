library(data.table)
library(stringr)

input_dir = "."
csv_dir = "."
chrsz = data.table(
  chr = paste0('chr',c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                       "11", "12", "13", "14", "15", "16", "17", "18", 
                       "19", "20", "21", "22", "X", "Y")), 
  label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
            "12", "13", "14", "15", "16", "17", "18", "19", "20", 
            "21", "22", "X", "Y"),
  starts = c(0, 253902921, 500669562, 703689723, 898025604, 1084280925, 
             1259870526, 1424144247, 1574191848, 1717288929, 1855896810, 
             1995944331, 2134165212, 2253485013, 2363996694, 2470839615, 
             2565902256, 2654091777, 2739309138, 2802869979, 2871941700, 
             2923598421, 2979326382, 3140008743),
  length = c(248902921, 241766641, 198020161, 189335881, 181255321, 
             170589601, 159273721, 145047601, 138097081, 133607881, 
             135047521, 133220881, 114319801, 105511681, 101842921, 
             90062641, 83189521, 80217361, 58560841, 64071721, 46656721,
             50727961, 155682361, 56827081)
  )


reference_genome = 'GRCh38'
baf = "the_BAF_file_bugger"
logr = "the_logR_file_bugger"

cat("Shaping tumor genes dataframe\n")
tumorgenes=fread('~/reports/reference_data/small.cancer_gene_census.csv',key='Gene Symbol')
tumorgenes$chr=paste0('chr',str_replace(string=tumorgenes$`Genome Location`,pattern = ':.*',replacement = ''))
temp=str_replace(string=tumorgenes$`Genome Location`,pattern = '.*:',replacement = '')
tumorgenes$start=as.numeric(str_replace(string=temp,pattern = '-.*',replacement = ''))
tumorgenes$end=as.numeric(str_replace(string=temp,pattern = '.*-',replacement = ''))
tumorgenes$cumstart=NA; tumorgenes$cumend=NA
for (i in 1:nrow(chrsz)) {
  ix <- tumorgenes$chr==chrsz$chr[i]
  tumorgenes$cumstart[ix] <- tumorgenes$start[ix]+chrsz$starts[i]
  tumorgenes$cumend[ix] <- tumorgenes$end[ix]+chrsz$starts[i]
}
browser()
cat("Merging locals\n")
#local_tumorgenes=fread('~/reports/reference_data/2018_gene_list_tere_ref.csv',key='Gene')[Gene!='',1:2]
local_tumorgenes = fread('~/reports/reference_data/2018_gene_list_tere_ref.csv',key='Gene')
# get rid of rows without a gene name
local_tumorgenes = subset.data.frame(local_tumorgenes,local_tumorgenes$Gene!='')
# selet only the gene names and the tier (a.k.a. score) 
local_tumorgenes = local_tumorgenes[,1:2]
alltumorgenes=unique(c(tumorgenes$`Gene Symbol`,local_tumorgenes$Gene))
cat("Creating tiers\n")
# select tier 1 genes from both
# then get a list of gene names as the union of the two sets
tumorgenes_tier1 = subset.data.frame(tumorgenes,tumorgenes$Tier==1)
local_tumorgenes_tier1 = subset.data.frame(local_tumorgenes,local_tumorgenes$`Tier 1 and 2 for pediatric cancers final`==1)
# select the first column (gene names)
tumorgenes_tier1_genes = tumorgenes_tier1[,1] 
local_tumorgenes_tier1_genes = local_tumorgenes_tier1[,1]
alltier1 = union(tumorgenes_tier1_genes,local_tumorgenes_tier1_genes)

# DRY but whatever, do the same for tier 2 TODO make it nice without repetition
tumorgenes_tier2 = subset.data.frame(tumorgenes,tumorgenes$Tier==2)
local_tumorgenes_tier2 = subset.data.frame(local_tumorgenes,local_tumorgenes$`Tier 1 and 2 for pediatric cancers final`==2)
alltier2 = union(tumorgenes_tier2[,1],local_tumorgenes_tier2[,1])

#alltier1=union(tumorgenes[Tier==1,tumorgenes$`Gene Symbol`],local_tumorgenes[`Tier 1 and 2 for pediatric cancers final`==1,local_tumorgenes$Gene])
#alltier2=union(tumorgenes[Tier==2,tumorgenes$`Gene Symbol`],local_tumorgenes[`Tier 1 and 2 for pediatric cancers final`==2,local_tumorgenes$Gene])
