<!---
seokho92/seokho92 is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
# SAIGE T2D Anaylsis
## Conda Environment Setup
```
conda create -n saige -c conda-forge -c bioconda "r-base>=4.0" r-saige
conda activate saige
```

## R code for data preperation

```
library(dplyr)
library(data.table)
pheno = fread('/media/leelabsg_storage01/DATA/UKBB/PheCode/UKB_Phenome/PEDMASTER_WhiteBritish_20180612.txt.gz')
pheno = cbind(pheno[,c(1:17)],pheno$X250.2)

map = fread("/media/leelabsg_storage01/DATA/UKBB/Mapping/mapping.csv", header=T)
sample = fread("/home/lee7801/DATA/UKBB/imp/ukb45227_imp_chr1_v3_s487296.sample")

join = left_join(pheno, map, by=c("IID"="Goncalo"))
join2 = join[,c(18,19)]
join2 = inner_join(sample, join2, by=c("ID_1"="Shawn"))
colnames(join2)[5] = "pheno"
join2 = rbind(sample[1,],join2,fill=T)
join2$pheno[1] = "B"
join2$pheno[is.na(join2$pheno)] = -9
head(join2)
out_path = "/home/leelabsg/seokho/SAIGE_T2D/T2D_sample.txt"
write.table(join2, out_path, row.names=F, col.names=T, quote=F)

#join = left_join(pheno, map, by=c("IID"="Goncalo"))
join = join[,-c(1,2)]
colnames(join)[ncol(join)] = "IID"
join = cbind(FID = join$IID,IID = join$IID,join)
ncol(join)
join = join[,-19]
head(join)
out_path = "/home/leelabsg/seokho/SAIGE_T2D/forcovar.txt"
write.table(join, out_path, row.names=F, col.names=T, quote=F)

```

