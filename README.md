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

## Data preperation
### Re-identification of phenotypes 
Phenotype of whitebritish was imported from UKB_Phenome and mapped with mapping.csv
Covariates(Sex, birthYear, PC1~4) was extracted as forcovar.txt

```
library(dplyr)
library(data.table)
pheno = fread('/media/leelabsg_storage01/DATA/UKBB/PheCode/UKB_Phenome/PEDMASTER_WhiteBritish_20180612.txt.gz')
pheno = cbind(pheno[,c(1:17)],pheno$X250.2)

map = fread("/media/leelabsg_storage01/DATA/UKBB/Mapping/mapping.csv", header=T)
join = left_join(pheno, map, by=c("IID"="Goncalo"))
join = join[,-c(1,2)]
colnames(join)[ncol(join)] = "IID"
join = cbind(FID = join$IID,IID = join$IID,join)
ncol(join)
join = join[,-19]
head(join)
out_path = "/home/leelabsg/seokho/SAIGE_T2D/forcovar.txt"
write.table(join, out_path, row.names=F, col.names=T, quote=F)
```
### get qualified GRM_markers from called genotype qc information
![image](https://user-images.githubusercontent.com/22064612/127957044-5e7537af-516b-428a-b07e-4ad9242a19c8.png)
![image](https://user-images.githubusercontent.com/22064612/127956977-e3d6a09d-170c-4c1f-8138-aab3e336df67.png)

93,511 markers were gathered to build GRM matrix

```
qc = fread('ukb_snp_qc.txt')
dim(qc)
names(qc)
qc2 = qc[abs(qc$PC1_loading)<0.003&abs(qc$PC2_loading)<0.003&abs(qc$PC3_loading)<0.003,]
dim(qc2)
write.table(qc2$rs_id,"/home/leelabsg/seokho/SAIGE_T2D/step0/GRM_markers.txt",row.names=F, col.names=F, quote=F)
```
Following the process resulted in 93,488 markers

### Extract Qualified Markers

```
for i in $(seq 1 22)
do
    nohup ~/plink \
    --bed "/media/leelabsg_storage01/DATA/UKBB/cal/ukb_cal_chr$i""_v2.bed" \
    --bim "/media/leelabsg_storage01/DATA/UKBB/cal/ukb_snp_chr$i""_v2.bim" \
    --fam /media/leelabsg_storage01/DATA/UKBB/cal/ukb45227_cal_chr20_v2_s488264.fam \
    --extract /home/leelabsg/seokho/SAIGE_T2D/step0/GRM_markers.txt \
    --make-bed \
    --out "/home/leelabsg/seokho/SAIGE_T2D/step0/extracted_chr"$i &
done
```

### Merge into genomwide plink binaries
```
nohup ~/plink \
--bed /home/leelabsg/seokho/SAIGE_T2D/step0/extracted_chr1.bed \
--bim /home/leelabsg/seokho/SAIGE_T2D/step0/extracted_chr1.bim \
--fam /home/leelabsg/seokho/SAIGE_T2D/step0/extracted_chr1.fam \
--merge-list /home/leelabsg/seokho/SAIGE_T2D/step0/allfiles.txt \
--make-bed --out /home/leelabsg/seokho/SAIGE_T2D/step0/merged &
```
## Step 1 fitNULLGLMM
```
nohup step1_fitNULLGLMM.R     \
	--plinkFile=/home/leelabsg/seokho/SAIGE_T2D/step0/merged \
  --phenoFile=/home/leelabsg/seokho/SAIGE_T2D/forcovar.txt \
  --phenoCol=V2 \
  --covarColList=Sex,birthYear,PC1,PC2,PC3,PC4 \
  --sexCol=Sex \
  --FemaleCode='2' \
  --MaleCode='1' \
  --sampleIDColinphenoFile=IID \
  --traitType=binary        \
  --outputPrefix=/home/leelabsg/seokho/SAIGE_T2D/output/T2D_binary_step1 \
  --nThreads=24 \
  --LOCO=FALSE \
  --IsOverwriteVarianceRatioFile=TRUE &

```
![image](https://user-images.githubusercontent.com/22064612/127957606-a7c0a120-d43c-4bf8-9542-971a35b6a097.png)

## Step 2 Run SPA tests
SPA tests were run on imputed chromosome 10 first. (GPU02)

### Current progress
![image](https://user-images.githubusercontent.com/22064612/127957803-f1f8ae0a-bc21-4e51-ad8e-51cf97ef1184.png)

### Comparing pvalues with pheweb results
From up to estimated markers, simple comparison was made with pheweb gwas summary
![image](https://user-images.githubusercontent.com/22064612/127957963-27fe990c-2754-42fe-8274-8e51c65aa524.png)

