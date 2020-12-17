## Library of custom functions. Last updated 12/17/2020

Basic function for loop template:
```
FUNCTION_NAME <- function(INPUT_FILES){
  INDEX_LIST = list()
  for(i in INDICES){
    OPPERATION
    INDEX_LIST[[i]] <- x
  }
  OUTPUT_FILE <- do.call(rbind, INDEX_LIST)
}
```


Replace nan in data frame with 0
```
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
```


Calculate the correlation between variables in a matrix
```
CalcGeneCor <- function(kal.data, h2.data, geneIDs){
  GeneList = list()
  for(i in geneIDs){
    x <- cor(kal.data[,i], h2.data[,i])
    GeneList[[i]] <- x
  }
  big_genes <- do.call(rbind, GeneList)
}
```


Calculate linear regression between variables in a matrix and record summary outputs
```
CalcPhenCor <- function(ind.data, survival, geneIDs){
  GeneList = list()
  for(i in geneIDs){
    x <- summary(lm(ind.data[,i] ~ Survival))
    y <- data.frame(Estimate = x$coefficients[2,1],
                    tvalue = x$coefficients[2,3],
                    rsquare = x$r.squared,
                    Fvalue = x$fstatistic[1],
                    Pvalue = x$coefficients[2,4])
    GeneList[[i]] <- y
  }
  big_genes <- do.call(rbind, GeneList)
}
```


Calculate founder effects for positions given a particular phenotype in the A DSPR panel.
This version includes a covariate (subpanel) which affects the Estimate calculation.
```
CalcFounderEffectsA_covariate <- function(MATRIX_DATA, CHR_POS){
  PosList = list()
  for(i in CHR_POS){
    i.position <- droplevels(subset(MATRIX_DATA, chr_pos == i))
    rr <- lm(survival ~ subpop, i.position)$resid
    sub.position <- i.position[,c(1,3:11,13,14)]
    sub.matrix <- as.matrix(sub.position[,-12])
    qtl.mod <- as.formula(rr ~  + sub.matrix[,3] + sub.matrix[,4] + sub.matrix[,5] + sub.matrix[,6] + sub.matrix[,7] + sub.matrix[,8] + sub.matrix[,9] + sub.matrix[,10] - 1)
    qtlmean <- lm(qtl.mod, data = sub.position)
    qtlmean.sum <- summary(qtlmean)
    ee <- summary(qtlmean)$coef[,1]
    lm.m <- lm(survival ~ 1, i.position)$coefficients[1]
    hap.data <- data.frame(chr_pos = rep(i, 8),
                           Estimate = ee + lm.m,
                           SE = qtlmean.sum$coefficients[1:8,"Std. Error"],
                           Founder = c("A1","A2","A3","A4","A5","A6","A7","AB8"))
    print(PosList[[i]] <- hap.data)
  }
  big.Effects <- do.call(rbind, PosList)
}
```


Calculate founder effects for positions given a particular phenotype in the B panel.
This version does not include a covariate.
```
CalcFounderEffectsB <- function(MATRIX_DATA, CHR_POS){
  PosList = list()
  for(i in CHR_POS){
    i.position <- droplevels(subset(MATRIX_DATA, chr_pos == i))
    sub.position <- i.position[,c(1,3:11,13)]
    sub.matrix <- as.matrix(sub.position)
    qtl.mod <- as.formula(paste(deparse(survival ~ 1), "+ sub.matrix[,3] + sub.matrix[,4] + sub.matrix[,5] + sub.matrix[,6] + sub.matrix[,7] + sub.matrix[,8] + sub.matrix[,9] + sub.matrix[,10] - 1", sep = ""))
    qtlmean <- lm(qtl.mod, data = sub.position)
    qtlmean.sum <- summary(qtlmean)
    hap.data <- data.frame(chr_pos = rep(i, 8),
                           Estimate = qtlmean.sum$coefficients[,"Estimate"],
                           SE = qtlmean.sum$coefficients[,"Std. Error"],
                           Founder = c("B1","B2","B3","B4","B5","B6","B7","AB8"))
    print(PosList[[i]] <- hap.data)
  }
  big.Effects <- do.call(rbind, PosList)
}
```


Calculate percent variance for multiple genes at multiple positions.
The first function is called by the second function to loop through many genes.
The strains must be sorted for all input files the same way.
```
CalcPerctVar <- function(phenotype, HAP_DF, CHR_POS, g){
  PosList = list()
  for(i in CHR_POS){
    i.position <- droplevels(subset(HapProb_MI_RILsubset_sorted, chr_pos == i))
    sub.position <- cbind(i.position[,c(3:11)], phenotype)
    sub.matrix <- as.matrix(sub.position)
    qtl.mod <- as.formula(paste(deparse(phenotype ~ 1), "+ sub.matrix[,2] + sub.matrix[,3] + sub.matrix[,4] + sub.matrix[,5] + sub.matrix[,6] + sub.matrix[,7] + sub.matrix[,8]", sep = ""))
    qtlmean <- lm(qtl.mod, data = sub.position)
    qtlmean.sum <- summary.aov(qtlmean)
    sumqtl <- as.data.frame(qtlmean.sum[[1]])
    ss <- sumqtl["sub.matrix[, 2]", "Sum Sq"] + sumqtl["sub.matrix[, 3]", "Sum Sq"] + sumqtl["sub.matrix[, 4]", "Sum Sq"] + sumqtl["sub.matrix[, 5]", "Sum Sq"] + sumqtl["sub.matrix[, 6]", "Sum Sq"] + sumqtl["sub.matrix[, 7]", "Sum Sq"] + sumqtl["sub.matrix[, 8]", "Sum Sq"]
    ndf <- 7
    ddf <- sumqtl["Residuals", "Df"]
    mss <- ss/ndf
    mse <- sumqtl["Residuals", "Mean Sq"]
    test.f <- mss/mse
    pctvar <- 100 * (test.f/(test.f + (ddf/ndf)))
    pctvar.data <- data.frame(PercentVariance = pctvar,
                              Gene = g)
    PosList[[i]] <- pctvar.data
  }
  big.PctVar <- do.call(rbind, PosList)
}


eQTL_PercentVariance <- function(GENE_LIST, EXP_DF, HAP_DF, CHR_POS){
  PosList = list()
  for(g in GENE_LIST){
    gene <- g
    print(gene)
    phenotype <- subset(EXP_DF, select = gene)
    colnames(phenotype)[1] <- "phenotype"
    chr_pos <- CHR_POS
    hap_df <- HAP_DF
    x <- CalcPerctVar(phenotype = phenotype, CHR_POS = chr_pos, HAP_DF = hap_df, g = gene)
    PosList[[g]] <- x
  }
  big.Effects <- do.call(rbind, PosList)
}
```


Function for combining sequencing files (shell)
```
for file in *L00{1,5}*.gz
do
  sample=${file%_S*_*.fastq.gz}
  zcat ${sample}_*.fastq.gz >> ../combined/${sample}.fastq.gz
done
```


To subsample RNAseq files build this .py script:
```
import sys, random, itertools
import gzip
import HTSeq

fraction = float( sys.argv[1] )
fileR1 = gzip.open( sys.argv[2], "rb" )
fileR2 = gzip.open( sys.argv[3], "rb" )
in1 = iter( HTSeq.FastqReader( fileR1 ) )
in2 = iter( HTSeq.FastqReader( fileR2 ) )
out1 = gzip.open( sys.argv[4], "wb" )
out2 = gzip.open( sys.argv[5], "wb" )

for read1, read2 in itertools.izip( in1, in2 ):
   if random.random() < fraction:
      read1.write_to_fastq_file( out1 )
      read2.write_to_fastq_file( out2 )

out1.close()
out2.close()



# To Implement, use:
python subsamplePE.py 0.25 ../data/filtered/Cu_22325_S20.R1.filt.fastq.gz ../data/filtered/Cu_22325_S20.R2.filt.fastq.gz ../data/subset8mil/Cu_22325_S20.sub.R1.filt.fastq.gz ../data/subset8mil/Cu_22325_S20.sub.R2.filt.fastq.gz


Files are called subsamplePE.py and finishingSubsetting.sh
```


Regress out PCs for eQTL file prep:
```
ExtractResiduals <- function(QNORM, TARGETPCs){
  index_list = list()
  for(i in 1:ncol(QNORM)) {
    # EXTRACT vector of expression values for single gene
    single.gene.exp <- QNORM[,i]
    
    # GENERATE data frame of the expression vector and the target PCs
    data.tmp <- data.frame(geneexp = single.gene.exp, TARGETPCs)
    
    # USE lm to fit PCs
    # GET residuals
    # ADD to output 'exp.resids' object and
    # ADD column name (as gene symbol)
    PCfit <- lm(geneexp~.,data=data.tmp)
    PCfit.resids <- as.data.frame(PCfit$residuals)
    colnames(PCfit.resids) <- colnames(QNORM)[i]
    index_list[[i]] <- PCfit.resids
  }
  OUTPUT_FILE <- do.call(cbind, index_list)
}
```
