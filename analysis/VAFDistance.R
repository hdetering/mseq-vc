df_mut_sample$chrom <- sapply(df_mut_sample$id_mut, function(x)  df_mut$chrom[match(x, df_mut$id_mut)])
df_mut_sample$pos <- sapply(df_mut_sample$id_mut, function(x)  df_mut$pos[match(x, df_mut$id_mut)])
df_mut_sample$ref <- sapply(df_mut_sample$id_mut, function(x)  df_mut$ref[match(x, df_mut$id_mut)])
df_mut_sample$alt <- sapply(df_mut_sample$id_mut, function(x)  df_mut$alt[match(x, df_mut$id_mut)])
# create a new column mutinfo with the four columns collapsed together
cols <- c( 'chrom' , 'pos' , 'ref', 'alt' )
df_mut_sample$mut_info <- apply( df_mut_sample[ , cols ] , 1 , paste , collapse = "_" )
cols <- c( 'chrom' , 'pos')
df_mut_sample$chrom_pos <- apply( df_mut_sample[ , cols ] , 1 , paste , collapse = "_" )


getVAFs <- function(SOFTWARE, TAG, REPLICATE_ID){
  print(SOFTWARE)
  print(REPLICATE_ID)
  ## OBTAIN EMPIRICAL VAFs
  if (SOFTWARE=="SNooPer") {
    SOFTWARE="SNooPerCorrected"
  } else if (SOFTWARE=="Strelka2"){
    SOFTWARE="Strelka2Corrected"}
  vcf <- read.vcfR(paste("/Users/tama/Downloads/VAFs/", REPLICATE_ID,".",SOFTWARE,".PASS.vcf", sep=""), verbose = F)
  if (is.na(TAG)){
    return(NA)
  }else if (TAG == "BCOUNT-DP"){
    AlternativeAlleles <- sub(",.*","",getALT(vcf))
    # BCOUNT: Occurrence count for each base at this site (A,C,G,T)
    BCOUNT <- extract.gt(vcf,element = "BCOUNT", as.numeric = FALSE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
    # DP: Total read depth
    DP <- extract.gt(vcf,element = "DP", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
    
    # letter: 1>A 2>C 3>G 4>T            
    getAlleles <- function(x, letter) {
      as.numeric(sub("^([0-9]+),([0-9]+),([0-9]+),([0-9]+)", paste("\\",letter,sep=""), x))
    }
    A_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=1)
    C_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=2)
    G_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=3)
    T_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=4)
    
    
    GetAlternativeCounts <- function(pos){
      ALT_ALL <- AlternativeAlleles[pos]
      alternative_counts_matrix <- get(paste(ALT_ALL, "_counts", sep=""))
      freq <- alternative_counts_matrix[pos,]/DP[pos,]
      return(freq)}
    a<-lapply(1:nrow(DP), FUN = GetAlternativeCounts)
    freq_matrix <- as.data.frame(do.call(rbind, a))
    freq_matrix
  } else if  (TAG == "FREQ" ) {
    freq_matrix <- extract.gt(
      vcf,
      element = TAG,
      mask = FALSE,
      as.numeric = FALSE,
      return.alleles = FALSE,
      IDtoRowNames = TRUE,
      extract = TRUE,
      convertNA = TRUE)
    #remove % symbol and divide percentage to get proportion
    tmp = freq_matrix %>% data.frame(stringsAsFactors = FALSE) %>% mutate(mut = rownames(freq_matrix)) %>% gather(key = "sample", value = "vaf", 1:ncol(freq_matrix)) %>%
      mutate(vaf = as.numeric(sub("%","", vaf))/100)      %>%
      spread(key = sample, value = vaf)
    freq_matrix = as.matrix.data.frame(tmp %>% select(-c(mut)))
    rownames(freq_matrix) = tmp$mut
  } else if (TAG=="VR-DP") {
    vr_matrix <- extract.gt(
      vcf,
      element = "VR",
      mask = FALSE,
      as.numeric = TRUE,
      return.alleles = FALSE,
      IDtoRowNames = TRUE,
      extract = TRUE,
      convertNA = TRUE)
    dp_matrix <- extract.gt(
      vcf,
      element = "DP",
      mask = FALSE,
      as.numeric = TRUE,
      return.alleles = FALSE,
      IDtoRowNames = TRUE,
      extract = TRUE,
      convertNA = TRUE)
    freq_matrix = as.matrix(vr_matrix / dp_matrix)
    
  } else if (TAG=="U") {
    AlternativeAlleles <- sub(",.*","",getALT(vcf))
    # DP: Total read depth
    DP <- extract.gt(vcf,element = "DP", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)  
    AU <- extract.gt(vcf,element = "AU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
    CU <- extract.gt(vcf,element = "CU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
    GU <- extract.gt(vcf,element = "GU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
    TU <- extract.gt(vcf,element = "TU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
    GetAlternativeCounts <- function(pos){
      ALT_ALL <- AlternativeAlleles[pos]
      alternative_counts_matrix <- get(paste(ALT_ALL, "U", sep=""))
      freq <- alternative_counts_matrix[pos,]/DP[pos,]
      return(freq)}
    a<-lapply(1:nrow(DP), FUN = GetAlternativeCounts)
    freq_matrix <- as.data.frame(do.call(rbind, a))
    
  } else if (SOFTWARE=="Mutect2_multi_F") {
    
    # REMOVE FREQ IN POSITIONS WITHOUT ALTERNATIVE READS
    AD_matrix <- extract.gt(
      vcf,
      element = 'AD',
      mask = FALSE,
      as.numeric = FALSE,
      return.alleles = FALSE,
      IDtoRowNames = TRUE,
      extract = TRUE,
      convertNA = TRUE
    )
    
    getAltAlleleFromAD <- function(x) {
      ifelse(as.numeric(sub("^([0-9]+),([0-9]+)", paste("\\2"), x))==0, FALSE, TRUE)
    }
    Pres_Abs <- apply(AD_matrix, 2, FUN = getAltAlleleFromAD)
    ToKeepRemove<-reshape2::melt(Pres_Abs[,-1])
    
    freq_matrix <- extract.gt(
      vcf,
      element = "AF",
      mask = FALSE,
      as.numeric = FALSE,
      return.alleles = FALSE,
      IDtoRowNames = TRUE,
      extract = TRUE,
      convertNA = TRUE)
    
  } else {
    freq_matrix <- extract.gt(
      vcf,
      element = TAG,
      mask = FALSE,
      as.numeric = FALSE,
      return.alleles = FALSE,
      IDtoRowNames = TRUE,
      extract = TRUE,
      convertNA = TRUE
    )        
  }
  if ( TAG =="AD" ) {
    freq_matrix <- AD_frequency(freq_matrix, delim = ",", allele = 2L, sum_type = 1L, decreasing = 1L)
  }
  freq_data <- data.frame(freq_matrix)
  freq_data$mut_info <- apply(getFIX(vcf)[,c("CHROM","POS","REF","ALT")],1,paste,collapse="_"    )
  freq_data$chrom_pos <- apply(getFIX(vcf)[,c("CHROM","POS")],1,paste,collapse="_")
  colnames(freq_data) <- sub("\\.S.*","",colnames(freq_data))
  emp_data<- melt(freq_data, id.vars = c("mut_info","chrom_pos"))
  colnames(emp_data) <- c("mut_info","chrom_pos","id_sample","AF")
  emp_data <- emp_data[which(emp_data$id_sample!="Chines"),]
  emp_data <- emp_data[which(emp_data$id_sample!="N"),]
  
  if (SOFTWARE == "Mutect2_multi_F"){
    emp_data<-emp_data[which(ToKeepRemove$value==TRUE),]
  }
  
  if (SOFTWARE == "MultiSNV"){
    emp_data<-emp_data[which(!is.na(emp_data$AF)),]
  }
  
  # remove FP not in chromosome 21 (In bcftools there are variants in chromosome 1, don't know why)
  emp_data <- emp_data[grepl( "^21" , emp_data$chrom_pos) ,  ]

  
  # OBTAIN SIMULATED VAFS   
  # Select only the Replicate correspondent with the empirical data
  VcfEquivalentTruereplicate <-     df_rep[which(df_rep$name_rep==REPLICATE_ID),]
  id_equiv <- as.character(VcfEquivalentTruereplicate["id_rep"])
  df_mut_sample_equivalentreplicate <-     df_mut_sample[which(df_mut_sample$id_rep==id_equiv),]
  
  emp_data_small <- emp_data[,c("chrom_pos","id_sample","AF")]
  
  emp_data_small_noNA <- emp_data_small[!is.na(emp_data_small$AF),]
  df_mut_sample_equivalentreplicate_small <-     df_mut_sample_equivalentreplicate[,c("vaf_exp","chrom_pos", "id_sample")]
  df_mut_sample_selected <- join(emp_data_small_noNA, df_mut_sample_equivalentreplicate_small,     by = c("chrom_pos", "id_sample"), type="left", match="first")

  
  # substitute NA values by zeros
  df_mut_sample_selected[is.na(df_mut_sample_selected$vaf_exp),"vaf_exp"] <- 0
  df_mut_sample_selected <- df_mut_sample_selected[!is.na(df_mut_sample_selected$AF),]
  
  matrixfordist <- rbind(df_mut_sample_selected$AF, df_mut_sample_selected$vaf_exp)
  distance <- dist(matrixfordist, method = "euclidean") 
  return(distance)
}

# Create a table with all the combinations of factors we want to get the distance of allele frequencies
df_caller$name_caller
TAGS <-  c(
  "AD", # HaplotypeCaller
  "PM", # CaVEMan
  "FA", # MuTect1
  "AF", # Mutect2_single
  "VR-DP", # SNooPer
  "BCOUNT-DP", # SomaticSniper
  "U", # Strelka2
  "AF", # VarDict
  NA, # Shimmer
  "BCOUNT-DP", # MultiSNV
  NA, # SNV-PPILP
  NA, # MuClone
  "AF", # NeuSomatic
  "FREQ", # VarScan
  "AD", # Bcftools
  "AF") # Mutect2_multi_F

caller_AFtags <- cbind.data.frame(df_caller$name_caller, TAGS)
colnames(caller_AFtags) <- c("name_caller", "AF_TAG")


Replicates_row<- as.data.frame(t(df_rep$name_rep))
library(dplyr)
replicates_rep <- Replicates_row %>% slice(rep(1:n(), each = nrow(caller_AFtags)))
caller_AFtags_repli <- cbind.data.frame(caller_AFtags, replicates_rep)
caller_AFtags_repli_info <- reshape2::melt(caller_AFtags_repli, id.vars=c("name_caller", "AF_TAG"))
caller_AFtags_repli_info <- caller_AFtags_repli_info[,-3]
colnames(caller_AFtags_repli_info) <- c("name_caller", "AF_TAG", "Replicate")
Strategy <- sub("\\..*","",caller_AFtags_repli_info$Replicate)

# Calculate the distance for each row
Distance <- apply(caller_AFtags_repli_info, 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate']))

AF_Results <- cbind.data.frame(caller_AFtags_repli_info, Strategy, Distance)
saveRDS(AF_Results, file = "/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/RRSV/Results/RRSV.AFdistances.rds")

AF_Results$type <- sapply(AF_Results$name_caller, function(x) df_caller$class[match(x, df_caller$name_caller)])

# remove MultiSNV distance (Inf result)
AF_Results_final <- subset(AF_Results, name_caller != "MultiSNV")

AF_Results_final$name_caller<- sub("_","\n",AF_Results_final$name_caller)

AF_Results_final$name_caller_sorted <- with(AF_Results_final, reorder(name_caller , Distance, median , na.rm=T))

# Select the same two colors than for the rest of plots even with a levels missing
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 3
cols = gg_color_hue(n)
mycols=c(cols[1],cols[3])

AF_Results_final$type <- factor(AF_Results_final$type, levels = c("marginal","joint"))

ggplot(na.omit(AF_Results_final)) +
  geom_boxplot(aes(name_caller, Distance, fill=type)) +
  facet_grid(. ~ name_caller_sorted, scales = "free") +
  labs(x="caller", y="Euclidean distance between\nempirical and true alternative allele fractions\nof both true and false positives") +
  scale_fill_manual(values = mycols) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "bottom")

```