data_dir="/Users/tama/Google\ Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/RRSV/Results/"
data_root="/Users/tama/Downloads/VAFs/"


####################################
#       Obtain reported VAFs      #
####################################

getVAFs <- function(SOFTWARE, TAG, REPLICATE_ID){
      print(SOFTWARE)
      print(REPLICATE_ID)
      ## OBTAIN EMPIRICAL VAFs
      if (SOFTWARE=="SNooPer") {
        SOFTWARE="SNooPerCorrected"
        } else if (SOFTWARE=="Strelka2"){
        SOFTWARE="Strelka2Corrected"}
      fn_vcf <- paste(data_root, REPLICATE_ID,".",SOFTWARE,".PASS.vcf", sep="")
      vcf <- read.vcfR(fn_vcf, verbose = F)
      if (is.na(TAG)){
        return(NULL)
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
      } else if (TAG == "AD") {
    freq_matrix <- extract.gt(vcf,element = TAG,
                              mask = FALSE, as.numeric = FALSE,
                              return.alleles = FALSE, IDtoRowNames = TRUE,
                              extract = TRUE, convertNA = TRUE) 
    freq_matrix <- AD_frequency(freq_matrix, delim = ",", allele = 2L, sum_type = 1L, decreasing = 1L)
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

    freq_data <- data.frame(freq_matrix)
    freq_data$mut_info <- apply(getFIX(vcf)[,c("CHROM","POS","REF","ALT")],1,paste,collapse="_"    )
    freq_data$chrom_pos <- apply(getFIX(vcf)[,c("CHROM","POS")],1,paste,collapse="_")
    colnames(freq_data) <- gsub("\\.variant","",sub("\\.S.*","",colnames(freq_data)))
    if (SOFTWARE!="Bcftools"){
    freq_data <- freq_data %>%
      dplyr::select(-one_of("RN","N","ChineseSon.H","NORMAL","NORMAL2","NORMAL3","NORMAL4","NORMAL5"))
    }
  #--------------------------------------------#    
  # Filter VAFs depending on the caller        #
  #--------------------------------------------#
    
    if (SOFTWARE=="Bcftools") {
      vcf_table <- read.vcfR(fn_vcf, verbose = F) %>% extract_gt_tidy()
      vcf_table <- reshape2::melt(read.vcfR(fn_vcf, verbose = F) %>% 
                                    extract.gt(element = 'GT'))
      colnames(vcf_table) <- c("chrom_pos","Indiv","gt_GT")
    } else {
          vcf_table <- read.vcfR(fn_vcf, verbose = F) %>% vcfR2tidy()
          FILTER  <-vcf_table$fix$FILTER
          freq_data <- cbind.data.frame(freq_data, FILTER)
    }
    # --- MuTect1 Mutect2 Strelka1 Strelka2 NeuSomatic SNooPer'------------------
    if ( SOFTWARE %in% c('MuTect1', 'Mutect2', 'Strelka1', 'Strelka2', 'NeuSomatic','SNooPer') ) {
      freq_data <- freq_data %>%
      dplyr::filter( FILTER == 'PASS' ) %>%
      select (-FILTER)
    # --- VarDict ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarDict' ) {
      STATUS <- vcf_table$fix$STATUS
      freq_data <- cbind.data.frame(freq_data, STATUS) %>%
      dplyr::filter( FILTER == 'PASS' & grepl('Somatic$', STATUS) ) %>%
      select (-FILTER,-STATUS)
    # --- VarScan ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarScan' ) {
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_SS,Indiv) %>%
        dplyr::mutate(ChromKey=21) %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        inner_join(freq_data, by=c('chrom_pos')) %>%
        dplyr::filter( FILTER == 'PASS') %>%
        dplyr::filter( gt_SS == 2 ) %>%
        select (-FILTER,-gt_SS,-ChromKey,-POS,-Indiv) %>%
        select (ChineseSon.T1,ChineseSon.T2,ChineseSon.T3,ChineseSon.T4,ChineseSon.T5,mut_info,chrom_pos)
    # --- MultiSNV --------------------------------------------------------------
    } else if ( SOFTWARE == 'MultiSNV' ) {
        freq_data <- freq_data %>%
        reshape2::melt(id.vars=c("chrom_pos","mut_info","FILTER"))
        colnames(freq_data)[4] <- "Indiv"
      SS <- extract.gt(vcf,element = "SS", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)[,-1]
      SS <- cbind(rownames(SS), data.frame(SS, row.names=NULL))
      colnames(SS)[1] <- "chrom_pos"
      freq_data <- SS %>%
        reshape2::melt(variable.name="Indiv",value.name="gt_SS") %>%
        dplyr::inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter( FILTER == 'PASS' ) %>%
        dplyr::filter( gt_SS == 2 ) %>%
        dplyr::select (-FILTER,-gt_SS)  %>%
        spread(key = Indiv,value = value, convert = FALSE) %>%
        dplyr::select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    # --- HaplotypeCaller  --------------------------------------------
    } else if (SOFTWARE %in% c('HaplotypeCaller')){
      INDEL<-is.indel(vcf)
      freq_data <-cbind.data.frame(freq_data,INDEL) %>%
        reshape2::melt(variable.name="Indiv")
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,Indiv) %>%
        dplyr::mutate(ChromKey=21) %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        dplyr::filter(Indiv!="ChineseSon.H") %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>%
        inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter(gt_GT!='0/0') %>%
        dplyr::filter( INDEL == FALSE ) %>%
        select (-FILTER,-gt_GT,-ChromKey,-POS,-INDEL) %>%
        spread(key = Indiv,value = value, convert = FALSE) %>%
        select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    # --- Bcftools -------------------------------------------------------------
    } else if (SOFTWARE %in% c('Bcftools')) {
      freq_data <- reshape2::melt(freq_data,id.vars=c("mut_info","chrom_pos"),variable.name="Indiv")
      freq_data <- vcf_table %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>%
        inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        #mutate(value = ifelse (gt_GT!='0/0',NA,value)) %>%
        dplyr::filter( gt_GT != '0/0' ) %>%
        dplyr::filter( grepl('^21', chrom_pos)) %>%
        select(-gt_GT) %>%
        spread(key = Indiv,value=value) %>%
        select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    } else if ( SOFTWARE == 'Mutect2_multi_F' ) {
      freq_data <- reshape2::melt(freq_data,id.vars=c("mut_info","chrom_pos","FILTER"),variable.name="Indiv")
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,gt_AD,Indiv) %>%
        dplyr::mutate(ChromKey=21) %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        dplyr::filter(Indiv!="ChineseSon.H") %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>%
        inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter( FILTER == 'PASS' ) %>%
        dplyr::filter( gt_GT == '0/1' ) %>%
        rowwise() %>% 
        mutate( rc_alt = as.numeric(str_split(gt_AD, ',', simplify=T)[2]) ) %>%
        dplyr::filter( rc_alt > 0 ) %>%
        select (-FILTER,-gt_GT,-gt_AD,-ChromKey,-POS,-rc_alt) %>%
        spread(key = Indiv,value=value) %>%
        select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    } else {
      freq_data <- freq_data %>% select(-FILTER)
    }
    output <- cbind.data.frame(freq_data, rep(SOFTWARE,nrow(freq_data)), rep(REPLICATE_ID,nrow(freq_data)))
    colnames(output) <- c("T1","T2","T3","T4","T5","mut_info","chrom_pos","software","replicate")
    return(output)
}

# Create a table with all the combinations of factors we want to get the distance of allele frequencies
CallersInfo <- readRDS(paste(data_dir,"RRSV.callers.rds",sep=""))
CallersInfo$name_caller
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
caller_AFtags <- cbind.data.frame(CallersInfo$name_caller, TAGS)
colnames(caller_AFtags) <- c("name_caller", "AF_TAG")
Replicates_row<- as.data.frame(t(Replicates_Info$name_rep))
library(dplyr)
replicates_rep <- Replicates_row %>% slice(rep(1:n(), each = nrow(caller_AFtags)))
caller_AFtags_repli <- cbind.data.frame(caller_AFtags, replicates_rep)
caller_AFtags_repli_info <- reshape2::melt(caller_AFtags_repli, id.vars=c("name_caller", "AF_TAG"))
caller_AFtags_repli_info <- caller_AFtags_repli_info[,-3]
colnames(caller_AFtags_repli_info) <- c("name_caller", "AF_TAG", "Replicate")
Strategy <- sub("\\..*","",caller_AFtags_repli_info$Replicate)


# Obtain the AFs
AFs <- do.call("rbind", apply(caller_AFtags_repli_info, 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
saveRDS(AFs, file = paste(data_dir,"RRSV.AFs.rds",sep=""))


####################################
#      Obtain simulated VAFs       #
####################################

True_VAFs_muts<- readRDS(paste(data_dir,"RRSV.mutations_samples.rds",sep=""))
mutations <- readRDS(paste(data_dir,"RRSV.mutations.rds",sep=""))
Replicates_Info <- readRDS(paste(data_dir,"RRSV.replicates.rds",sep=""))


df_vaf_exp <- True_VAFs_muts %>%
  inner_join( mutations, by = c('id_rep', 'id_mut') ) %>%
  inner_join( Replicates_Info, by = c('id_rep') ) %>%
  unite( 'chrom_pos', chrom, pos) %>%
  mutate( vaf_exp = as.numeric(vaf_exp) ) %>%
  select( name_rep, id_sample, chrom_pos, vaf_exp)



#####################################
# Merge reported and simulated VAFs #
#####################################

AFs <-readRDS(paste(data_dir,"RRSV.AFs.rds",sep=""))
AFs$name_rep <- as.character( AFs$replicate )
df_vaf_obs <- AFs %>% pivot_longer( cols = paste0('T', 1:5), names_to = 'id_sample', values_to = 'vaf' ) %>%
  mutate( vaf_obs = as.numeric(vaf) ) %>%
  select( -mut_info, -replicate, -vaf )

df_vaf <- df_vaf_obs %>% 
  dplyr::left_join( df_vaf_exp, by = c('name_rep', 'id_sample', 'chrom_pos') )


getDistances <-  function(SOFTWARE, REPLICATE_ID){
  
  print (SOFTWARE)
  print (REPLICATE_ID)
  True_VAFs_muts_selected <- df_vaf %>%
    dplyr::filter(grepl(REPLICATE_ID,name_rep)) %>%
    dplyr::filter(grepl(SOFTWARE,software))
  
  ##############################################
  # Substitute NAs by 0s at FPs and remove TNs #
  ##############################################
  True_VAFs_muts_selected[is.na(True_VAFs_muts_selected$vaf_exp),"vaf_exp"] <- 0
  True_VAFs_muts_selected <- True_VAFs_muts_selected[!is.na(True_VAFs_muts_selected$vaf_obs),]
  
  
  matrixfordist <- rbind(True_VAFs_muts_selected$vaf_obs, True_VAFs_muts_selected$vaf_exp)
  distance <- dist(matrixfordist, method = "euclidean") 
  return(distance)
  
}


# Calculate distances
Distance <- apply(caller_AFtags_repli_info, 1, function(x) getDistances(SOFTWARE = x['name_caller'], REPLICATE_ID = x['Replicate']))

AF_Distance_Results <- cbind.data.frame(caller_AFtags_repli_info, Strategy, Distance)
saveRDS(AF_Distance_Results, file = paste(data_dir,"RRSV.AFdistances.rds",sep=""))

