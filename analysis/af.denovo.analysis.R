data_dir="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/de-novo/data/"
data_root="/Users/tama/Downloads/VCFs_SRSV/"


########################
# OBTAIN REPORTED VAFs #
########################
# Create a function to obtain AFs per caller per replicate

getVAFs <- function(SOFTWARE, TAG, REPLICATE_ID){
      print(SOFTWARE)
      print(REPLICATE_ID)
  # Obtain the vcf file name    
      fn_vcf <- file.path( data_root, REPLICATE_ID, 
                                switch(SOFTWARE,
                                       HaplotypeCaller = 'haplotypecaller.filt.vcf',
                                       Bcftools = 'bcftools.vcf',
                                       CaVEMan = 'caveman.vcf',
                                       MuTect1 = 'mutect1.vcf',
                                       Mutect2_single = 'mutect2.vcf',
                                       NeuSomatic = 'neusomatic.vcf',
                                       SNooPer = 'snooper.corrected.vcf',
                                       SomaticSniper = 'somaticsniper.vcf',
                                       Strelka2 = 'strelka2.corrected.vcf',
                                       Shimmer = 'shimmer.vcf',
                                       VarDict = 'vardict.vcf',
                                       VarScan = 'varscan.vcf',
                                       MultiSNV = 'multisnv.vcf',
                                       `SNV-PPILP` = 'snv-ppilp.vcf',
                                       MuClone = 'muclone.vcf',
                                       `Mutect2_multi_F` = 'mutect2m.vcf')
                           )
        if (SOFTWARE=="Mutect2_multi_F") {
        SOFTWARE="Mutect2_mseq"}
  # Read VCF
      vcf <- read.vcfR(fn_vcf, verbose = F)
  #--------------------------------------------#    
  # Obtain empirical VAFs depending on the TAG #
  #--------------------------------------------#    
    # --- TAG = NA --------------------------------------------------------------
      if (is.na(TAG)){
        return(NA)
    # --- TAG = BCOUNT-DP --------------------------------------------------------------
      } else if (TAG == "BCOUNT-DP"){
        AlternativeAlleles <- sub(",.*","",getALT(vcf))
        BCOUNT <- extract.gt(vcf,element = "BCOUNT", as.numeric = FALSE, return.alleles = FALSE,   IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE) # BCOUNT: Occurrence count for each base at this site (A,C,G,T)
        DP <- extract.gt(vcf,element = "DP", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames =   TRUE,extract = TRUE,convertNA = TRUE)  # DP: Total read depth
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
    # --- TAG = FREQ --------------------------------------------------------------
      } else if  (TAG == "FREQ" ) {
        freq_matrix <- extract.gt(vcf,element = TAG,
                                  mask = FALSE,
                                  as.numeric = FALSE,
                                  return.alleles = FALSE,
                                  IDtoRowNames = TRUE,
                                  extract = TRUE,
                                  convertNA = TRUE)
        #remove % symbol and divide percentage to get proportion
        tmp = freq_matrix %>% data.frame(stringsAsFactors = FALSE) %>% 
          mutate(mut = rownames(freq_matrix)) %>% 
          gather(key = "sample", value = "vaf", 1:ncol(freq_matrix)) %>%
          mutate(vaf = as.numeric(sub("%","", vaf))/100)      %>%
          spread(key = sample, value = vaf)
        freq_matrix = as.matrix.data.frame(tmp %>% select(-c(mut)))
        rownames(freq_matrix) = tmp$mut
    # --- TAG = VR-DP --------------------------------------------------------------
        } else if (TAG=="VR-DP") {
          vr_matrix <- extract.gt(vcf, element = "VR",
                                  mask = FALSE, as.numeric = TRUE,
                                  return.alleles = FALSE,IDtoRowNames = TRUE,
                                  extract = TRUE,convertNA = TRUE)
          dp_matrix <- extract.gt(vcf,element = "DP",
                                  mask = FALSE, as.numeric = TRUE,
                                  return.alleles = FALSE,IDtoRowNames = TRUE,
                                  extract = TRUE,convertNA = TRUE)
          freq_matrix = as.matrix(vr_matrix / dp_matrix)
    # --- TAG = U --------------------------------------------------------------
       } else if (TAG=="U") {
          AlternativeAlleles <- sub(",.*","",getALT(vcf))
          DP <- extract.gt(vcf,element = "DP", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames =     TRUE,extract = TRUE,convertNA = TRUE)        # DP: Total read depth
          AU <- extract.gt(vcf,element = "AU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames =     TRUE,extract = TRUE,convertNA = TRUE)
          CU <- extract.gt(vcf,element = "CU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames =     TRUE,extract = TRUE,convertNA = TRUE)
          GU <- extract.gt(vcf,element = "GU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames =     TRUE,extract = TRUE,convertNA = TRUE)
          TU <- extract.gt(vcf,element = "TU", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames =     TRUE,extract = TRUE,convertNA = TRUE)
          GetAlternativeCounts <- function(pos){
          ALT_ALL <- AlternativeAlleles[pos]
            alternative_counts_matrix <- get(paste(ALT_ALL, "U", sep=""))
            freq <- alternative_counts_matrix[pos,]/DP[pos,]
            return(freq)}
          a<-lapply(1:nrow(DP), FUN = GetAlternativeCounts)
          freq_matrix <- as.data.frame(do.call(rbind, a))
    # --- TAG = AD --------------------------------------------------------------
      } else if ( TAG == "AD" ) {
            freq_matrix <- extract.gt(vcf,element = TAG,
                                      mask = FALSE,as.numeric = FALSE,
                                      return.alleles = FALSE,IDtoRowNames = TRUE,
                                      extract = TRUE,convertNA = TRUE)
            freq_matrix <- AD_frequency(freq_matrix, delim = ",",
                                        allele = 2L, sum_type = 1L,
                                        decreasing = 1L)
      } else {
            freq_matrix <- extract.gt(vcf,element = TAG,
                                      mask = FALSE,as.numeric = FALSE,
                                      return.alleles = FALSE,IDtoRowNames = TRUE,
                                      extract = TRUE,convertNA = TRUE)
      }
      
    # Obtain freq_data
    freq_data <- data.frame(freq_matrix)
    freq_data$mut_info <- apply(getFIX(vcf)[,c("CHROM","POS","REF","ALT")],1,paste,collapse="_"    )
    freq_data$chrom_pos <- apply(getFIX(vcf)[,c("CHROM","POS")],1,paste,collapse="_")
    colnames(freq_data) <- gsub("\\.variant.*","",sub("\\.S.*","",colnames(freq_data)))
    freq_data <- freq_data %>% select(-one_of("RN","N"))

  #--------------------------------------------#    
  # Filter VAFs depending on the caller        #
  #--------------------------------------------#
  # https://github.com/hdetering/mseq-vc/blob/master/analysis/SRSV_create_analysis_db.util.R
  # https://github.com/hdetering/mseq-vc/blob/master/analysis/SRSV_create_analysis_db.R
    
    vcf_table <- read.vcfR(fn_vcf, verbose = F) %>% vcfR2tidy()
    FILTER  <-vcf_table$fix$FILTER
    freq_data <- cbind.data.frame(freq_data, FILTER)
    # --- MuTect1 Mutect2 Strelka1 Strelka2 NeuSomatic SNooPer'------------------
    if ( SOFTWARE %in% c('MuTect1', 'Mutect2', 'Strelka1', 'Strelka2', 'NeuSomatic','SNooPer') ) {
      freq_data <- freq_data %>%
      dplyr::filter( FILTER == 'PASS' ) %>%
      select (-FILTER)  %>%
      select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    # --- VarDict ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarDict' ) {
      STATUS <- vcf_table$fix$STATUS
      freq_data <- cbind.data.frame(freq_data, STATUS) %>%
      dplyr::filter( FILTER == 'PASS' & grepl('Somatic$', STATUS) ) %>%
      select (-FILTER,-STATUS)  %>%
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    # --- VarScan ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarScan' ) {
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_SS,Indiv) %>% 
        dplyr::mutate(chrom_pos=paste("chr",ChromKey,"_",POS,sep="")) %>%
        inner_join(freq_data, by=c('chrom_pos')) %>%
        dplyr::filter( FILTER == 'PASS') %>%
        dplyr::filter( gt_SS == 2 ) %>%
        select (-FILTER,-gt_SS,-ChromKey,-POS,-Indiv)  %>%
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    # --- MultiSNV --------------------------------------------------------------
    } else if ( SOFTWARE == 'MultiSNV' ) {
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_SS,Indiv) %>% 
        dplyr::mutate(chrom_pos=paste("chr",ChromKey,"_",POS,sep="")) %>%
        inner_join(freq_data, by=c('chrom_pos')) %>%
        dplyr::filter( FILTER == 'PASS' ) %>%
        dplyr::filter( gt_SS == 2 ) %>%
        select (-FILTER,-gt_SS,-ChromKey,-POS,-Indiv) %>%
        select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    # --- HaplotypeCaller, Bcftools  --------------------------------------------
    } else if (SOFTWARE %in% c('HaplotypeCaller','Bcftools')){
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,Indiv) %>%
        dplyr::mutate(chrom_pos=paste("chr",ChromKey,"_",POS,sep="")) %>%
        inner_join(freq_data, by=c('chrom_pos')) %>%
        dplyr::filter( gt_GT != '0/0' ) %>%
        select (-FILTER,-gt_GT,-ChromKey,-POS,-Indiv)  %>%
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    # --- Mutect2 multi-sample --------------------------------------------------
    } else if ( SOFTWARE == 'Mutect2_mseq' ) {
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,gt_AD,Indiv) %>%
        dplyr::mutate(chrom_pos=paste("chr",ChromKey,"_",POS,sep="")) %>%
        inner_join(freq_data, by=c('chrom_pos')) %>%
        dplyr::filter( FILTER == 'PASS' ) %>%
        dplyr::filter( gt_GT == '0/1' ) %>%
        rowwise() %>% 
        mutate( rc_alt = as.numeric(str_split(gt_AD, ',', simplify=T)[2]) ) %>%
        dplyr::filter( rc_alt > 0 ) %>%
        select (-FILTER,-gt_GT,-gt_AD,-ChromKey,-POS,-Indiv,-rc_alt) %>%
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    } else {
      freq_data <- freq_data %>% select(-FILTER)
    }
    
    output <- cbind.data.frame(freq_data, rep(SOFTWARE,nrow(freq_data)), rep(REPLICATE_ID,nrow(freq_data)))
    colnames(output) <- c("T1","T2","T3","T4","T5","mut_info","chrom_pos","software","replicate")
    return(output)
    
}


# Create a table with all the combinations of factors we want to get the distance of allele frequencies
CallersInfo <- readRDS(paste(data_dir,"df_caller.rds",sep=""))
name_caller_pub <- CallersInfo$name_caller[which(CallersInfo$name_caller!="MuClone_perf")]

TAGS <-  c(
  "AD", # Bcftools
  "PM", # CaVEMan
  "AD", # HaplotypeCaller
  NA, # MuClone
  "BCOUNT-DP", # MultiSNV
  "FA", # MuTect1
  "AF", # Mutect2_multi_F
  "AF", # Mutect2_single
  "AF", # NeuSomatic
  NA, # Shimmer
  "VR-DP", # SNooPer
  NA, # SNV-PPILP
  "BCOUNT-DP", # SomaticSniper
  "U", # Strelka2
  "AF", # VarDict
  "FREQ" # VarScan
)

caller_AFtags <- cbind.data.frame(name_caller_pub, TAGS)
colnames(caller_AFtags) <- c("name_caller_pub","AF_TAG")
Replicates_row<- as.data.frame(t(Replicates_Info$name_rep))
library(dplyr)
replicates_rep <- Replicates_row %>% slice(rep(1:n(), each = nrow(caller_AFtags)))
caller_AFtags_repli <- cbind.data.frame(caller_AFtags, replicates_rep)
caller_AFtags_repli_info <- reshape2::melt(caller_AFtags_repli, id.vars=c("name_caller_pub", "AF_TAG"))
caller_AFtags_repli_info <- caller_AFtags_repli_info[,-3]
colnames(caller_AFtags_repli_info) <- c("name_caller_pub","AF_TAG", "Replicate")
Strategy <- sapply(caller_AFtags_repli_info$Replicate, function(x)  Replicates_Info$ttype[match(x, Replicates_Info$name_rep)])


getVAFs(SOFTWARE = as.character(caller_AFtags_repli_info$name_caller_pub[2]),  TAG = as.character(caller_AFtags_repli_info$AF_TAG[2]), REPLICATE = as.character(caller_AFtags_repli_info$Replicate[2]))

AFs <- do.call("rbind", apply(caller_AFtags_repli_info, 1, function(x) getVAFs(SOFTWARE = x['name_caller_pub'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
saveRDS(AFs, file = paste(data_dir,"SRSV.AFs.rds",sep=""))

########################################
# OBTAIN TABLE WITH ALL SIMULATED VAFs #
########################################

True_VAFs_muts <- readRDS(paste(data_dir,"df_mut_sample.rds",sep=""))
mutations <- readRDS(paste(data_dir,"df_mut.rds",sep=""))
Replicates_Info <- readRDS(paste(data_dir,"df_rep.rds",sep=""))


df_vaf_exp <- True_VAFs_muts %>%
  inner_join( mutations, by = c('id_rep', 'id_mut') ) %>%
  inner_join( Replicates_Info, by = c('id_rep') ) %>%
  unite( 'chrom_pos', chrom, pos) %>%
  select( name_rep, id_sample, chrom_pos, vaf_exp, rc_tot, rc_alt )

########################################
#   MERGE SIMULATED AND REPORTED VAFs  #
########################################

AFs <-readRDS(paste(data_dir,"SRSV.AFs.rds",sep=""))
AFs$name_rep <- as.character( AFs$replicate )
AFs <- AFs %>% rename( R1 = T1, R2 = T2, R3 = T3, R4 = T4, R5 = T5 )
df_vaf_obs <- AFs %>% pivot_longer( cols = paste0('R', 1:5), names_to = 'id_sample', values_to = 'vaf' ) %>%
  mutate( vaf_obs = as.numeric(vaf) ) %>%
  select( -mut_info, -replicate, -vaf )

df_vaf <- df_vaf_obs %>% 
  dplyr::left_join( df_vaf_exp, by = c('name_rep', 'id_sample', 'chrom_pos') )


####
# EQUIVALENT PLOTS TO FIGURE S8
#SimulatedvsReportedVAFs %>%
#  mutate(TP=!is.na(SimulatedvsReportedVAFs$vaf_exp))
####


getDistances <- function(SOFTWARE,REPLICATE_ID) {
    
    print (SOFTWARE)
    print (REPLICATE_ID)
    if (SOFTWARE=="Mutect2_multi_F") {SOFTWARE="Mutect2_mseq"}
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


# Calculate the distance for each row
Distance <- apply(caller_AFtags_repli_info[,-2], 1, function(x) getDistances(SOFTWARE = x['name_caller_pub'],REPLICATE_ID = x['Replicate']))

AF_Distance_Results <- cbind.data.frame(caller_AFtags_repli_info, Distance)
saveRDS(AF_Distance_Results, file =  paste(data_dir,"df_AFdistances.rds",sep=""))

