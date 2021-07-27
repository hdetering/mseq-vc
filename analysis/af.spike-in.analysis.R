data_dir="/Users/tama/Google\ Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/RRSV/Results/"
data_root="/Users/tama/Downloads/VAFs/"

data_dir = "/Users/laura/GoogleDrive/LAB FOLDERS/M-seq Variant Calling Benchmarking/RRSV/Results/"
data_root = "/Users/laura/GoogleDrive/PROJECTS/PROYECTOS/TestingVC_withTama/07_WORKDIR/vcfs_vaf_tama/"


library(tidyverse)
library(vcfR)
library(magrittr)

####################################
#       Obtain reported VAFs      #
####################################

SOFTWARE="Bcftools"
REPLICATE_ID="S1.R1"
TAG="AD"
# freq_data["21_9484558",]


SOFTWARE="SomaticSniper"
REPLICATE_ID="S1.R1"
TAG="BCOUNT-DP"

######################################
### FOR SOMATIC SNIPER, WHEN THERE TWO ALTERNATIVE ALLELES CHECK THE REFERENCE GENOTYPE AND IF IT
# CONTAINS THE first alternative allele (genotype 0/1 and 1/1) The I will take VAF for the second allele
# TAKING ONLY SS=2 AND USE ALSO SET INFO  tag(variant)


### APPLY ALL CHANGES TO DE NOVO TOO

### Use SnooPerCorrected after Laura fix them  (gemres?)


### VarScan: remove % from FREQ!!!not working?  9483746

### SUBSTITUTE FUNCTION AD BY ANOTHER ONE

SOFTWARE="VarScan"
REPLICATE_ID="S1.R1"
TAG="FREQ"

varscan_af = getVAFs(SOFTWARE = SOFTWARE, TAG = TAG, REPLICATE_ID = REPLICATE_ID)

# (Laura) for MultiSNV, modify the vcfs removing foo

SOFTWARE="MultiSNV"
REPLICATE_ID="S1.R1"
TAG="BCOUNT-DP"
getVAFs(SOFTWARE = SOFTWARE, TAG = TAG, REPLICATE_ID = REPLICATE_ID)


msnv = getVAFs(SOFTWARE = SOFTWARE, TAG = TAG, REPLICATE_ID = REPLICATE_ID)
head(vcf2)
head(msnv)
vtama = msnv %>% separate(chrom_pos, into = c("chrom", "pos"), sep = "_") %>%
  gather(key = "id_sample", value = "vaf", T1:T5) %>%
  select(chrom, pos, id_sample, vaf) %>%
  mutate(tama = 1) %>% 
  filter(!is.na(vaf))   #### IMPORTANT. Otherwise, I keep 0-vaf variants (inferred as present by MultiSNV) and Tama doesn't
vlau = vcf2 %>% select(-c(id_caller, id_rep)) %>% mutate(lau = 1, chrom = as.character(chrom), pos = as.character(pos))
joinedmsnv = full_join(vtama, vlau)
joinedmsnv %>% group_by(tama, lau) %>% summarise(n())
joinedmsnv %>% filter(is.na(tama))

SOFTWARE="SNooPer"
REPLICATE_ID="S1.R1"
REPLICATE_ID="S2.R1"
TAG="VR-DP"
getVAFs(SOFTWARE = SOFTWARE, TAG = TAG, REPLICATE_ID = REPLICATE_ID)



SOFTWARE="HaplotypeCaller"
REPLICATE_ID="S1.R1"
TAG="AD"
hc = getVAFs(SOFTWARE = SOFTWARE, TAG = TAG, REPLICATE_ID = REPLICATE_ID)

SOFTWARE="VarDict"
REPLICATE_ID="S1.R1"
TAG="AF"

#######################################

getVAFs <- function(SOFTWARE, TAG, REPLICATE_ID){
      print(SOFTWARE)
      print(REPLICATE_ID)
      ## OBTAIN EMPIRICAL VAFs
      if (SOFTWARE=="SNooPer") {
        SOFTWARE="SNooPerGermresCorrected"
        } else if (SOFTWARE=="Strelka2"){
        SOFTWARE="Strelka2Corrected"
        } else if(SOFTWARE == "MultiSNV"){
          SOFTWARE= "MultiSNVCorrected"
        }
      fn_vcf <- paste(data_root, REPLICATE_ID,".",SOFTWARE,".PASS.vcf", sep="")
      if(SOFTWARE=="MultiSNVCorrected"){
        fn_vcf <- paste(data_root, REPLICATE_ID,".",SOFTWARE,".PASS_and_lowqual.vcf", sep="")
      }
      vcf <- read.vcfR(fn_vcf, verbose = F)
      if (is.na(TAG)){
        return(NULL)
      }else if (TAG == "BCOUNT-DP"){
        
old_way_laura = function(){
        allele_for_vaf = data.frame(alt_alleles = getALT(vcf),
                              ref_alleles = getREF(vcf),
                              healthy_genotype = extract.gt(vcf, element = "GT")[,1],
                              gt_T1 = extract.gt(vcf, element = "GT")[,2],
                              gt_T2 = extract.gt(vcf, element = "GT")[,3],
                              gt_T3 = extract.gt(vcf, element = "GT")[,4],
                              gt_T4 = extract.gt(vcf, element = "GT")[,5],
                              gt_T5 = extract.gt(vcf, element = "GT")[,6]) %>% 
          rowid_to_column("mut_id") %>%
          gather(key = "sample", value = "gt", gt_T1 ,gt_T2, gt_T3 ,gt_T4 ,gt_T5) %>% 
          filter(!is.na(gt)) %>% 
          group_by(mut_id, alt_alleles, ref_alleles, healthy_genotype) %>%
          dplyr::summarise(g00 = sum(gt == "0/0"),
                    g01 = sum(gt == "0/1"),
                    g02 = sum(gt == "0/2"),
                    g03 = sum(gt == "0/3"),
                    g11 = sum(gt == "1/1"),
                    g12 = sum(gt == "1/2"),
                    g13 = sum(gt == "1/3"),
                    g22 = sum(gt == "2/2"),
                    g23 = sum(gt == "2/3"),
                    g33 = sum(gt == "3/3"),
                    n_gts = n(),
                    n_dis_gts = n_distinct(gt)) %>% 
          arrange(mut_id) %>%
          gather(key = "gt", value = "gt_count", g00:g33) %>% 
          filter(gt_count>0) %>% 
          group_by(mut_id) %>%
          filter(gt_count == max(gt_count)) %>%
          mutate(id_gt_in_mut = row_number()) %>%
          filter(id_gt_in_mut == 1) %>%
          mutate(gt = str_remove(gt, "g")) %>%
          mutate(gt=gsub('(?<=.)(?=.)', '/', gt, perl=TRUE)) %>%
          separate(gt, into = c("allele_left", "allele_right"), sep = "/", remove = FALSE) %>%
          gather(key = "allele_type", value = "allele_id", allele_left:allele_right) %>% 
          filter(!str_detect(healthy_genotype, allele_id)) %>%
          group_by(mut_id) %>%
          filter(allele_id == min(allele_id)) %>%
          mutate(id_allele_in_mut = row_number()) %>%
          filter(id_allele_in_mut == 1) %>%
          select(mut_id, alt_alleles, ref_alleles, allele_id) %>%
          mutate(all_alleles = paste(ref_alleles, alt_alleles, sep = ",")) %>%
          separate_rows(all_alleles, sep = ",") %>%
          group_by(mut_id) %>%
          mutate(allele_id_temp = row_number()) %>%
          mutate(n_of_alleles = n()) %>%
          filter(allele_id_temp-1 == allele_id) %>%
          select(mut_id, allele_to_get_vaf = all_alleles)
          
        
        # # Which genotype combinations do we have?
        # 
        # all_the_samples = colnames(extract.gt(vcf, element = "GT"))
        # all_genotypes = as.data.frame(extract.gt(vcf, element = "GT"))
        # colnames(all_genotypes) = c("H", "T1", "T2", "T3", "T4", "T5")
        # all_genotypes %>% 
        #   filter(H != "0/0") %>% 
        #   filter(!is.na(T1) &!is.na(T2) &!is.na(T3) &!is.na(T4) &!is.na(T5)) %>%
        #   group_by(H, T1, T2, T3, T4, T5) %>%
        #   summarise(n())
        
        AlternativeAlleles = allele_for_vaf$allele_to_get_vaf
}   
      # AlternativeAlleles <- sub(",.*","",getALT(vcf))

laura_including_other_healthy_genotypes = function(){


      
      # EXTRACT COUNTS FOR EACH ALLELE

      # BCOUNT: Occurrence count for each base at this site (A,C,G,T)
      BCOUNT <- extract.gt(vcf,element = "BCOUNT", as.numeric = FALSE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
      A_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=1) %>%
        as.data.frame() %>%
        gather(key = "sample", value = "A") 
      C_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=2) %>%
        as.data.frame() %>%
        gather(key = "sample2", value = "C") 
      G_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=3)%>%
        as.data.frame() %>%
        gather(key = "sample3", value = "G")
      T_counts <- apply(BCOUNT, 2, FUN = getAlleles, letter=4)%>%
        as.data.frame() %>%
        gather(key = "sample4", value = "T")
      
      allele_counts = cbind(A_counts, C_counts, G_counts, T_counts) %>% 
        group_by(sample) %>%
        mutate(mut_id = row_number()) %>%
        ungroup() %>%
        select(mut_id,sample, A, C, G, T) %>%
        mutate(sample = str_remove(sample, "ChineseSon.")) %>%
        filter(sample!= "H") %>%
        separate(sample, into = c("sample", NA, NA)) %>%
        gather(key = "allele", value = "count", A:T)
      
      # EXTRACT DEPTH FOR EACH POSITION
      # DP: Total read depth
      DP <- extract.gt(vcf,element = "DP", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
      
      depth = DP %>%
        as.data.frame() %>% 
        rowid_to_column("mut_id") %>% 
        gather(key = "sample", value = "depth", 2:(ncol(DP)+1) ) %>% 
        mutate(sample = str_remove(sample, "ChineseSon.")) %>%
        filter(sample!= "H") %>%
        separate(sample, into = c("sample", NA, NA)) 
      # DECIDE WHICH ALLELES ARE MUTATIONS AND MERGE WITH COUNTS AND DEPTH
      
      allele_for_vaf = data.frame(alt_alleles = getALT(vcf),
                                  ref_alleles = getREF(vcf),
                                  healthy_genotype = extract.gt(vcf, element = "GT")[,1],
                                  gt_T1 = extract.gt(vcf, element = "GT")[,2],
                                  gt_T2 = extract.gt(vcf, element = "GT")[,3],
                                  gt_T3 = extract.gt(vcf, element = "GT")[,4],
                                  gt_T4 = extract.gt(vcf, element = "GT")[,5],
                                  gt_T5 = extract.gt(vcf, element = "GT")[,6])   %>%
        rowid_to_column("mut_id") %>%
        gather(key = "sample", value = "tumor_genotype", gt_T1: gt_T5) %>%
        filter(!is.na(tumor_genotype)) %>%
        mutate(sample = str_remove(sample, "gt_")) %>%
        mutate(tumor_allele = tumor_genotype) %>%
        separate_rows(tumor_allele, sep = "/") %>%
        filter(str_count(healthy_genotype, tumor_allele) < str_count(tumor_genotype, tumor_allele)) %>%
        mutate(counts_allele_in_healthy = str_count(healthy_genotype, tumor_allele) ) %>%
        mutate(all_alleles = paste(ref_alleles, alt_alleles, sep = ",")) %>%
        separate_rows(all_alleles, sep = ",") %>%
        group_by(mut_id, alt_alleles, ref_alleles, healthy_genotype, sample, tumor_genotype, tumor_allele) %>%
        mutate(allele_id = row_number()-1)  %>%
        filter(tumor_allele == allele_id) %>% arrange(mut_id) %>%
        rename(allele_to_get_vaf = all_alleles) %>%
        ungroup() %>%
        select(mut_id, sample, allele_to_get_vaf,healthy_genotype, tumor_genotype, counts_allele_in_healthy ) %>%
        rename(allele = allele_to_get_vaf) %>%
        left_join(allele_counts) %>%
        left_join(depth) %>%
        mutate(VAF = count/depth ) %>%
        mutate(somatic_VAF = VAF - (0.5*counts_allele_in_healthy))  %>%# counts_allele_in_healthy always zero in the test file, check with others
        mutate(somatic_VAF = if_else(somatic_VAF<0, 0, somatic_VAF))  # In case it becomes negative 
}
# letter: 1>A 2>C 3>G 4>T

AlternativeAlleles <- sub(",.*","",getALT(vcf))
# BCOUNT: Occurrence count for each base at this site (A,C,G,T)
BCOUNT <- extract.gt(vcf,element = "BCOUNT", as.numeric = FALSE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
# DP: Total read depth
DP <- extract.gt(vcf,element = "DP", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
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
      
      # Keep only 0/0 in the healthy   (moved to later otherwise it fails when wiltering for PASS I think)
      
      # healthy_genotype = extract.gt(vcf, element = "GT")[,1] 
      # freq_matrix = freq_matrix[which(healthy_genotype=="0/0"),]
      
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
     tmp = freq_matrix %>% data.frame(stringsAsFactors = FALSE) %>%
       dplyr::mutate(mut = rownames(freq_matrix)) %>%
       tidyr::gather(key = "sample", value = "vaf", 1:ncol(freq_matrix)) %>%
        dplyr::mutate(vaf = as.numeric(sub("%","", vaf))/100)      %>%
        tidyr::spread(key = sample, value = vaf)
      freq_matrix <- as.matrix.data.frame(tmp %>% select(-c(mut)))
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
    
    
    ###### THIS FUNCTION BELOW IS NOT DOING WHAT I WAS EXPECTING
    # CHANGE ALSO IN de novo!!!!
    # freq_matrix <- AD_frequency(freq_matrix, delim = ",", allele = 2L, sum_type = 1L, decreasing = 1L)
    
    # Customized extraction (Laura)  # works for HC. Bcftools?
    
    if(SOFTWARE=="HaplotypeCaller"){
    freq_df = freq_matrix %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "mut_id") %>%
      gather(key = "sample", value = "ads", 2:7) %>%
      filter(!str_detect(sample, "H")) %>%
      filter(!str_detect(sample, "RN")) %>%
      separate(ads, into = c("ref_ad", "alt_ad"), extra = "drop")  %>%
      mutate(ref_ad = as.numeric(ref_ad),
             alt_ad = as.numeric(alt_ad)) %>%
      mutate(vaf = if_else(alt_ad+ref_ad >0, 
                           alt_ad/(alt_ad+ref_ad), 0)) %>%
      select(-c(ref_ad, alt_ad)) %>%
      spread(key= sample, value = vaf) 
    
    }else if(SOFTWARE == "Bcftools"){
      freq_df = freq_matrix %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "mut_id") %>%
        gather(key = "sample", value = "ads", 2:6) %>%
        separate(ads, into = c("ref_ad", "alt_ad"), extra = "drop")  %>%
        mutate(ref_ad = as.numeric(ref_ad),
               alt_ad = as.numeric(alt_ad)) %>%
        mutate(vaf = if_else(alt_ad+ref_ad >0, 
                             alt_ad/(alt_ad+ref_ad), 0)) %>%
        select(-c(ref_ad, alt_ad)) %>%
        spread(key= sample, value = vaf) 
      
    }
    
    freq_matrix = as.matrix(freq_df)
    rownames(freq_matrix) = freq_df$mut_id
    freq_matrix = freq_matrix[,2:ncol(freq_matrix)]
      
    ###
    
    ######
  
    
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
    if(SOFTWARE=="VarScan" | SOFTWARE == "HaplotypeCaller"){   # CHECK IF MORE CALLERS NEED THIS. MAYBE BCFTOOLS
      freq_data = freq_data %>% mutate(chrom_pos = rownames(freq_data)) %>%
        separate(chrom_pos, into = c("chrom", "pos"), remove = FALSE) %>% 
        mutate(pos = as.numeric(pos)) %>%
        arrange(pos) %>%
        select(-c(chrom_pos, chrom, pos))
    }
    freq_data$mut_info <- apply(getFIX(vcf)[,c("CHROM","POS","REF","ALT")],1,paste,collapse="_"    ) #### WRONG FOR VARSCAN and HaplotypeCaller. Workaround above
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
    if ( SOFTWARE %in% c('MuTect1', 'Mutect2', 'Strelka1', 'Strelka2', 'NeuSomatic','SNooPerGermresCorrected') ) {
      freq_data <- freq_data %>%
      dplyr::filter( FILTER == 'PASS' ) %>%
      select (-FILTER)
    # --- VarDict ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarDict' ) {
      STATUS <- vcf_table$fix$STATUS
      TYPE <- vcf_table$fix$TYPE
      freq_data <- cbind.data.frame(freq_data, STATUS) %>%
        cbind.data.frame(TYPE) %>%
      dplyr::filter( FILTER == 'PASS' & grepl('Somatic$', STATUS) & TYPE =="SNV") %>%
      select (-FILTER,-STATUS, -TYPE)
    # --- VarScan ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarScan' ) {
      freq_data <- vcf_table$gt %>% dplyr::select(ChromKey,POS,gt_SS,Indiv) %>%
        dplyr::mutate(ChromKey=21) %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        dplyr::inner_join(freq_data, by=c('chrom_pos')) %>%
        dplyr::filter( FILTER == 'PASS') %>%
        dplyr::filter( gt_SS == 2 ) %>%
        dplyr::select(-FILTER,-gt_SS,-ChromKey,-POS,-Indiv) %>%
        dplyr::select(ChineseSon.T1,ChineseSon.T2,ChineseSon.T3,ChineseSon.T4,ChineseSon.T5,mut_info,chrom_pos)
    # --- MultiSNV --------------------------------------------------------------
    } else if ( SOFTWARE == 'MultiSNVCorrected' ) {
        freq_data <- freq_data %>%
        reshape2::melt(id.vars=c("chrom_pos","mut_info","FILTER"))
        colnames(freq_data)[4] <- "Indiv"
      SS <- extract.gt(vcf,element = "SS", as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)[,-1]
      SS <- cbind(rownames(SS), data.frame(SS, row.names=NULL))
      colnames(SS)[1] <- "chrom_pos"
      freq_data <- SS %>%
        reshape2::melt(variable.name="Indiv",value.name="gt_SS") %>%
        dplyr::inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter( FILTER == 'PASS' | FILTER =="LOW_QUAL") %>%   # Laura added LOW_QUAL
        dplyr::filter( gt_SS == 2 ) %>%
        dplyr::select (-FILTER,-gt_SS)  %>%
        spread(key = Indiv,value = value, convert = FALSE) %>%
        dplyr::select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    # --- HaplotypeCaller  --------------------------------------------
    } else if (SOFTWARE %in% c('HaplotypeCaller')){
      INDEL<-is.indel(vcf)
      freq_data2 <-cbind.data.frame(freq_data,INDEL) %>%
        gather(key = "Indiv", value = "value", T1:T5)    # (?) Laura 
        # reshape2::melt(variable.name="Indiv")
      
      
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,Indiv) %>%
        dplyr::mutate(ChromKey=21) %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        dplyr::filter(Indiv!="ChineseSon.H") %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>% 
        inner_join(freq_data2, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter(gt_GT!='0/0' & gt_GT!='0|0') %>%   # Laura: change | to &
        dplyr::filter( INDEL == FALSE ) %>%
        # select (-c(FILTER,gt_GT,ChromKey,POS,INDEL))
        select (-c(gt_GT,ChromKey,POS,INDEL))  %>%  # Laura
        spread(key = Indiv,value = value, convert = FALSE) %>%
        select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    # --- Bcftools -------------------------------------------------------------
    } else if (SOFTWARE %in% c('Bcftools')) {
      freq_data <- reshape2::melt(freq_data,id.vars=c("mut_info","chrom_pos"),variable.name="Indiv")
      freq_data <- vcf_table %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>%
        inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        #mutate(value = ifelse (gt_GT!='0/0',NA,value)) %>%
        dplyr::filter( gt_GT != '0/0' & gt_GT!='0|0' ) %>%    # Laura: change | to &
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
        dplyr::filter( gt_GT == '0/1' | gt_GT == '0|1' | gt_GT == '0/1/2' | gt_GT == '0|1/2' | gt_GT == '0/1|2') %>% # We don't want 0/2 or 0|2
        rowwise() %>% 
        mutate( rc_alt = as.numeric(str_split(gt_AD, ',', simplify=T)[2]) ) %>%
        dplyr::filter( rc_alt > 0 ) %>%
        select (-FILTER,-gt_GT,-gt_AD,-ChromKey,-POS,-rc_alt) %>%
        spread(key = Indiv,value=value) %>%
        select (T1,T2,T3,T4,T5,mut_info,chrom_pos)
    } else {
      freq_data <- freq_data %>% select(-FILTER)
      
      if(SOFTWARE=="SomaticSniper"){
      # Keep only 0/0 in the healthy   (moved to later otherwise it fails when wiltering for PASS I think)
      
      healthy_genotype = extract.gt(vcf, element = "GT")[,1]
      freq_data = freq_data[which(healthy_genotype=="0/0"),]
        
      }
      
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
Replicates_Info <- readRDS(paste(data_dir,"RRSV.replicates.rds",sep=""))
Replicates_row<- as.data.frame(t(Replicates_Info$name_rep))
replicates_rep <- Replicates_row %>% slice(rep(1:n(), each = nrow(caller_AFtags)))
caller_AFtags_repli <- cbind.data.frame(caller_AFtags, replicates_rep)
caller_AFtags_repli_info <- reshape2::melt(caller_AFtags_repli, id.vars=c("name_caller", "AF_TAG"))
caller_AFtags_repli_info <- caller_AFtags_repli_info[,-3]
colnames(caller_AFtags_repli_info) <- c("name_caller", "AF_TAG", "Replicate")
Strategy <- sub("\\..*","",caller_AFtags_repli_info$Replicate)


# Obtain the AFs
AFs <- do.call("rbind", apply(caller_AFtags_repli_info, 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
AFs_multisnv <- do.call("rbind", apply(caller_AFtags_repli_info %>% filter(name_caller=="MultiSNV"), 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
AFs_haplotypecaller <- do.call("rbind", apply(caller_AFtags_repli_info %>% filter(name_caller=="HaplotypeCaller"), 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
AFs_vardict <- do.call("rbind", apply(caller_AFtags_repli_info %>% filter(name_caller=="VarDict"), 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
AFs_bcftools <- do.call("rbind", apply(caller_AFtags_repli_info %>% filter(name_caller=="Bcftools"), 1, function(x) getVAFs(SOFTWARE = x['name_caller'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))

# saveRDS(AFs, file = paste(data_dir,"RRSV.AFs.rds",sep=""))
# saveRDS(AFs, file = paste(data_dir,"RRSV.AFs.test20201022.rds",sep=""))
# saveRDS(AFs, file = paste(data_dir,"RRSV.AFs.test20201022v2.rds",sep=""))
saveRDS(AFs, file = paste(data_dir,"RRSV.AFs.test20201022v3.rds",sep=""))
####################################
#      Obtain simulated VAFs       #
####################################

True_VAFs_muts<- readRDS(paste(data_dir,"RRSV.mutations_samples.rds",sep=""))
mutations <- readRDS(paste(data_dir,"RRSV.mutations.rds",sep=""))
#Replicates_Info <- readRDS(paste(data_dir,"RRSV.replicates.rds",sep=""))


df_vaf_exp <- True_VAFs_muts %>%
  inner_join( mutations, by = c('id_rep', 'id_mut') ) %>%
  inner_join( Replicates_Info, by = c('id_rep') ) %>%
  tidyr::unite( 'chrom_pos', chrom, pos) %>%
  mutate( vaf_exp = as.numeric(vaf_exp) ) %>%
  select( name_rep, id_sample, chrom_pos, vaf_exp)



#####################################
# Merge reported and simulated VAFs #
#####################################

AFs <-readRDS(paste(data_dir,"RRSV.AFs.rds",sep=""))
AFs$name_rep <- as.character( AFs$replicate )
df_vaf_obs <- AFs %>% tidyr::pivot_longer( cols = paste0('T', 1:5), names_to = 'id_sample', values_to = 'vaf' ) %>%
  mutate( vaf_obs = as.numeric(vaf) ) %>%
  select( -mut_info, -replicate, -vaf )

df_vaf <- df_vaf_obs %>% 
  dplyr::full_join( df_vaf_exp, by = c('name_rep', 'id_sample', 'chrom_pos') )


getDistances <-  function(SOFTWARE, REPLICATE_ID){
  
  print (SOFTWARE)
  print (REPLICATE_ID)
  True_VAFs_muts_selected <- df_vaf %>%
    dplyr::filter(grepl(REPLICATE_ID,name_rep)) %>%
    dplyr::filter(grepl(SOFTWARE,software))
  
  #######################################################################
  # Substitute simulated NAs by 0s at FPs and observed NAs by 0s at TNs #
  #######################################################################
  True_VAFs_muts_selected[is.na(True_VAFs_muts_selected$vaf_exp),"vaf_exp"] <- 0
  # Line to remove TNs
  #True_VAFs_muts_selected <- True_VAFs_muts_selected[!is.na(True_VAFs_muts_selected$vaf_obs),]
  # At the end we decided to keep TNs
  True_VAFs_muts_selected[is.na(True_VAFs_muts_selected$vaf_obs),"vaf_obs"] <- 0
  
  
  matrixfordist <- rbind(True_VAFs_muts_selected$vaf_obs, True_VAFs_muts_selected$vaf_exp)
  distance <- dist(matrixfordist, method = "euclidean") 
  return(distance)
  
}


# Calculate distances
Distance <- apply(caller_AFtags_repli_info, 1, function(x) getDistances(SOFTWARE = x['name_caller'], REPLICATE_ID = x['Replicate']))

AF_Distance_Results <- cbind.data.frame(caller_AFtags_repli_info, Strategy, Distance)
saveRDS(AF_Distance_Results, file = paste(data_dir,"RRSV.AFdistances.rds",sep=""))

