data_dir  = "/Users/laura/GoogleDrive/LAB FOLDERS/M-seq Variant Calling Benchmarking/de-novo/data/"
data_root="/Users/laura/GoogleDrive/PROJECTS/PROYECTOS/TestingVC_withTama/07_WORKDIR/VCFs_SRSV/"
library(vcfR)

########################
# OBTAIN REPORTED VAFs #
########################
# Create a function to obtain AFs per caller per replicate

SOFTWARE="VarDict"
caller =SOFTWARE
REPLICATE_ID="2a1730b2-bb2c-11e8-96d1-80000048fe80"
TAG="AF"
problematic_pos = "chr1_1545663"

SOFTWARE="Bcftools"
caller =SOFTWARE
REPLICATE_ID="28ca39a2-bb2c-11e8-9dcd-80000048fe80"
TAG="AD"

SOFTWARE="SNooPer"
caller =SOFTWARE
REPLICATE_ID="e9213d38-c241-11e8-a1c5-80000048fe80"
TAG="VR-DP"

getVAFs(SOFTWARE, TAG, REPLICATE_ID)
  


tag_caller = caller_AFtags_repli_info %>% 
  group_by(name_caller_pub, AF_TAG) %>% 
  dplyr::summarise()

df = data.frame(caller = character, columns = integer())

for(caller in tag_caller$name_caller_pub){
  print(caller)
t = getVAFs(SOFTWARE = caller,
            TAG = tag_caller[which(name_caller_pub == caller),"AF_TAG"][[1]], 
            REPLICATE_ID = REPLICATE_ID)
columns = ncol(t)
if(is.null(columns)){columns = 0}
df = rbind(df, data.frame(caller = caller, columns = columns))
}

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
                                       SNooPer = 'snooper.fix.vcf',  # Changed from snooper.corrected.vcf to snooper.fix.vcf 2/4/21
                                       SomaticSniper = 'somaticsniper.vcf',
                                       Strelka2 = 'strelka2.fix.vcf', # Changed from strelka2.corrected.vcf to strelka2.fix.vcf 2/4/21
                                       Shimmer = 'shimmer.vcf',
                                       VarDict = 'vardict.vcf',
                                       VarScan = 'varscan.vcf',
                                       MultiSNV = 'multisnv.vcf',  # already includes LOW_QUAL
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
            freq_matrix <- extract.gt(vcf,element = TAG,
                                      mask = FALSE,as.numeric = FALSE,
                                      return.alleles = FALSE,IDtoRowNames = TRUE,
                                      extract = TRUE,convertNA = TRUE)
      }
      
    # Obtain freq_data
      freq_data <- data.frame(freq_matrix)
      if(SOFTWARE=="VarScan" | SOFTWARE == "HaplotypeCaller" | SOFTWARE == "Bcftools"){   # CHECK IF MORE CALLERS NEED THIS. MAYBE BCFTOOLS
        freq_data= freq_data %>% mutate(chrom_pos = rownames(freq_data)) %>%
          separate(chrom_pos, into = c("chrom", "pos"), remove = FALSE) %>% 
          mutate(pos = as.numeric(pos)) %>%
          # mutate(chrom_pos_backup = chrom_pos
          arrange(pos) %>%
          select(-c(chrom_pos, chrom, pos))
      }
      freq_data$mut_info <- apply(getFIX(vcf)[,c("CHROM","POS","REF","ALT")],1,paste,collapse="_"    ) #### WRONG FOR VARSCAN and HaplotypeCaller. Workaround above
      freq_data$chrom_pos <- apply(getFIX(vcf)[,c("CHROM","POS")],1,paste,collapse="_")
      if(SOFTWARE == "NeuSomatic"){ # Added  2/4/21
        colnames(freq_data) <- gsub("\\.variant","",sub("\\.S.*","",colnames(freq_data)))
        colnames(freq_data) <- gsub("TUMOR$","TUMOR1",colnames(freq_data))
        colnames(freq_data) <- gsub("TUMOR","R",colnames(freq_data))
        
        
      }else{
      colnames(freq_data) <- gsub("\\.variant.*","",sub("\\.S.*","",colnames(freq_data)))
      }
      
      if (SOFTWARE!="Bcftools" & SOFTWARE != "NeuSomatic" & SOFTWARE != "SNooPer" ){ # Added NeuSomatic, SNooPer  2/4/21
        freq_data <- freq_data %>%
          dplyr::select(-one_of("RN","N","ChineseSon.H","NORMAL","NORMAL2","NORMAL3","NORMAL4","NORMAL5"))
      }

  #--------------------------------------------#    
  # Filter VAFs depending on the caller        #
  #--------------------------------------------#
  # https://github.com/hdetering/mseq-vc/blob/master/analysis/SRSV_create_analysis_db.util.R
  # https://github.com/hdetering/mseq-vc/blob/master/analysis/SRSV_create_analysis_db.R
    
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
      callers_without_set = c("Bcftools", "HaplotypeCaller", "MultiSNV", "Mutect2_multi_F", "Mutect2_mseq")
      
      if( ! SOFTWARE %in% callers_without_set){
      set <- vcf_table$fix$set
      freq_data_with_set  =  cbind.data.frame(freq_data, set) 
      }
      
    # --- MuTect1 Mutect2 Strelka1 Strelka2 NeuSomatic SNooPer'------------------
    if ( SOFTWARE %in% c('MuTect1', 'Mutect2', 'Mutect2_single', 'Strelka1', 'Strelka2', 'NeuSomatic','SNooPer') ) {
      freq_data <- freq_data %>%
      dplyr::filter( FILTER == 'PASS' ) %>%
      select (-FILTER)  %>%
      select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    # --- VarDict ---------------------------------------------------------------
    } else if ( SOFTWARE == 'VarDict' ) {
      
      STATUS <- vcf_table$fix$STATUS
      TYPE <- vcf_table$fix$TYPE
      
      
      freq_data <- cbind.data.frame(freq_data, STATUS) %>%
        cbind.data.frame(TYPE) %>%
        dplyr::filter( FILTER == 'PASS' & grepl('Somatic$', STATUS) & TYPE =="SNV") %>%
        select (-FILTER,-STATUS, -TYPE) %>%
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
        gather(key = "Indiv", value = "value", R1:R5)    # (?) Laura 
      # reshape2::melt(variable.name="Indiv")
      

      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,Indiv) %>%
        dplyr::mutate(ChromKey="chr1") %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        dplyr::filter(Indiv!="RN") %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>% 
        inner_join(freq_data2, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter(gt_GT!='0/0' & gt_GT!='0|0') %>%   # Laura: change | to &
        dplyr::filter( INDEL == FALSE ) %>%
        # select (-c(FILTER,gt_GT,ChromKey,POS,INDEL))
        select (-c(gt_GT,ChromKey,POS,INDEL))  %>%  # Laura
        spread(key = Indiv,value = value, convert = FALSE) %>% 
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
      # --- Bcftools -------------------------------------------------------------
    } else if (SOFTWARE %in% c('Bcftools')) {
      freq_data <- reshape2::melt(freq_data,id.vars=c("mut_info","chrom_pos"),variable.name="Indiv")
      freq_data <- vcf_table %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>% 
        inner_join(freq_data, by=c('chrom_pos','Indiv')) %>% 
        #mutate(value = ifelse (gt_GT!='0/0',NA,value)) %>%
        dplyr::filter( gt_GT != '0/0' & gt_GT!='0|0' ) %>%    # Laura: change | to &
        # dplyr::filter( grepl('^21', chrom_pos)) %>% 
        select(-gt_GT) %>%
        spread(key = Indiv,value=value) %>% 
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    # --- Mutect2 multi-sample --------------------------------------------------
    } else if ( SOFTWARE == 'Mutect2_mseq' ) {
      freq_data <- reshape2::melt(freq_data,id.vars=c("mut_info","chrom_pos","FILTER"),variable.name="Indiv")
      freq_data <- vcf_table$gt %>% select(ChromKey,POS,gt_GT,gt_AD,Indiv) %>%
        dplyr::mutate(ChromKey="chr1") %>%
        dplyr::mutate(chrom_pos=paste(ChromKey,"_",POS,sep="")) %>%
        dplyr::filter(Indiv!="RN") %>%
        mutate(Indiv=sub("\\..*","",Indiv)) %>%
        inner_join(freq_data, by=c('chrom_pos','Indiv')) %>%
        dplyr::filter( FILTER == 'PASS' ) %>%
        dplyr::filter( gt_GT == '0/1' | gt_GT == '0|1' | gt_GT == '0/1/2' | gt_GT == '0|1/2' | gt_GT == '0/1|2') %>% # We don't want 0/2 or 0|2
        rowwise() %>% 
        mutate( rc_alt = as.numeric(str_split(gt_AD, ',', simplify=T)[2]) ) %>%
        dplyr::filter( rc_alt > 0 ) %>%
        select (-FILTER,-gt_GT,-gt_AD,-ChromKey,-POS,-rc_alt) %>%
        spread(key = Indiv,value=value) %>%
        select (R1,R2,R3,R4,R5,mut_info,chrom_pos)
    } else {
      freq_data <- freq_data %>% select(-FILTER)
      
      if(SOFTWARE=="SomaticSniper"){
        # Keep only 0/0 in the healthy   (moved to later otherwise it fails when wiltering for PASS I think)
        
        healthy_genotype = extract.gt(vcf, element = "GT")[,6:10]%>%
          as.data.frame() %>%
          rowid_to_column("mut_id") %>%
          gather(key = "healthy_sample", value = "genotype", RN.variant:RN.variant5) %>%
          group_by(mut_id)   %>%
          dplyr::summarise(healthy_genotype_sum = min(genotype, na.rm = TRUE)) %>%
          mutate(healthy_genotype_sum = if_else(is.na(healthy_genotype_sum), "0/0", healthy_genotype_sum))
          
        freq_data = freq_data[which(healthy_genotype$healthy_genotype_sum =="0/0"),]
        
      }
    }
      
      
      if( (! SOFTWARE %in% callers_without_set) & nrow(freq_data) > 0){
      freq_data = freq_data %>% left_join(freq_data_with_set %>% select(mut_info, set)) %>%
        unique() %>%
        filter(!set == "FilteredInAll") %>%
        gather(key = "sample", value = "vaf", 1:5) %>%
        mutate(samplename_to_detect = if_else(sample == "R1", "variant", paste0("variant", str_remove(sample, "R")))) %>%
        mutate(vaf = as.character(vaf)) %>%
        mutate(vaf = case_when(set  == "Intersection" ~ vaf,
                               set == "FilteredInAll" ~ "NA",
                               str_detect(set, paste0("-",samplename_to_detect, "-")) | 
                                 str_detect(set, paste0("^",samplename_to_detect, "-")) |
                                 str_detect(set, paste0("-",samplename_to_detect, "$")) |
                                 str_detect(set, paste0("^",samplename_to_detect, "$")) ~ vaf,
                               is.na(vaf) ~ "NA",
                               TRUE ~ "NA"
        )
        ) %>% 
        select(-c(set, samplename_to_detect)) %>%
        mutate(vaf = na_if(vaf, "NA")) %>%
        spread(key = sample, value = vaf) %>% 
        select("R1","R2","R3","R4","R5","mut_info","chrom_pos")
      }
    
    if(nrow(freq_data) == 0){
      return(data.frame(T1 = character(), 
                        T2 = character(), 
                        T3=character(), 
                        T4 = character(), 
                        T5=character(), 
                        mut_info = character(), 
                        chrom_pos = character(), 
                        software = character(), 
                        replicate = character()))
    }
    output <- cbind.data.frame(freq_data, rep(SOFTWARE,nrow(freq_data)), rep(REPLICATE_ID,nrow(freq_data)))
    colnames(output) <- c("T1","T2","T3","T4","T5","mut_info","chrom_pos","software","replicate")
    return(output)
    
}

 
# Create a table with all the combinations of factors we want to get the distance of allele frequencies
CallersInfo <- readRDS(paste(data_dir,"df_caller.rds",sep=""))
CallersInfo = rbind(rbind(CallersInfo[1:6,], data.frame(id_caller = 7, name_caller = "Mutect2_single", class= "marginal")), CallersInfo[7:nrow(CallersInfo), ])
name_caller_pub <- CallersInfo$name_caller[which(CallersInfo$name_caller!="MuClone_perf")]
    ####### Mutect2_single missing in df_caller.rds for some reason!!!!
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
Replicates_Info <- readRDS(paste(data_dir,"df_rep.rds",sep=""))
Replicates_row<- as.data.frame(t(Replicates_Info$name_rep))
replicates_rep <- Replicates_row %>% slice(rep(1:n(), each = nrow(caller_AFtags)))
caller_AFtags_repli <- cbind.data.frame(caller_AFtags, replicates_rep)
caller_AFtags_repli_info <- reshape2::melt(caller_AFtags_repli, id.vars=c("name_caller_pub", "AF_TAG"))
caller_AFtags_repli_info <- caller_AFtags_repli_info[,-3]
colnames(caller_AFtags_repli_info) <- c("name_caller_pub","AF_TAG", "Replicate")
Strategy <- sapply(caller_AFtags_repli_info$Replicate, function(x)  Replicates_Info$ttype[match(x, Replicates_Info$name_rep)])


AFs <- do.call("rbind", 
               apply(caller_AFtags_repli_info, 1, function(x) getVAFs(SOFTWARE = x['name_caller_pub'],TAG = x['AF_TAG'], REPLICATE_ID = x['Replicate'])))
saveRDS(AFs, file = paste(data_dir,"SRSV.AFs.20210402.rds",sep=""))

########################################
# OBTAIN TABLE WITH ALL SIMULATED VAFs #
########################################

True_VAFs_muts <- readRDS(paste(data_dir,"df_mut_sample.rds",sep=""))
mutations <- readRDS(paste(data_dir,"df_mut.rds",sep=""))
#Replicates_Info <- readRDS(paste(data_dir,"df_rep.rds",sep=""))


df_vaf_exp <- True_VAFs_muts %>%
  inner_join( mutations, by = c('id_rep', 'id_mut') ) %>%
  inner_join( Replicates_Info, by = c('id_rep') ) %>%
  tidyr::unite( 'chrom_pos', chrom, pos) %>%
  select( name_rep, id_sample, chrom_pos, vaf_exp )

########################################
#   MERGE SIMULATED AND REPORTED VAFs  #
########################################

AFs <-readRDS(paste(data_dir,"SRSV.AFs.rds",sep=""))
AFs$name_rep <- as.character( AFs$replicate )
AFs <- AFs %>% rename( R1 = T1, R2 = T2, R3 = T3, R4 = T4, R5 = T5 )
df_vaf_obs <- AFs %>% tidyr::pivot_longer( cols = paste0('R', 1:5), names_to = 'id_sample', values_to = 'vaf' ) %>%
  mutate( vaf_obs = as.numeric(vaf) ) %>%
  select( -mut_info, -replicate, -vaf )

df_vaf <- df_vaf_obs %>% 
  dplyr::full_join( df_vaf_exp, by = c('name_rep', 'id_sample', 'chrom_pos') )



getDistances <- function(SOFTWARE,REPLICATE_ID) {
    
    print (SOFTWARE)
    print (REPLICATE_ID)
    if (SOFTWARE=="Mutect2_multi_F") {SOFTWARE="Mutect2_mseq"}
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


# Calculate the distance for each row
Distance <- apply(caller_AFtags_repli_info[,-2], 1, function(x) getDistances(SOFTWARE = x['name_caller_pub'],REPLICATE_ID = x['Replicate']))

AF_Distance_Results <- cbind.data.frame(caller_AFtags_repli_info, Distance)
saveRDS(AF_Distance_Results, file =  paste(data_dir,"df_AFdistances.rds",sep=""))

