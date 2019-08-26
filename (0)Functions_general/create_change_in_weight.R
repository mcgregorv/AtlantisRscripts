create_change_in_weight <- function(outdir, fgfile, biolprm, ncout, startyear, toutinc){
  
  
  # outdir=outPath
  # startyear=modelStartYear
  # toutinc=daysTimeStep
  # diet = TRUE
  # ncout<-"output"
  # fgfile<-paste(thisPath,"TBGB_Groups.csv",sep="")
  # biolprm<-paste(thisPath,"TBGB_biol",this_out,".prm",sep="")
  # 
  # contants
  nsecs <- 86400
  ndays <- 365
  g_per_ton <- 1e6
  species <- c("BIRD", "FISH", "MAMMAL", "SHARK")
  tons <- function(mgN) return(mgN * 5.7 * 20 / 1e9)
  
  cat("### ------------ Reading in data                                         ------------ ###\n")
  nc_out <- nc_open(paste(outdir, ncout, ".nc", sep = ""))
  # prod_out <- nc_open(paste(outdir, ncout, "PROD.nc", sep = ""))
  # bio_agg <- read.table(paste(outdir, ncout, "BoxBiomass.txt", sep = ""), header = T)
  # ssb <- read.table(paste(outdir, ncout, "SSB.txt", sep = ""), header = TRUE)
  # yoy <- read.table(paste(outdir, ncout, "YOY.txt", sep = ""), header = TRUE)
  bgm <- readLines(paste(outdir, grep(".bgm",dir(outdir), value = T), sep = ""))
  # biomass <- read.table(paste(outdir, ncout, "BiomIndx.txt", sep = ""), header = T)
  # rel_bio <- biomass[,c(1, grep("Rel",colnames(biomass)))]
  # tot_bio <- biomass[,c(1:(grep("Rel",colnames(biomass))[1]-1))]
  fun_group <- read.csv(fgfile, header = T, stringsAsFactors = FALSE)#[, c(1,3, 4,5,6, 9,16, 12)]
  fun_group <- subset(fun_group, fun_group$IsTurnedOn == 1)
  
  biolprm <- readLines(biolprm)
  
  # diet_l <- NULL
  # tot_pred <- NULL
  
  if(sum(names(fun_group) == "InvertType") > 0)
    names(fun_group)[names(fun_group) == "InvertType"] <- "GroupType"
  
  # Subset invertebrates and vertebrates
  rs_names <- str_trim(fun_group[fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD"), "Name"]) # trim white space
  rs_codes <- str_trim(fun_group[fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD"), "Code"])
  invert_names <- fun_group[!(fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD")),]
  invert_codes<- str_trim(fun_group[!(fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD")), "Code"])
  
  cat("### ------------ Creating dynamic labels for vat                         ------------ ###\n")
  max_tracer <- nc_out$nvars
  max_layers <- length(nc_out$dim$z$vals)
  max_time <- length(nc_out$dim$t$vals)
  var_names <- names(nc_out$var)
  
  cat("### ------------ Creating map from BGM file                              ------------ ###\n")
  
  # Find number of boxes
  numboxes <- length(grep("# Box number", bgm))
  
  
  cat("### ------------ Setting up disaggregated spatial plots                  ------------ ###\n")
  nums <- grep("Nums", var_names, value = TRUE)
  N <- grep("_N", var_names, value = TRUE)
  N <- N[-grep("_Nums", N, value = FALSE)]
  tot_num <- c(nums)
  
  ResNs <- grep("ResN", var_names, value = TRUE)
  StructNs <- grep("StructN", var_names, value = TRUE)
  
  # extract tracers from the ncdf4 object
  vars <- list()
  for (i in 1:length(tot_num)){
    vars[[i]] <- ncvar_get(nc = nc_out, varid = tot_num[i])
  }
  names(vars) <- tot_num
  
  
  varsResN <- list()
  for (i in 1:length(ResNs)){
    varsResN[[i]] <- ncvar_get(nc = nc_out, varid = ResNs[i])
  }
  names(varsResN) <- ResNs
  
  
  varsStrN <- list()
  for (i in 1:length(StructNs)){
    varsStrN[[i]] <- ncvar_get(nc = nc_out, varid = StructNs[i])
  }
  names(varsStrN) <- StructNs
  
  
  # Create Erla's plots
  nominal_dz <- ncvar_get(nc = nc_out, varid = "nominal_dz")
  depth_layers <- nominal_dz[,which.max(colSums(nominal_dz))]
  depth_layers <- depth_layers[-c(1, length(depth_layers))]
  depth_layers <- cumsum(rev(depth_layers))
  
  depth_labels <- rep(NA, (length(depth_layers) + 1))
  for(i in 1:(length(depth_layers) + 1)){
    if(i == 1){
      depth_labels[i] <- paste("0 - ", depth_layers[i], sep = "")
    } else if(i == (length(depth_layers) + 1))
    {
      depth_labels[i] <- paste(depth_layers[i - 1], " + ", sep ="")
    } else depth_labels[i] <- paste(depth_layers[i - 1], " - ",  depth_layers[i], sep ="")
  }
  depth_labels <- c(depth_labels, "Sediment")
  
  vert_names <- fun_group[fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD"), "Code"]
  mat_age <- grep("_age_mat", biolprm, value = T)
  species_ids <- str_split_fixed(mat_age, "_age_mat", n = 2)[,1]
  #   juvenile_age <- as.numeric(gsub("[^\\d]+", "", mat_age, perl=TRUE))[which(species_ids %in% vert_names)]
  #test how many numbers there are. Sometimes there is an example value as well as the actual value, 
  #in which case we only want the first one
  get_first_number<-function(x){
    yy<-gsub("([^\\d])","#",x,perl=TRUE)
    yyy<-unlist(str_split(yy,"#"))
    xPos<-grep("[^\\d]",yyy)[1]
    thisNum<-as.numeric(yyy[xPos])
  }
  temp<-lapply(mat_age,FUN=get_first_number)
  juvenile_age<-as.numeric(unlist(temp))
  species_ids <- species_ids[which(species_ids %in% vert_names)]
  
  erla_plots <- list()
  for(i in 1:length(species_ids)){
    spp <- fun_group[fun_group$Code == species_ids[i],c("Name", "NumCohorts")]
    spp <- str_trim(spp)
    if(juvenile_age[i] != 1){
      juv <- paste(spp[[1]], 1:(juvenile_age[i] - 1), "_Nums", sep = "")
      ad <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_Nums", sep = "")
      
      juvR <- paste(spp[[1]], 1:(juvenile_age[i] - 1), "_ResN", sep = "")
      adR <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_ResN", sep = "")
      
      juvS <- paste(spp[[1]], 1:(juvenile_age[i] - 1), "_StructN", sep = "")
      adS <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_StructN", sep = "")
      
      
      # Create the juveniles data
      juv_tmp <- NULL
      juvR_tmp <- NULL
      juvS_tmp <- NULL
      for(j in juv){
        x <- adply(vars[[j]], c(1, 3))
        juv_tmp <- rbind(juv_tmp, x)
      }
      for(j in juvS){
        x <- adply(varsStrN[[j]], c(1, 3))
        juvS_tmp <- rbind(juvS_tmp, x)
      }
      for(j in juvR){
        x <- adply(varsResN[[j]], c(1, 3))
        juvR_tmp <- rbind(juvR_tmp, x)
      }
      juv_tmp <- juv_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(juv_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(juv_tmp)-3), sep =" "))
      juv_tmp$Layer <- factor(juv_tmp$Layer,levels(juv_tmp$Layer)[c(((length(unique(juv_tmp$Layer)))-1):1, length(unique(juv_tmp$Layer)))])
      levels(juv_tmp$Layer) <- depth_labels  
      juv_tmp<-data.frame(juv_tmp) #turn to df so gather works
      juv_tmp <- gather(juv_tmp, Box, value = number, 3:ncol(juv_tmp))
      
      ### ResN mean
      juvR_tmp <- juvR_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(juvR_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(juvR_tmp)-3), sep =" "))
      juvR_tmp$Layer <- factor(juvR_tmp$Layer,levels(juvR_tmp$Layer)[c(((length(unique(juvR_tmp$Layer)))-1):1, length(unique(juvR_tmp$Layer)))])
      levels(juvR_tmp$Layer) <- depth_labels  
      juvR_tmp<-data.frame(juvR_tmp) #turn to df so gather works
      juvR_tmp <- gather(juvR_tmp, Box, value = Weight, 3:ncol(juvR_tmp))
      
      ### StructN mean
      ### ResN mean
      juvS_tmp <- juvS_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(juvS_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(juvS_tmp)-3), sep =" "))
      juvS_tmp$Layer <- factor(juvS_tmp$Layer,levels(juvS_tmp$Layer)[c(((length(unique(juvS_tmp$Layer)))-1):1, length(unique(juvS_tmp$Layer)))])
      levels(juvS_tmp$Layer) <- depth_labels  
      juvS_tmp<-data.frame(juvS_tmp) #turn to df so gather works
      juvS_tmp <- gather(juvS_tmp, Box, value = Weight, 3:ncol(juvS_tmp))
      ###
      
      ##add ResN and StructN to get TotalN
      juvTN_tmp<-juvS_tmp
      juvTN_tmp$Weight<-juvS_tmp$Weight+juvR_tmp$Weight
      
      erla_plots[[paste(spp[[1]], "Juvenile")]] <- juv_tmp
      
      erla_plots[[paste(spp[[1]], "Juvenile Weight")]] <- juvTN_tmp
      
      # Create the adults data
      ad_tmp <- NULL
      for(j in ad){
        x <- adply(vars[[j]], c(1, 3))
        ad_tmp <- rbind(ad_tmp, x)
      }
      
      adR_tmp <- NULL
      for(j in adR){
        x <- adply(varsResN[[j]], c(1, 3))
        adR_tmp <- rbind(adR_tmp, x)
      }
      
      adS_tmp <- NULL
      for(j in adS){
        x <- adply(varsStrN[[j]], c(1, 3))
        adS_tmp <- rbind(adS_tmp, x)
      }
      
      ad_tmp <- ad_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(ad_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(ad_tmp)-3), sep =" "))
      ad_tmp$Layer <- factor(ad_tmp$Layer,levels(ad_tmp$Layer)[c(((length(unique(ad_tmp$Layer)))-1):1, length(unique(ad_tmp$Layer)))])
      levels(ad_tmp$Layer) <- depth_labels
      ad_tmp<-data.frame(ad_tmp)
      ad_tmp <- gather(ad_tmp, Box, value = number, 3:ncol(ad_tmp))
      
      ###
      adR_tmp <- adR_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(adR_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(adR_tmp)-3), sep =" "))
      adR_tmp$Layer <- factor(adR_tmp$Layer,levels(adR_tmp$Layer)[c(((length(unique(adR_tmp$Layer)))-1):1, length(unique(adR_tmp$Layer)))])
      levels(adR_tmp$Layer) <- depth_labels
      adR_tmp<-data.frame(adR_tmp)
      adR_tmp <- gather(adR_tmp, Box, value = Weight, 3:ncol(adR_tmp))
      
      ###
      adS_tmp <- adS_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(adS_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(adS_tmp)-3), sep =" "))
      adS_tmp$Layer <- factor(adS_tmp$Layer,levels(adS_tmp$Layer)[c(((length(unique(adS_tmp$Layer)))-1):1, length(unique(adS_tmp$Layer)))])
      levels(adS_tmp$Layer) <- depth_labels
      adS_tmp<-data.frame(adS_tmp)
      adS_tmp <- gather(adS_tmp, Box, value = Weight, 3:ncol(adS_tmp))
      ###
      
      ##add ResN and StructN to get TotalN
      adTN_tmp<-adS_tmp
      adTN_tmp$Weight<-adS_tmp$Weight+adR_tmp$Weight
      
      erla_plots[[paste(spp[[1]], "Adult")]] <- ad_tmp
      
      erla_plots[[paste(spp[[1]], "Adult Weight")]] <- adTN_tmp
      
    } 
    else {
      ad <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_Nums", sep = "")
      ad_tmp <- NULL
      for(j in ad){
        x <- adply(vars[[j]], c(1, 3))
        ad_tmp <- rbind(ad_tmp, x)
      }
      ad_tmp <- ad_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(ad_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(ad_tmp)-3), sep =" "))
      ad_tmp$Layer <- factor(ad_tmp$Layer,levels(ad_tmp$Layer)[c(((length(unique(ad_tmp$Layer)))-1):1, length(unique(ad_tmp$Layer)))])
      levels(ad_tmp$Layer) <- depth_labels
      ad_tmp<-data.frame(ad_tmp)
      ad_tmp <- gather(ad_tmp, Box, value = number, 3:ncol(ad_tmp))
      erla_plots[[paste(spp[[1]], "Adult")]] <- ad_tmp
    }
  }
  
  # --- End Erla Plots -- #
  
  #for each vertebrate, get the relative mean weight by adult/juvenile
  relWeight <- list()
  for(i in 1:length(species_ids)){
    spp <- fun_group[fun_group$Code == species_ids[i],c("Name", "NumCohorts")]
    spp <- str_trim(spp)
    if(juvenile_age[i] != 1){
      adWeight<-erla_plots[[paste(spp[[1]], "Adult Weight")]]
      juvWeight<-erla_plots[[paste(spp[[1]], "Juvenile Weight")]]
      xx<-tapply(adWeight$Weight,adWeight$Time,mean)
      adRelWeight<-xx/xx[1]
      names(adRelWeight)<-startyear + as.double(names(adRelWeight))-1
      
      xx<-tapply(juvWeight$Weight,juvWeight$Time,mean)
      juvRelWeight<-xx/xx[1]
      names(juvRelWeight)<-startyear + as.double(names(juvRelWeight))-1
      
      relWeight[[paste(spp[[1]],"Adult")]] <- adRelWeight
      relWeight[[paste(spp[[1]],"Juvenile")]] <- juvRelWeight
    }
  }
  
  
  
  
  output <- list(disagg = vars, var_names = tot_num, max_layers = max_layers, max_time = max_time,  rs_names = rs_names, islands = islands, fun_group = fun_group, invert_names = invert_names,erla_plots = erla_plots, toutinc = toutinc, startyear = startyear, "relWeight"=relWeight)
  cat("### ------------ vat object created, you can now run the vat application ------------ ###\n") 
  return(output)
  class(output) <- "vadt"
}