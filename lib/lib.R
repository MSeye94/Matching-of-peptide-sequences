# options(repos = BiocManager::repositories())

usePackage <- function(p)
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, force = TRUE)
  require(p, character.only = TRUE)
}



if(!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  if(get_r_version()>"4.1.3") {
    message(paste0("R version ",get_r_version()," is detected \n",
                   "Installing package BiocManager..."))
    BiocManager::install()
  } else if(get_r_version()=="4.1.3") {
    message(paste0("R version ",get_r_version()," is detected \n",
                   "Installing package BiocManager 3.14..."))
    BiocManager::install(version = "3.14", update = FALSE)
  } else {
    stop("Pleas Install R version >= 4.1.3")
  }
}





if(!require(xcms))
  devtools::install_local("lib/packagesR/xcms")
require(xcms)

usePackage("shiny")
usePackage("shinyjs")
usePackage("shinythemes")
usePackage("shinyWidgets")
usePackage("shinycssloaders")
usePackage("shinyalert")
usePackage("shinydashboard")
usePackage("shinyFiles")
usePackage("progressr")
usePackage("devtools")

usePackage("np")
usePackage("MASS")
usePackage("DT")


### Volume for shinyFiles


### required packages server
initPackages <- function(session){
  
  pkgsBiocManager <- c("RCurl","tools","RSQLite","colourpicker", "yaml","utils", "BiocParallel")
  
  pkgsCRAN<-c("stringr", "bigstatsr", "bigreadr", "readxl", "signal","ggplot2","plotly",
              "tidyverse", "dplyr", "ggalt", "parallel", "grid", "gridExtra", "cowplot", "patchwork", "MsCoreUtils")
  
  allPkgs<-c(pkgsBiocManager, pkgsCRAN)
  
  message("Installing packages of Biconductor...\n")
  if(!all(pkgsBiocManager %in% installed.packages()[,1])){
    withProgress(message = 'Installing packages of Biconductor...', value = 0, {
      for (i in 1:length(pkgsBiocManager)) {
        incProgress(1/length(pkgsBiocManager), detail = paste("Installing ",pkgsBiocManager[i],collapse=""))
        if(!is.element(pkgsBiocManager[i], installed.packages()[,1])){
          BiocManager::install(pkgsBiocManager[i])
        }
      }
    })
  }
  
  
  message("Installing packages of CRAN...\n")
  if(!all( pkgsCRAN %in% installed.packages()[,1])){
    withProgress(message = 'Installing packages of CRAN...', value = 0, {
      for (i in 1:length(pkgsCRAN)) {
        incProgress(1/length(pkgsCRAN), detail = paste("Installing ",pkgsCRAN[i],collapse=""))
        if(!is.element(pkgsCRAN[i], installed.packages()[,1])){
          install.packages(pkgsCRAN[i],dep = TRUE, force = TRUE)
          
        }
      }
    })
  }
  
  
  message("Installing packages finished Ok!\n")
  message("Loading packages...\n")
  withProgress(message = 'Loading packages...', value = 0, {
    for (i in 1:length(allPkgs)) {
      incProgress(1/length(allPkgs), detail = paste("Loading ",allPkgs[i],collapse=""))
      require(allPkgs[i],character.only = TRUE)
    }
  })
  
}







volumes <- c(Home = fs::path_home(), getVolumes()())








###### Fonction R pour matcher deux jeux de données

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("MsCoreUtils")

# library(MsCoreUtils)
# library(xcms)

# This function allows to match two mass lists containing two columns "mz" and "rt" or specify by mzcol and rtcol.
# It searches for correspondences of the elements of x in the elements of the array table.
# It returns a list containing the indices of the rows in table which match well with the rows of x,
# the percentage of the values which match well
# and the two tables with the rows which match and the rows which do not match its replace by NA in table


matchRtMz <- function(x,
                      table,
                      nomatch = NA_integer_,
                      rt_tolerance = 2,
                      mzcol = "mz",
                      rtcol = "rt",
                      session) {
  
  if(!require(MsCoreUtils))
    stop("R package \"MsCoreUtils\" is not found, please install R package \"MsCoreUtils\" ")
  
  
  if (!is.numeric(nomatch) || length(nomatch) != 1L)
    stop("'nomatch' has to be a 'numeric' of length one.")
  
  if (length(dim(x)) != 2 || length(dim(table)) != 2)
    stop("'x' and 'table' have to be two data frames")
  if (!all(c(mzcol, rtcol) %in% colnames(x)) ||
      !all(c(mzcol, rtcol) %in% colnames(table))){
    
    # stop("Columns '", mzcol, "' and '", rtcol, "' not found in 'x' and ",
    #      "'table'")
    
    sendSweetAlert(
      session = session,
      title = "Warning !",
      text = paste("Required columns :", "'mz', 'rt' not found at the same time in reference file and file to align.", sep = " "),
      type = "warning",
      width = "60%"
    )
    
    return(NULL)
  } else{
    
    if(nrow(x)<2)
      stop("'x' must have two or more rows!")
    
    
    table<- table[order(table[,mzcol]), ]
    mz1 <- x[, mzcol]
    rt1 <- x[, rtcol]
    mz2 <- table[, mzcol]
    rt2 <- table[, rtcol]
    idxl <- vector("list", length = nrow(x))
    
    
    withProgress(message = 'Matching sample', value = 0, {
      
      pb_match <- txtProgressBar(min=1, max = length(seq_along(idxl)), style = 3)
      cat("matching rt1~rt2 and mz1~mz2 .....!\n")
      
      
      for (i in seq_along(idxl)) {
        setTxtProgressBar(pb_match, i)
        
        incProgress(1/length(seq_along(idxl)), detail = paste(round(i*100/length(seq_along(idxl)),0)," %..."))
        
        matches <- which(abs(rt2 - rt1[i]) <= rt_tolerance)
        
        
        if (length(matches)) {
          
          #matches <- matches[which(abs(mz2[matches]-mz1[i])<=MsCoreUtils::ppm(mz1[i], 100))]
          
          if(mz1[i]<4000){
            matches <- matches[which(abs(mz2[matches]-mz1[i])<=MsCoreUtils::ppm(mz1[i], 50))]
          } else if(mz1[i]>=4000 & mz1[i]<=6000){
            ##fite a linear model between 4KDA et 6KDA to get ppm (50-150)
            Data.model<-data.frame(KDA = seq(from = 4, to = 6, by=0.1)*1000, ppm.data = seq(from = 50, to = 150, by=5))
            model.ppm<-lm(ppm.data~KDA, data = Data.model)
            ppm.pred<-as.double(predict.lm(model.ppm, newdata = data.frame(KDA=mz1[i])))
            matches <- matches[which(abs(mz2[matches]-mz1[i])<=MsCoreUtils::ppm(mz1[i],ppm.pred))]
          } else {
            matches <- matches[which(abs(mz2[matches]-mz1[i])<=MsCoreUtils::ppm(mz1[i], 150))]
          }
          
          
          if (length(matches)) {
            matches <- matches[abs(mz2[matches]-mz1[i])==min(abs(mz2[matches]-mz1[i]))][1]
            idxl[[i]] <- matches
          }
          else idxl[[i]] <- nomatch
        } else idxl[[i]] <- nomatch
      }
      close(pb_match)
      
    })
    
    
    cat("OK...!\n")
    rownames(table)<-1:nrow(table)
    # x<-x[,c("mz","rt")]
    # table<-table[,c("mz","rt")]
    colnames_x<-colnames(x)
    colnames_table<-colnames(table)
    colnames(x)[which(colnames_x %in% c("mz","rt"))]<-c("mz1","rt1")
    colnames(table)[colnames_table %in% c("mz","rt")]<-c("mz2", "rt2")
    MatchTable = cbind(x[seq_along(x[,"mz1"]),],table[as.integer(idxl),])
    
    rownames(MatchTable)<-rownames(x)
    
    percentageMatch = round((nrow(MatchTable)-sum(is.na(MatchTable$mz2)))*100/nrow(MatchTable), digits = 2)
    numberMatch = (nrow(MatchTable)-sum(is.na(MatchTable$mz2)))
    # 
    #   View(MatchTable)
    #   write.table(MatchTable, file = "MatchTable.csv")
    #percentageMatch = round((nrow(table)-sum(is.na(MatchTable$mz2)))*100/nrow(table), digits = 2)
    return(ResMatch = list(idxl = as.integer(idxl),
                           percentageMatch = percentageMatch,
                           numberMatch = numberMatch,
                           MatchTable = MatchTable))
    
  }
    
 
}






theme_ben <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # L'ensemble de la figure
      plot.title = element_text(size = rel(0.85), face = "bold", color = "blue", margin = margin(0,0,5,0), hjust = 0.5),
      # Zone où se situe le graphique
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      # Les axes
      axis.title = element_text(size = rel(0.85), face = "bold.italic"),
      axis.text = element_text(size = rel(0.70), face = "bold.italic"),
      #axis.text.y = element_text(margin = margin(r = 0.8 *(100/2) / 2), hjust = 1),
      axis.line = element_line(color = "black", arrow = arrow(length = unit(0.3, "lines"), type = "closed")),
      # La légende
      legend.title = element_text(size = rel(0.85), face = "bold.italic", hjust = 0.5),
      legend.text = element_text(size = rel(0.70), face = "bold.italic"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(0.5, "cm"),
      legend.key.width = unit(0.5,"cm"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      # Les étiquettes dans le cas d'un facetting
      strip.background = element_rect(fill = "#17252D", color = "#17252D"),
      strip.text = element_text(size = rel(1), face = "bold.italic", color = "white", margin = margin(5,0,5,0))
    )
}





#################### Regroupement des pics similaire en temps et en masse dans et entre les runs #############################




#' Low level function to group chromatographic peaks within a m/z slice.
#'
#' @param x `matrix` such as the one returned by `chromPeaks,XCMSnExp`, just
#'     with the peaks within one m/z slice. Note that we require in addition
#'     a column `"index"` with the index of the peak within the full peak table.
#'
#' @param return `data.frame`
#'
#' @author Johannes Rainer
#'
#' @noRd
group_peaks_density_new <- function(x, bw, densFrom, densTo, densN, sampleGroups,
                                 sampleGroupTable, minFraction,
                                 minSamples, maxFeatures, sleep = 0) {
  den <- density(x[, "rt"], bw = bw, from = densFrom, to = densTo,
                 n = densN)
  maxden <- max(den$y)
  deny <- den$y
  sampleGroupNames <- names(sampleGroupTable)
  nSampleGroups <- length(sampleGroupNames)
  col_nms <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
               "npeaks", sampleGroupNames)
  res_mat <- matrix(nrow = 0, ncol = length(col_nms),
                    dimnames = list(character(), col_nms))
  res_idx <- list()
  while (deny[maxy <- which.max(deny)] > maxden / 20 && nrow(res_mat) <
         maxFeatures) {
    grange <- xcms:::descendMin(deny, maxy)
    deny[grange[1]:grange[2]] <- 0
    gidx <- which(x[,"rt"] >= den$x[grange[1]] &
                    x[,"rt"] <= den$x[grange[2]])
    ## Determine the sample group of the samples in which the peaks
    ## were detected and check if they correspond to the required limits.
    tt <- table(sampleGroups[unique(x[gidx, "sample"])])
    if (!any(tt / sampleGroupTable[names(tt)] >= minFraction &
             tt >= minSamples))
      next
    gcount <- rep(0, length(sampleGroupNames))
    names(gcount) <- sampleGroupNames
    gcount[names(tt)] <- as.numeric(tt)
    res_mat <- rbind(res_mat,
                     c(median(x[gidx, "mz"]),
                       range(x[gidx, "mz"]),
                       median(x[gidx, "rt"]),
                       range(x[gidx, "rt"]),
                       length(gidx),
                       gcount)
    )
    res_idx <- c(res_idx, list(unname(sort(x[gidx, "index"]))))
  }
  if (sleep > 0) {
    ## Plot the density
    plot(den, main = paste(round(min(x[,"mz"]), 2), "-",
                           round(max(x[,"mz"]), 2)))
    ## Highlight peaks per sample group.
    for (j in seq_len(nSampleGroups)) {
      ## Which peaks belong to this sample group.
      cur_group_samples <- which(sampleGroups == sampleGroupNames[j])
      idx <- x[, "sample"] %in% cur_group_samples
      points(x[idx, "rt"], x[idx, "into"] /
               max(x[, "into"]) * maxden,
             col = j, pch=20)
    }
    for (j in seq_len(nrow(res_mat)))
      abline(v = res_mat[j, 5:6], lty = "dashed", col = j)
    Sys.sleep(sleep)
  }
  res <- as.data.frame(res_mat)
  res$peakidx <- res_idx
  res
}


## Correspondence functions.
#' @include functions-Params.R

#' @title Core API function for peak density based chromatographic peak
#' grouping
#'
#' @description
#'
#' The `do_groupChromPeaks_density` function performs chromatographic peak
#' grouping based on the density (distribution) of peaks, found in different
#' samples, along the retention time axis in slices of overlapping mz ranges.
#'
#' @details For overlapping slices along the mz dimension, the function
#' calculates the density distribution of identified peaks along the
#' retention time axis and groups peaks from the same or different samples
#' that are close to each other. See (Smith 2006) for more details.
#'
#' @note The default settings might not be appropriate for all LC/GC-MS setups,
#' especially the `bw` and `binSize` parameter should be adjusted
#' accordingly.
#'
#' @param peaks A `matrix` or `data.frame` with the mz values and
#' retention times of the identified chromatographic peaks in all samples of an
#' experiment. Required columns are `"mz"`, `"rt"` and
#' `"sample"`. The latter should contain `numeric` values representing
#' the index of the sample in which the peak was found.
#'
#' @inheritParams groupChromPeaks-density
#'
#' @param sleep `numeric(1)` defining the time to *sleep* between
#'     iterations and plot the result from the current iteration.
#'
#' @return
#'
#' A `data.frame`, each row representing a (mz-rt) feature (i.e. a peak group)
#' with columns:
#'
#' - `"mzmed"`: median of the peaks' apex mz values.
#' - `"mzmin"`: smallest mz value of all peaks' apex within the feature.
#' - `"mzmax"`:largest mz value of all peaks' apex within the feature.
#' - `"rtmed"`: the median of the peaks' retention times.
#' - `"rtmin"`: the smallest retention time of the peaks in the group.
#' - `"rtmax"`: the largest retention time of the peaks in the group.
#' - `"npeaks"`: the total number of peaks assigned to the feature.
#' - `"peakidx"`: a `list` with the indices of all peaks in a feature in the
#'   `peaks` input matrix.
#'
#' Note that this number can be larger than the total number of samples, since
#' multiple peaks from the same sample could be assigned to a feature.
#'
#' @references
#'
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' Anal. Chem. 2006, 78:779-787.
#'
#' @author Colin Smith, Johannes Rainer
#'
#' @family core peak grouping algorithms
#'
#' @md
#'
#' @examples
#' ## Load the test file
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Extract the matrix with the identified peaks from the xcmsSet:
#' pks <- chromPeaks(faahko_sub)
#'
#' ## Perform the peak grouping with default settings:
#' res <- do_groupChromPeaks_density(pks, sampleGroups = rep(1, 3))
#'
#' ## The feature definitions:
#' head(res)
do_groupChromPeaks_density_new <- function(peaks, sampleGroups,
                                       bw = 30, minFraction = 0.5, minSamples = 1,
                                       ppm = TRUE,
                                       binSize = 0.25, maxFeatures = 50,
                                       sleep = 0) {

  if (missing(sampleGroups))
    stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
         "length equal to the number of samples specifying the group ",
         "assignment of the samples.")
  if (missing(peaks))
    stop("Parameter 'peaks' is missing!")
  if (!(is.matrix(peaks) | is.data.frame(peaks)))
    stop("'peaks' has to be a 'matrix' or a 'data.frame'!")
  ## Check that we've got all required columns
  .reqCols <- c("mz", "rt", "sample")
  if (sleep > 0)
    .reqCols <- c(.reqCols, "into")
  if (!all(.reqCols %in% colnames(peaks)))
    stop("Required columns ",
         paste0("'", .reqCols[!.reqCols %in% colnames(peaks)],"'",
                collapse = ", "), " not found in 'peaks' parameter")
  
  sampleGroups <- as.character(sampleGroups)
  sampleGroupNames <- unique(sampleGroups)
  sampleGroupTable <- table(sampleGroups)
  nSampleGroups <- length(sampleGroupTable)
  
  
  ## Check that sample groups matches with sample column.
  if (max(peaks[, "sample"]) > length(sampleGroups))
    stop("Sample indices in 'peaks' are larger than there are sample",
         " groups specified with 'sampleGroups'!")
  
  peaks <- cbind(peaks[, .reqCols, drop = FALSE],
                 index = seq_len(nrow(peaks)))
  
  
  ## Order peaks matrix by mz
  peaks <- peaks[order(peaks[, "mz"]), , drop = FALSE]
  rownames(peaks) <- NULL
  rtRange <- range(peaks[, "rt"])
  
  
  
  # # Define the mass slices and the index in the peaks matrix with an mz
  # # value >= mass[i].
  # if(ppm>0){
  #   message("ppm tolerance is used for mass grouping...\n")
  #   if(!require(MsCoreUtils))
  #     stop("Package MsCoreUtils is required !")
  #   mass<-peaks[1, "mz"]
  #   m<-peaks[1, "mz"]
  #   while (m<=peaks[nrow(peaks), "mz"]+ MsCoreUtils::ppm(peaks[nrow(peaks), "mz"], ppm)) {
  #     m<-m+MsCoreUtils::ppm(m,ppm = ppm)
  #     mass<-c(mass,m)
  #   }
  # } else {
  #   message("mz or m tolerance is used for mass grouping...\n")
  #   mass <- seq(peaks[1, "mz"], peaks[nrow(peaks), "mz"] + binSize,
  #               by = binSize/2)
  # }
  
  #masspos <- findEqualGreaterM(peaks[, "mz"], mass)
  
  if(!require(MsCoreUtils))
    stop("Package MsCoreUtils is required !")
  
  if(ppm){
    message("\n ppm tolerance is used for mass grouping...")
    
    idxl_matches<- vector("list", length = nrow(peaks))
    for(i in seq_along(idxl_matches)){
      if(peaks[i, "mz"]<4000){
        matches <- which(abs(peaks[, "mz"]-peaks[i, "mz"])<=MsCoreUtils::ppm(peaks[i, "mz"], 50))
      } else if(peaks[i, "mz"]>=4000 & peaks[i, "mz"]<=6000){
        ##fite a linear model between 4KDA et 6KDA to get ppm (50-150)
        Data.model<-data.frame(KDA = seq(from = 4, to = 6, by=0.1)*1000, ppm.data = seq(from = 50, to = 150, by=5))
        model.ppm<-lm(ppm.data~KDA, data = Data.model)
        ppm.pred<-as.double(predict.lm(model.ppm, newdata = data.frame(KDA=peaks[i, "mz"])))
        
        matches <- which(abs(peaks[, "mz"]-peaks[i, "mz"])<=MsCoreUtils::ppm(peaks[i, "mz"], ppm.pred))
      } else {
        matches <- which(abs(peaks[, "mz"]-peaks[i, "mz"])<=MsCoreUtils::ppm(peaks[i, "mz"], 150))
      }
      
      
      idxl_matches[[i]] <- matches
    }
  } else {
    message("\n mz or m tolerance is used for mass grouping...")
    idxl_matches<- vector("list", length = nrow(peaks))
    for(i in seq_along(idxl_matches)){
      matches <- which(abs(peaks[, "mz"]-peaks[i, "mz"])<=binSize)
      
      idxl_matches[[i]] <- matches
    }
  }
  
  
  
  
  densFrom <- rtRange[1] - 3 * bw
  densTo <- rtRange[2] + 3 * bw
  ## Increase the number of sampling points for the density distribution.
  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange) / (bw / 2)))))
  endIdx <- 0
  message("\n Processing ", nrow(peaks), " peaks... ",
          appendLF = FALSE)
  
  resL <- vector("list", nrow(peaks))
  for (i in seq_along(idxl_matches)) {
    ## That's identifying overlapping mz slices.
    # startIdx <- masspos[i]
    # endIdx <- masspos[i + 2] - 1
    # if (endIdx - startIdx < 0)
    #     next
    resL[[i]] <- group_peaks_density_new(#peaks[startIdx:endIdx, , drop = FALSE],
      peaks[idxl_matches[[i]], , drop = FALSE],
      bw = bw, densFrom = densFrom,
      densTo = densTo, densN = densN,
      sampleGroups = sampleGroups,
      sampleGroupTable = sampleGroupTable,
      minFraction = minFraction,
      minSamples = minSamples,
      maxFeatures = maxFeatures,
      sleep = sleep)
  }
  message("OK")
  res <- do.call(rbind, resL)
  
  if (nrow(res)) {
    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(
      as.matrix(res[, (match("npeaks", colnames(res)) +1):(ncol(res) -1),
                    drop = FALSE]))
    uorder <- order(-numsamp, res[, "npeaks"])
    
    uindex <- xcms:::rectUnique(
      as.matrix(res[, c("mzmin", "mzmax", "rtmin", "rtmax"),
                    drop = FALSE]), uorder)
    res <- res[uindex, , drop = FALSE]
    rownames(res) <- NULL
  }
  res
}















groupPing<-function(peaks, 
                    sampleGroups,
                    bw = 30, 
                    minFraction = 0.5, 
                    minSamples = 1,
                    ppm = TRUE,
                    binSize = 0.25, 
                    maxFeatures = 50,
                    sleep = 0)  {
  

  
  data_group<-do_groupChromPeaks_density_new(peaks = peaks, 
                                             sampleGroups = sampleGroups,
                                             bw = bw, 
                                             minFraction = minFraction, 
                                             minSamples = minSamples,
                                             ppm = ppm,
                                             binSize = binSize, 
                                             maxFeatures = maxFeatures,
                                             sleep = sleep) 
  return(data_group)
}


processing_Grouping_xcms<-function(peaks_sample_list,
                                   pheno_Data,
                                   bw = 30,
                                   minFraction = 0.1,
                                   minSamples = 1,
                                   ppm = TRUE,
                                   binSize = 0.25,
                                   maxFeatures = 50){

  
  message("\n --- Peak grouping within and accross sample  ---\n")

  
  withProgress(message = 'Peak grouping', value = 0, {
    
    massif_list_samples<-list()
    peptide_abondances_list<-list()
    peptide_list<-list()
    n<-length(pheno_Data$Filename)
    for (i in 1:n) {
      
      incProgress(1/n, detail = paste(round(i*100/n,0)," %..."))
      
      system.time(massif_list_samples[[i]]<-groupPing(peaks = peaks_sample_list[[i]],
                                                      sampleGroups = pheno_Data$Class[i],
                                                      bw = bw, 
                                                      minFraction = minFraction, 
                                                      minSamples = minSamples,
                                                      ppm = ppm,
                                                      binSize = binSize, 
                                                      maxFeatures = maxFeatures)) 
      
      
     
      
      
      
      
       
       


        peptide_abondances_list[[i]]<-xcms:::.feature_values(pks = peaks_sample_list[[i]],
                                                             fts = massif_list_samples[[i]],
                                                             colnames = pheno_Data$Filename[i],
                                                             value = "maxo",
                                                             method = "sum")
       


        peptide_list[[i]]<-data.frame(mz = massif_list_samples[[i]]$mzmed,
                                                 rt = massif_list_samples[[i]]$rtmed,
                                                 intensity = as.vector(peptide_abondances_list[[i]]),
                                                 sample = pheno_Data$Filename[i])
        }
    
    
  })
  
 
 names(peptide_list)<-pheno_Data$Filename

  
  message("\n --- Peak grouping within and accross sample  ---\n")

  cat("--- End peak grouping ...!\n")
  return(peptide_list)
}

