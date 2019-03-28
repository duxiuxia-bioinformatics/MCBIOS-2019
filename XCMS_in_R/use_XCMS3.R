# UAB Metabolomics Workshop 2018
# Prepared by Xiuxia Du
# July 2018

# Information in this file is mostly taken from:
# https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html

rm(list=ls())
graphics.off()

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("xcms", version = "3.8")


# --------------------------------------------
# get the libraries
# --------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("xcms")
library(xcms)

biocLite("faahKO")
library(faahKO)

library(RColorBrewer)
install.packages("pander")
library(pander)
library(magrittr)

# --------------------------------------------
# specify raw data folders and working directory
# --------------------------------------------
## Get the full path to the CDF files
all_raw_files <- dir(system.file("cdf", package = "faahKO"), 
                     full.names = TRUE, 
                     recursive = TRUE)

# # For your own raw data files, specify the folder
# path_to_raw_data <- "/Users/xdu4/Documents/Duxiuxia/Analysis/UAB_workshop_2018/raw_data/faahKO"
# all_raw_files <- list.files(path_to_raw_data, recursive=T, full.names=T)

# set up the working directory
setwd("/Users/xdu4/Documents/Duxiuxia/Analysis/UAB_workshop_2018/XCMS_analysis/")

# --------------------------------------------
# import raw data
# --------------------------------------------
# get file path and names

## Create a phenodata data.frame
pd <- data.frame(sample_name=sub(basename(all_raw_files), 
                                 pattern = ".CDF", 
                                 replacement = "", 
                                 fixed = TRUE),
                 sample_group = c(rep("KO", 6), rep("WT", 6)),
                 stringsAsFactors = FALSE) 

# read raw data
raw_data <- readMSData(files=all_raw_files, 
                       pdata=new("NAnnotatedDataFrame", pd), 
                       mode="onDisk")

# --------------------------------------------
# inspect the raw data
# --------------------------------------------
# *** examine raw_data structure
str(raw_data, max.level=2)
str(raw_data@phenoData)

str(raw_data@featureData, max.level=3)
II <- raw_data@featureData@data$fileIdx

str(raw_data@experimentData, max.level=2)
str(raw_data@protocolData)
str(raw_data@processingData)

# *** plot the base peak chromatograms.
BPCs <- chromatogram(raw_data, aggregationFun="max")

# define colors for the two groups
group_colors <- brewer.pal(3, "Set1")[1:2]
names(group_colors) <- c("KO", "WT")

# plot all BPCs
par(mfrow = c(1, 1))
plot(BPCs, col=group_colors[raw_data$sample_group],
     main="Base Peak Chromatograms")

BPC_1 <- BPCs[[1]]
plot(BPC_1, type='l', col='red')

head(rtime(BPC_1))
head(intensity(BPC_1))

# *** plot the TIC
TICs <- chromatogram(raw_data, aggregationFun="sum")

# plot all TICs
plot(TICs, col=group_colors[raw_data$sample_group],
     main="Total Ion Chromatograms")

TIC_1 <- TICs[[1]]
plot(TIC_1, type='l', col='red')

# or individual TICs can be obtained by
TIC_by_file <- split(tic(raw_data), 
                     f=fromFile(raw_data))

plot(TIC_by_file[[1]], type='l', col='red',
     xlab='retention time index',
     ylab='total intensity',
     main='One TIC')

# *** an alternative way to plot the TIC by directly accessing the raw_data
fileIdx_vec <- raw_data@featureData@data$fileIdx
file_index <- 1
II <- which(fileIdx_vec==file_index)
par(mfrow = c(1, 1))

figure_title <- paste("TIC: file=", file_index, sep="")

plot(raw_data@featureData@data$retentionTime[II], 
     raw_data@featureData@data$totIonCurrent[II], 
     type='l', 
     xlab="RT (sec)", ylab="total intensity",
     main=figure_title)

# *** extract retention time
# use the method rtime() to extract retention time
rtime(raw_data)
# returns a numeric vector of all the retention time

# *** extract mz
mzs <- mz(raw_data)
# mzs is a list

mzs_by_file <- split(mzs, f=fromFile(raw_data))
# don't display mzs_by_file

length(mzs_by_file)

# mz values in spectrum 1 in file 2
mzs_by_file[[2]][1]

# *** extract intensity
intensity_values <- intensity(raw_data)
intensity_by_file <- split(intensity_values, 
                           f=fromFile(raw_data))

# *** plot a spectrum
file_index <- 3
spectrum_index <- 250

figure_title <- paste("file=", 
                      file_index, 
                      ", scan=", 
                      spectrum_index, 
                      sep="")

plot(mzs_by_file[[file_index]][[spectrum_index]], 
     intensity_by_file[[file_index]][[spectrum_index]],
     type='h', 
     xlab='m/z', ylab='intensity', 
     main=figure_title)

# ---------------------------------------------
# chromatographic peak detection
# ---------------------------------------------
# centWave algorithm: The two most critical parameters for centWave are peakWidth and ppm.
# peakWidth: expected range of chromatographic peak widths
# ppm: maximum expected deviation of m/z values of centroids corresponding to one chromatographic peak

# *** To determine peakWidth, we plot the EIC for one peak.
# define the rt and m/z range of the peak area
rt_range <- c(2700, 2900)
mz_range <- c(334.9, 335.1)

# extract the chromatogram
chr_raw <- chromatogram(raw_data, mz=mz_range, rt=rt_range)
plot(chr_raw, col=group_colors[chr_raw$sample_group])
# the peak width is about 50 seconds

# we will set peakWidth that will contain the observed peak width above

# *** to determine the ppm parameter
raw_data %>%
    filterRt(rt=rt_range) %>%
    filterMz(mz=mz_range) %>%
    plot(type='XIC')

# *** specify centWave parameters
cwt_parameters <- CentWaveParam(peakwidth=c(30, 80),
                                ppm=50,
                                noise=1000)
# type ?CentWaveParam to understand other parameters

xdata <- findChromPeaks(raw_data, 
                        param=cwt_parameters)

# *** access the peak detection results
dim(chromPeaks(xdata))
head(chromPeaks(xdata))

# *** save the peak detection results
file_name <- "peak_detection_results.csv"
write.csv(file=file_name, chromPeaks(xdata))

# *** plot the detected peaks
par(mfrow = c(1, 1))
plotChromPeaks(xdata, file=12)

# ----------------------------------------------
# alignment
# ----------------------------------------------
# *** obiwarp alignment
obiwarp_parameters <- ObiwarpParam(binSize=0.6)
?ObiwarpParam

xdata <- adjustRtime(xdata, param=obiwarp_parameters)
# The obiwarp alignment is performed directly on the profile-matrix and is therefore independent of peak detection. 

# extract adjusted retention times
head(adjustedRtime(xdata))

# or simply use the rtime method
head(rtime(xdata))

# raw retention time can be extracted from an XCMSnExp containing aligned data with
rtime(xdata, adjusted=FALSE)

# to evaluate the impact of the alignment
BPCs_adjusted <- chromatogram(xdata, aggregationFun="max")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(BPCs_adjusted, col=group_colors[BPCs_adjusted$sample_group])
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group]) 
# too large differences between adjusted and raw retention times could indicate poorly performing samples or alignment.

# *** peak groups alignment
# This algorithm adjusts retention time by aligning previously identified hook peaks (chromatographic peaks present in most samples).
# Next, we drop the adjusted retention time and re-align.

# use the DROP method
# does the object have adjusted retention time
hasAdjustedRtime(xdata)

# drop the alignment results
xdata <- dropAdjustedRtime(xdata)

hasAdjustedRtime(xdata)

# correspondence: group peaks across samples
peak_density_parameters <- PeakDensityParam(sampleGroups=xdata$sample_group,
                                            minFraction = 0.8)
# minFraction: the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature).

xdata <- groupChromPeaks(xdata, param=peak_density_parameters)

# now the retention time correction based on house keeping peak groups that are present in most samples. First, the retention time deviation of these peaks groups is described by fitting either a polymonial or a linear model to the data points. These models are subsequently used to adjust the retention time.
peak_group_parameters <- PeakGroupsParam(minFraction = 0.85)
?PeakGroupsParam

# get the peak groups that would be used for alignment
xdata <- adjustRtime(xdata, param=peak_group_parameters)

## Plot the difference of adjusted to raw retention time.
par(mfrow=c(1,1))
plotAdjustedRtime(xdata, 
                  col = group_colors[xdata$sample_group],
                  peakGroupsCol = "grey", 
                  peakGroupsPch = 1)

# evaluate the impact of the alignmnet
par(mfrow = c(2, 1))

## Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, 
                        rt = rt_range, 
                        mz = mz_range)
plot(chr_adj, col = group_colors[chr_raw$sample_group]) 

# --------------------------------------------
# correspondence
# --------------------------------------------
# purpose: to match detected chromatographic peaks between samples

## Define the mz slice.
mz_range <- c(305.05, 305.15)
rt_range <- c(2500, 4000)

## Extract and plot the chromatograms
chr_mzr <- chromatogram(xdata, 
                        mz = mz_range, 
                        rt = rt_range)

par(mfrow = c(3, 1), mar = c(1, 4, 1, 0.5))

cols <- group_colors[chr_mzr$sample_group]
plot(chr_mzr, col = cols, 
     xaxt = "n", xlab = "")

## Highlight the detected peaks in that region.
highlightChromPeaks(xdata, 
                    mz = mz_range, 
                    col = cols, 
                    type = "point", 
                    pch = 16)

## Define the parameters for the peak density method
peak_density_parameters <- PeakDensityParam(sampleGroups = xdata$sample_group,
                                            minFraction = 0.4, 
                                            bw = 30)
# minFraction: minimum fraction of samples in at least one sample group in which hte peaks have to be present to be considered as a peak group
# bw: standard deviation of the smoothing kernel

par(mar = c(4, 4, 1, 0.5))

plotChromPeakDensity(xdata, 
                     mz = mz_range, 
                     col = cols, 
                     param = peak_density_parameters, 
                     pch = 16, 
                     xlim = c(2500, 4000))

## Use a different bw
peak_density_parameters <- PeakDensityParam(sampleGroups = xdata$sample_group,
                                            minFraction = 0.4, 
                                            bw = 20)

plotChromPeakDensity(xdata, 
                     mz = mz_range, 
                     col = cols, 
                     param = peak_density_parameters, 
                     pch = 16, 
                     xlim = c(2500, 4000)) 

## Perform the correspondence
peak_density_parameters <- PeakDensityParam(sampleGroups = xdata$sample_group,
                                            minFraction = 0.4, 
                                            bw = 20)
xdata <- groupChromPeaks(xdata, 
                         param = peak_density_parameters) 

## Extract the feature definitions
featureDefinitions(xdata)

## Extract the into column for each feature.
file_name <- "peaks_after_correspondence.csv"
write.csv(file=file_name, featureValues(xdata, value="into"))

head(featureValues(xdata, value = "into"))

# -------------------------------------------------
# fill missing peaks
# -------------------------------------------------
# defines intensity values of missing features by integrating the signal in the mz-rt region of the feature.

## Filling missing peaks using default settings. Alternatively we could
## pass a FillChromPeaksParam object to the method.
xdata <- fillChromPeaks(xdata)

?fillChromPeaks

file_name <- "peaks_after_filling_missing.csv"
write.csv(file=file_name, featureValues(xdata))

head(featureValues(xdata))

## Missing values before filling in peaks
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))

## Missing values after filling in peaks
apply(featureValues(xdata), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))

## Extract the features and log2 transform them
ft_ints <- log2(featureValues(xdata, value = "into"))

par(mfrow=c(1,1))
pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)

## Plot the PCA
cols <- group_colors[xdata$sample_group]
pcSummary <- summary(pc)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "", xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100, digits = 3), " % variance"), ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100, digits = 3), " % variance"), col = "darkgrey", bg = cols, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2], labels = xdata$sample_name, col = "darkgrey", pos = 3, cex = 2)

# ------------------------------------
# evaluate the process history
# ------------------------------------
processHistory(xdata)

ph <- processHistory(xdata, type = "Retention time correction")

ph 

