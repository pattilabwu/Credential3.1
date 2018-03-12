# Sample Script for XCMS data processing and credentialing

# Step0: Installation of pre-requisite R packages and setting working directory
source("https://bioconductor.org/biocLite.R")
biocLite("xcms")
install.packages("data.table")
install.packages("devtools")
install.pacakges("utils")
devtools::install_github("pattilabwu/Credential3.1")

# load required packages
library(xcms)
library(Credential3.1)
library(data.table)
library(utils)

# Set up working directory
setwd("/Users/Lingjue/Desktop/CredentialingDemo") # Mac OS format
setwd("C:/Users/Mike/Desktop/CredentialingDemo") # Windows format


# Step 1: Data-preprocessing with XCMS
    # suppose your mzXML files are stored in "1T1" and "1T2" folder under the working directory.

# For "1T1" condition

  # 1. peak detection with centWave algorithm
  xs_1 = xcmsSet("./1T1", method="centWave", ppm=20, peakwidth = c(10,30), snthresh=5, prefilter=c(3,100))
  # 2. initial peak grouping, 1st round
  xs_1= group(xs_1, bw=5, mzwid=.015, minfrac=0.5)
  # 3. retention time correction
  xs_1 = retcor(xs_1, method="obiwarp",profStep=1)
  # 4. feature grouping, 2nd round
  xs_1 = group(xs_1, bw=5, mzwid=.015, minfrac=0.5)
  # 5. filling missing peaks
  xs_1 = fillPeaks(xs_1)

# For "1T2" condition
  xs_2 = xcmsSet("./1T2", method="centWave", ppm=20, peakwidth=c(5,30), snthresh=5, prefilter=c(3,1000))
  xs_2 = group(xs_2, bw=5, mzwid=.015, minfrac=0.5)
  xs_2 = retcor(xs_2, method="obiwarp",profStep=1)
  xs_2 = group(xs_2, bw=5, mzwid=.015, minfrac=0.5)
  xs_2 = fillPeaks(xs_2)


# Step 2_a: Automated Processing of XCMS files by Credentialing

  featureGroup1 = xsGroupExtract(xs_1,intchoice = "into", sampling = 1)
  featureGroup2 = xsGroupExtract(xs_2,intchoice = "into", sampling = 1)
  features1T1 = data.table(featureGroup1$credTable)		# data.table object with 4 columns: "cc", "rt", "mz", "i"
  features1T2 = data.table(featureGroup2&credTable)		# data.table object with 4 columns: "cc", "rt", "mz", "i"

# Or Step 2_b: Manual Processing of Files from Other Data Processing Pipelines
  # The features tables should include only four columns: “cc”, “mz”, “rt” and “i”.

  features1T1 = data.table(read.csv("features1T1.csv"))
  features1T2 = data.table(read.csv("features1T2.csv"))

# Step 3: Credentialing features

  credential = credentialing(features1T1,features1T2,ppm = 20,rtwin = 2,rtcom = 5, ratio1 = 1/1, ratio2 = 1/2, ratio_tol = 0.1, ratios_tol=0.2)

# Check the results

  credential$CredentialedFeatureR2F # Final credentialed features with 2nd ratio filter
  credential$CredentialedFeatureR2 # Final credentialed features without 2nd ratio filter
