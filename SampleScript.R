
# Step0: Installation of pre-requisite R packages and setting working directory
source("https://bioconductor.org/biocLite.R")
biocLite("xcms")
install.packages("data.table")
install.packages("devtools")
devtools::install_github("pattilabwu/Credential3.1")

library(xcms)
library(Credential3.1)
library(data.table)

# Set up working directory
#format for Mac OS:
setwd("/Volumes/Users/Lingjue/CredentialingDemo")
#Format for Windows:
setwd("C:\Users\Mike\Desktop\CredentialingDemo")


# Step 1: Data-preprocessing with XCMS

# suppose your mzXML files are stored in "1T1" and "1T2" folder under the working directory.
# For "1T1" condition
# Chromatographic peak detection
xs_1 = xcmsSet("./Ratio1", method="centWave", ppm=20, peakwidth=c(5,30), snthresh=5, prefilter=c(3,1000))
# Grouping peaks across the samples.
xs_1 = group(xs_1, bw=5, mzwid=.015, minfrac=0.5)
# Retention time re-alignment(correction)
xs_1 = retcor(xs_1, method="obiwarp",profStep=1)
# Grouping again after rt alignments
xs_1 = group(xs_1, bw=5, mzwid=.015, minfrac=0.5)
# Integrate areas of missing peaks in specific samples
xs_1 = fillPeaks(xs_1)

# For "1T2" condition
xs_2 = xcmsSet("./Ratio2", method="centWave", ppm=20, peakwidth=c(5,30), snthresh=5, prefilter=c(3,1000))
xs_2 = group(xs_2, bw=5, mzwid=.015, minfrac=0.5)
xs_2 = retcor(xs_2, method="obiwarp",profStep=1)
xs_2 = group(xs_2, bw=5, mzwid=.015, minfrac=0.5)
xs_2 = fillPeaks(xs_2)



# Step 2: Feature extraction from XCMS

featureGroup1 = xsGroupExtract(xs_1,intchoice = "into", sampling = 1)
featureGroup2 = xsGroupExtract(xs_2,intchoice = "into", sampling = 1)

# Check the results

featureGroup1$credTable    # the format for credentialing input
featureGroup1$featuretable # detailed table for features


# Step 3: Credentialing extracted features and check results

credential = credentialing(featureGroup1$credTable,featureGroup2$credTable,ppm = 20,rtwin = 2,rtcom = 1, ratio1 = 1/1, ratio2 = 1/2, ratio_tol = 0.1, charges = 1:4)

# Check the results

credential$CredentialedFeatureR2F # Final credentialed features with 2nd ratio filter
credential$CredentialedFeatureR2 # Final credentialed features without 2nd ratio filter
