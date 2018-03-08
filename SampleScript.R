# Credential3.1
A bioinformatics tool for biologically-relavent feature extraction features and untargeted metabolomics method optimization.
 
# Installation
```
devtools::install_github("pattilabwu/Credential3.1",ref="test")
library(Credential3.1)
```

# Input data format
The input data is a peaktable that contains at least 4 columns of data: "cc"-index No, "mz"- m/z value, "rt"- retention time and "i"-intensity
This peaktable can be generated from XCMS or many other peak-picking algorithms. 
The name of the columns in the peaktable has to be "mz", "cc", "rt" and "i".

```

> features
      cc       mz       rt            i
 1: 4850 550.3067  851.529     5419.261
 2: 5013 558.3340  850.141     5854.929
 3: 4804 548.3008  850.391     5209.042
 4: 5029 560.3396  850.642     5998.703
 5: 4827 549.4886 1536.830 10604955.633
 6: 4852 550.4924 1536.830  4241038.103
 7: 4874 551.4955 1536.580  1206786.640
 8: 4896 552.4979 1536.590   139663.706
 9: 5282 580.5907 1536.580    28051.912
10: 5303 581.5949 1536.580   202953.667
11: 5315 582.5990 1536.580  1286416.952
12: 5324 583.6029 1536.330  5415255.371
```

# Usage

```
credential = credentialing(peaktable1,peaktable2,ppm = 20,rtwin = 2,rtcom = 1, ratio1 = 1/1, ratio2 = 1/2, ratio_tol = 0.1, charges = 1:4,export = F)
```

For detailed illustration of the parameters, see:  
```
help(credentialing) 
```

# data output
A list of tables will be generated and compiled in the 'credential' object shown above. It includes:

CredentialedFeatureR2F -- final round credentialed features from both peaktables

CredentialedFeatureR2 -- mz&rt matched features from both peaktables

CredentialedFeature1N2 -- features that fail to pass 2nd filtering from credentialed features #1

CredentialedFeature2N2 -- features that fail to pass 2nd filtering from credentialed features #2

CredentialedFeature1R1 -- first round credentialed features from peaktable #1

CredentialedFeature2R1 -- first round credentialed features from peaktable #2

knots1&knots2 -- possible isotope pairs from peaktable #1&#2

credentialedKnots1&credentialedKnots2 -- credentialed groups from isotope pairs table #1&#2

CredentialedFeatureGroups -- matched index between group #1&#2("quipu" is the index for credentialed groups)

# Understanding the output data

#65 credentialed feature groups:

```
credential7$CredentialedFeatureR2F[quipu_1==65,]
   cc_1     mz_1     rt_1      int_1 knot_1 tail_1 quipu_1  ratio_1 cc_2     mz_2    rt_2      int_2 knot_2 tail_2 quipu_2  ratio_2 combined_ratio
1: 8843 833.6024 1594.130 123599.468   1090      0      65 8.961008 8854 833.6031 1593.79  77270.810    596      0     679 3.532957       2.536404
2: 8856 834.6063 1594.360  66766.104   1090      0      65 8.961008 8869 834.6067 1594.04  38324.387    596      0     679 3.532957       2.536404
3: 9336 872.7280 1593.995   6910.652   2180      0      65 8.961008 9333 872.7280 1593.79   7342.721   1578      0     679 3.532957       2.536404
4: 9346 873.7326 1594.360  29194.090   2180      0      65 8.961008 9344 873.7334 1593.54  40671.305   1578      0     679 3.532957       2.536404
5: 9354 874.7378 1593.880  96785.679   2180      0      65 8.961008 9353 874.7382 1593.54 131249.741   1578      0     679 3.532957       2.536404
6: 9364 875.7416 1594.110 192621.899   2180      0      65 8.961008 9365 875.7426 1593.54 230820.709   1578      0     679 3.532957       2.536404
7: 9372 876.7447 1593.860  13793.032   2180      1      65 8.961008 9377 876.7464 1593.21  21871.426   1578      1     679 3.532957       2.536404

```
cc_1/2: peak index number

mz_1/2: measured mz value

rt_1/2: retention time in seconds

int_1/2: intensity

knot_1/2: index for isotope pairs

tail_1/2: U13C peak or not

quipu_1/2: index for 1st round credentialed groups within each ratio condition

ratio_1/2: U12C/U13C ratio

Peaks from group1 and group2 are aligned together based on the mz and rt time of the monoisotopic mass within each credentialed group. 2nd credentialing is achieved by filtering the combined_ratio of the matches.
