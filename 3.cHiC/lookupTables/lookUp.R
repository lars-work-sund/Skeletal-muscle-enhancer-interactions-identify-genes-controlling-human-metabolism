library(data.table)
load("combinedResults.RData")

###############################################################################
#    Example 1, acessing individual data
###############################################################################
# Enhancers response to palmitate
enhancerResults$Palmitate

# Enhancers response to TNFa
enhancerResults$TNFa

# RNA response to palmitate
rnaResults$Palmitate

# TNA response to TNFa
rnaResults$TNFa

# HiC interactiome response to palmitate
mergedHiC$Palmitate

# HiC interactiome response to TNFa
mergedHiC$TNFa

# HiC interactiome response to Passages (in the control condition)
mergedHiC$Passage


###############################################################################
#    Example 2, going from promoters to enhancers, or vice versa
###############################################################################
# The interactions are the same across all conditions
# Going from promoters to enhancers
promoters2enhancers

# Going from enhancers to promoters
enhancers2promoters

# You can look up one or more genes or enhancers
promoters2enhancers["ENSG00000000419"]
promoters2enhancers[c("ENSG00000000419", "ENSG00000000971")]

enhancers2promoters["peak.1000"]
enhancers2promoters[c("peak.1000", "peak.1010")]

# If there are too many interaction partners they may not get printed, you can 
# force the printing of all like this (compare with the results on line 39)
promoters2enhancers[c("ENSG00000000419", "ENSG00000000971"), enhancerID]

# Same for hicIDs and promoterIDs
enhancers2promoters[c("peak.1000", "peak.1010"), promoterID]
enhancers2promoters[c("peak.1000", "peak.1010"), hicIDs]


###############################################################################
#    Example 3, subsetting results
###############################################################################
# You can filter the results you get based on any information in the cloumns
# All palmitate results
rnaResults$Palmitate
# Palmitate results with a FDR less than 5%
rnaResults$Palmitate[FDR < 0.05, ]

# Palmitate results with a FDR less than 5% and a logFC greater than 1
rnaResults$Palmitate[FDR < 0.05, ][logFC > 1, ]

# Palmitate results with a FDR less than 5% and 
# an absolute logFC greater than 1
rnaResults$Palmitate[FDR < 0.05, ][abs(logFC) > 1, ]

# Palmitate results where the gene is named PDK4
rnaResults$Palmitate[SYMBOL %in% c("PDK4"), ]

# Palmitate results where the gene is named PDK4 or MYLIP
rnaResults$Palmitate[SYMBOL %in% c("PDK4", "MYLIP"), ]

# You can extract the ENSEMBL IDs of your selections
rnaResults$Palmitate[SYMBOL %in% c("PDK4", "MYLIP"), ENSEMBL]

# and store them in a variable (will be used below)
myEnsIDs <- rnaResults$Palmitate[SYMBOL %in% c("PDK4", "MYLIP"), ENSEMBL]

# Or you can write them manually
myEnsIDs2 <- c("ENSG00000004799", "ENSG00000007944")

###############################################################################
#    Example 4, putting it together
###############################################################################
# Which enhacers interact with the genes we sleceted in the example above?
promoters2enhancers[myEnsIDs, ]
promoters2enhancers[myEnsIDs, enhancerID]
# In order to select the enhancers the IDs need to be in one list, and not
# a list of lists. This is achieved with the unlist function
unlist(promoters2enhancers[myEnsIDs, enhancerID])
myPeakIDs <- unlist(promoters2enhancers[myEnsIDs, enhancerID])

# Finally lets look at these enhancers
enhancerResults$Palmitate[myPeakIDs]

# Lets keep only the significant ones
enhancerResults$Palmitate[myPeakIDs][FDR < 0.05, ]

# Only one enhancer, is it involved in other interactions?
# Save the enhancer name
myEnhancer <- enhancerResults$Palmitate[myPeakIDs][FDR < 0.05, rn]

# Which promoters is it interacting with?
enhancers2promoters[myEnhancer, ]

# Lets look at both the genes and the interactions (remember to unlist)
myInt <- unlist(enhancers2promoters[myEnhancer, hicIDs])
myEnsID3 <- unlist(enhancers2promoters[myEnhancer, promoterID])

# Lets have alook at the genes it interacts with:
rnaResults$Palmitate[myEnsID3]
# One is apparently not expressed in this tissue

# Lets have a look at the interactions
mergedHiC$Palmitate[myInt]
