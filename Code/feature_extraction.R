#######################################

Pos = read.fasta('Pos.fasta', seqtype="AA", as.string = TRUE)
Neg = read.fasta('Neg.fasta', seqtype="AA", as.string = TRUE)

Pos <- Pos[(sapply(Pos, protcheck))]
Neg <- Neg[(sapply(Neg, protcheck))]
length(Pos)
length(Neg)
DD = c(Pos, Neg)

##############################################
################ Feature extraction ##########
##############################################

AAC = t(sapply(DD, extractAAC))
DPC = t(sapply(DD, extractDC))
TPC = t(sapply(DD, extractTC))
CTDC = t(sapply(DD, extractCTDC))
CTDD = t(sapply(DD, extractCTDD))
CTDT = t(sapply(DD, extractCTDT))
CTD = cbind(CTDC, CTDD, CTDT)

lamda = 1
PAAC1 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC1[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 2
PAAC2 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC2[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 3
PAAC3 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC3[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 4
PAAC4 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC4[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 5
PAAC5 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC5[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 6
PAAC6 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC6[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 7
PAAC7 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC7[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

lamda = 8
PAAC8 <- matrix(nrow = length(DD), ncol = 20 + lamda)
for(i in 1 : length(DD)) { 
PAAC8[i, ] = extractPAAC(DD[[i]][1], lambda = lamda,  props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

Class = c(rep("Pos", times= length(Pos)), rep("Neg", times= length(Neg)))
