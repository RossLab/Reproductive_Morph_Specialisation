# Gene ontology (GO) enrichment analysis for differentially expressed autosomal genes #
With the list of differentially expressed autosomal genes, I wanted to see whether they were enriched for a specific gene ontology (GO) term, which would provide some hints as to the functional basis of their divergence. I do this using the `topGO` package. 
\
\
The necessary packages are
```
library(topGO) # v2.54.0
library(GO.db) # v3.18.0
library(dplyr) # v1.1.4
```

First I read in an annotation file that has scanned all the _B. coprophila_ proteins against the InterProScan database, which has GO terms associated. 
```
Bcop_genome_annotation <- read.csv("inputs/GO_Bcop_genome_annotation.csv", header = FALSE)
# Keep only the columns with the gene name and GO terms associated
Bcop_gene_GO <- Bcop_genome_annotation[,c(1,14)]
# Keep only the rows with actual GO terms
Bcop_gene_GO <- Bcop_gene_GO[grep("GO", Bcop_gene_GO$V14), ]
# Get rid of the .t transcript numbers
Bcop_gene_GO$V1 <- gsub("\\.t\\d+", "", Bcop_gene_GO$V1)
# Reassign column names 
colnames(Bcop_gene_GO) <- c("GeneID", "GO_Term")
head(Bcop_gene_GO)

# From the file there are some repeated rows, so use distinct form the dplyr package to keep only unique whole rows
Bcop_gene_GO <- distinct(Bcop_gene_GO)
# Reformat to group by GeneID so I have one column with unique GeneID and another column with all the GO terms associated 
Bcop_gene_GO <- Bcop_gene_GO %>%
  group_by(GeneID) %>%
  summarize(GO_Term = paste(GO_Term, collapse = ", "))
Bcop_gene_GO$GO_Term <- gsub("\\|", ", ", Bcop_gene_GO$GO_Term)
# write.csv(Bcop_gene_GO, file = "inputs/Bcop_gene_to_GO.csv", row.names = FALSE)
# write.table(Bcop_gene_GO, file = "inputs/Bcop_gene_to_GO.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
```
Now we have the first file ready, which is a gene-to-GO mapping. I also read in a file that has all the annotated genes in the genome and the chromosomes they belong on.
```
genes_and_chromosomes <- read.table("inputs/Bcop_v3_genes_and_chromosomes.tsv")
autosomal_genes <- genes_and_chromosomes[genes_and_chromosomes$V1 == "II" | genes_and_chromosomes$V1 == "III" | genes_and_chromosomes$V1 == "IV", ]
```
\
\
Now we can start GO-Term Enrichment Analysis, first on the somatic non-reproductive tissue. The principle of this analysis is that you have a "gene universe", which in this case is all the autosomal genes with GO terms. You also have your "genes of interest", which in this case are the autosomal genes that are differentially expressed between gynogenic and androgenic females in somatic non-reproductive tissue. You want to see if there are GO terms that appear in your "genes of interest" more commonly than expected by chance. 
\
\
With topGO, you first build a topGO object with all this information before you can apply statistical analysis. I'm filtering the "gene universe" here to only have autosomal genes. 
```
### Start loading in the files necessary to build a topGO object ###
# Loads in the mappings of geneID to GO term
geneID2GO <- readMappings("inputs/Bcop_gene_to_GO.tsv", sep = "\t", IDsep = ",")
keep <- names(geneID2GO) %in% autosomal_genes$V2
geneID2GO <- geneID2GO[keep]
length(names(geneID2GO))
# The gene universe in this case is given by the list names
geneNames <- names(geneID2GO)
# Load in the list of interesting genes
InterestingGenes <- read.csv("inputs/nonrepro_androvsgyno_sig.csv", header = TRUE)
nrow(InterestingGenes) #393
# Combine information so I have a list of genes, as well as information about whether they're interesting or not
geneList <- factor(as.integer(geneNames %in% InterestingGenes))
names(geneList) <- geneNames
```
\
\
We can then start the GO term enrichment analysis. There are three main categories of GO terms: Biological Processes (BP), Molecular Function (MF), and Cellular Component (CC). I first look at Biological Processes. 
```
### Build our GO object ###
# annFUN.gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping
GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

# Now that we have our topGO object, we can look at it 
graph(GOdataBP)
termStat(GOdataBP)

# GO analysis 
# The main function is getSigGroups() which takes two parameters.The first parameter is of class topGOdata and the second parameter is of class groupStats
# So you have to make the second parameter, specifying what kind of stats you want to be applied to all the GO 
# Here we're using Fisher's exact test 
# Specify Fisher's exact test
test.stat.Fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

# Use it with the GO object we've created
resultFisher <- getSigGroups(GOdataBP, test.stat.Fisher)
## weight01 is the default algorithm, since I've not specified a specific algorithm 
resultFisher

# Keep in mind that these p values are not adjusted for multiple testing 
allResBP <- GenTable(GOdataBP, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 50, numChar=1000)
allResBP
write.csv(allResBP, file = "outputs/nonrepro_BP_Fisher_top50.csv")
png(filename = "outputs/nonrepro_BP_GOTerms.png", width = 1000, height = 1000, res = 300)
graphBP<- showSigOfNodes(GOdataBP, score(resultFisher), firstSigNodes = 5, useInfo = "all")
dev.off()
```
You can save your significally enriched GO terms as a table or as a graph.
\
\
I now do the same with the other categories (MF and CC)
```
## MF (molecular function) ##
GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
graph(GOdataMF)
termStat(GOdataMF)
resultFisherMF <- getSigGroups(GOdataMF, test.stat.Fisher)
resultFisherMF
# Keep in mind that these p values are not adjusted for multiple testing 
allResMF <- GenTable(GOdataMF, classicFisher = resultFisherMF, orderBy = "classicFisher", topNodes = 50, numChar=1000)
allResMF
write.csv(allResMF , file = "outputs/nonrepro_MF_Fisher_top50.csv")
graphMF<- showSigOfNodes(GOdataMF, score(resultFisherMF), firstSigNodes = 7, useInfo = "all")

## CC (cellular component) ##
GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
graph(GOdataCC)
termStat(GOdataCC)
resultFisherCC <- getSigGroups(GOdataCC, test.stat.Fisher)
resultFisherCC
# Keep in mind that these p values are not adjusted for multiple testing 
allResCC <- GenTable(GOdataCC, classicFisher = resultFisherCC, orderBy = "classicFisher", topNodes = 50, numChar=1000)
allResCC
write.csv(allResCC , file = "outputs/nonrepro_CC_Fisher_top50.csv")
graphCC<- showSigOfNodes(GOdataCC, score(resultFisherCC), firstSigNodes = 5, useInfo = "all")
```
\
\
And I repeat the same process for genes that are significantly differentially expressed between gynogenic and androgenic females in the germline.
```
#### GERMLINE TISSUE ####
#### GO analysis ####
### Start loading in the files necessary to build a topGO object ###
# Loads in the mappings of geneID to GO term
geneID2GO <- readMappings("inputs/Bcop_gene_to_GO.tsv", sep = "\t", IDsep = ",")
str(head(geneID2GO))
length(names(geneID2GO))
keep <- names(geneID2GO) %in% autosomal_genes$V2
geneID2GO <- geneID2GO[keep]
length(names(geneID2GO))
# The gene universe in this case is given by the list names
geneNames <- names(geneID2GO)
head(geneNames)
# Load in the list of interesting genes, which is the consensus list from the somatic non-reproductive data.
InterestingGenes <- read.csv("inputs/germline_androvsgyno_sig.csv", header = TRUE)
nrow(InterestingGenes) #11
InterestingGenes <- unique(InterestingGenes$X) # make sure there are no duplicates
length(InterestingGenes) #11
# Combine information so I have a list of genes, as well as information about whether they're interesting or not
geneList <- factor(as.integer(geneNames %in% InterestingGenes))
names(geneList) <- geneNames
str(geneList)

### Now that we have everything, we can build our GO object ###
# annFUN.gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping
GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

# Now that we have our topGO object, we can look at it 
graph(GOdataBP)
termStat(GOdataBP)

# GO analysis 
# The main function is getSigGroups() which takes two parameters.The first parameter is of class topGOdata and the second parameter is of class groupStats
# So you have to make the second parameter, specifying what kind of stats you want to be applied to all the GO 
# Here we're using Fisher's exact test 
# Specify Fisher's exact test
test.stat.Fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

# Use it with the GO object we've created
resultFisher <- getSigGroups(GOdataBP, test.stat.Fisher)
## weight01 is the default algorithm, since I've not specified a specific algorithm 
resultFisher

# Keep in mind that these p values are not adjusted for multiple testing 
allResBP <- GenTable(GOdataBP, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 50, numChar=1000)
allResBP
write.csv(allResBP, file = "outputs/germline_BP_Fisher_top50.csv")
graphBP<- showSigOfNodes(GOdataBP, score(resultFisher), firstSigNodes = 7, useInfo = "all")

# This should give the same result as above because it's just a more user-friendly interface to run the same test with the same algorithm
weight01.fisher.BP <- runTest(GOdataBP, statistic = "fisher", algorithm = "weight01")
allResBP2 <- GenTable(GOdataBP, P_Value = weight01.fisher.BP, orderBy = P_Value, topNodes = 20)
allResBP2
# But it gives different p values, which is weird 
graphBP2<- showSigOfNodes(GOdataBP, score(weight01.fisher.BP), firstSigNodes = 7, useInfo = "all")

## Trying this with other gene ontology categories
## MF (molecular function) ##
GOdataMF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
graph(GOdataMF)
termStat(GOdataMF)
resultFisherMF <- getSigGroups(GOdataMF, test.stat.Fisher)
resultFisherMF
# Keep in mind that these p values are not adjusted for multiple testing 
allResMF <- GenTable(GOdataMF, classicFisher = resultFisherMF, orderBy = "classicFisher", topNodes = 50, numChar=1000)
allResMF
write.csv(allResMF , file = "outputs/germline_MF_Fisher_top50.csv")
graphMF<- showSigOfNodes(GOdataMF, score(resultFisherMF), firstSigNodes = 7, useInfo = "all")

## CC (cellular component) ##
GOdataCC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
graph(GOdataCC)
termStat(GOdataCC)
resultFisherCC <- getSigGroups(GOdataCC, test.stat.Fisher)
resultFisherCC
# Keep in mind that these p values are not adjusted for multiple testing 
allResCC <- GenTable(GOdataCC, classicFisher = resultFisherCC, orderBy = "classicFisher", topNodes = 50, numChar=1000)
allResCC
write.csv(allResCC , file = "outputs/germline_CC_Fisher_top50.csv")
graphCC<- showSigOfNodes(GOdataCC, score(resultFisherCC), firstSigNodes = 5, useInfo = "all")
```
