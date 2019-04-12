# ontoNW

### What is ontoNW

**ontoNW** builds gene-ontology relationship static and dynamic network.

# Installation

The package can be installed and loaded with

```r
install_github("nahtan/ontoNW")
library(ontoNW)
```
# Using ontoNW

The main function in the **ontoNW** package is `ontoNW(onto, cutoff=0.05, pdf=T, html=T, dynamic=T)`.

See example below.

# Parameters

**onto** data.frame. Includes an 'id' column with term id (e.g. G0:0008150), a 'type' colum (e.g. BP, KEGG), a 'name' column (e.g. biological_process), a 'genes' column with string of genes of the query belonging to the ontology term (e.g. IL1RL1, IL5RA, TP53), and a 'fdr' column with FDR values.

**cutoff** numeric. A fdr value as cutoff threshold to plot edges/nodes below cutoff.

**dynamic** boolean. Produce a dynamic network.

**pdf** boolean. Saves the static network as pdf.

**html** boolean. Saves the dynamic network as html. Applies only if dynamic=T.

**palette** character. Rcolobrewer palette name. Default is 'Pastel2'.


# Example
```r
onto=data.frame(id  = c("GO:0002376", "GO:0032634")
                type= c("BP","BP")
                name= c("immune system process","interleukin-5 production")
                gene= c(c(IL5RA,SMPD3,CLC,IL1RL1,CYSLTR2,ALOX15,PIK3R6), c(IL5RA,IL1RL1))
                fdr = c(0.005, 0.00005)
                )
ontoNW(onto, cutoff=0.05, dynamic=T, pdf=F, html=T, palette="Accent")       
```

![ontoNWstatic](https://github.com/nahtan/ontoNW/blob/master/StaticNetwork_reingold_tilford.png)
