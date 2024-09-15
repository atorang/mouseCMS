# mouseCMS
CMS classifier for mouse models. 

The complete code used to generate figures and perform analyses in the article: **"Enterocyte-like differentiation defines the metabolic gene signatures of CMS3 colorectal cancers and provides a therapeutic vulnerability."** , along with the processed data required to run the scripts, can be accessed on GitHub.com/atorang/mouseCMS_Scripts, and synapse.org (id = syn53040345); respectively. 
[GitHub.com/atorang/mouseCMS_Scripts](https://github.com/atorang/mouseCMS_Scripts)


An example data set is provided to explain how to work with ```mouseCMS``` R package and ```mouseCMSclassifier``` function. 

To install the mouseCMS package use the following code:
```
devtools::install_github("atorang/mouseCMS")
```
The example data ```exprs.test``` is a data.frame object with log2-transformed expression levels. Rows are genes and columns are samples. Ensembl IDs are used.
```
library("mouseCMS")
res<-mouseCMSclassifier(exprs.test, perform.log2=FALSE)
```
