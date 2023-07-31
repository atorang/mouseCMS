# mouseCMS
CMS classifier for mouse models. 

An example data set is provided to explain how to work with mouseCMS R package and mouseCMSclassifier function. 

To install the mouseCMS package use the following code:
```
devtools::install_github("atorang/mouseCMS")
```
The example data ```exprs.test``` is a data.frame object with log2-transformed expression levels. Rows are genes and columns are samples. Ensembl IDs are used.
```
library("mouseCMS")
res<-mouseCMSclassifier(exprs.test, perform.log2=FALSE)
```
