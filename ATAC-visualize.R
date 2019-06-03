#ATAC-seq pipeline
#1.venndigram
library(VennDiagram)

venn1 <- read.xlsx("")
venn2 <- read.xlsx("")
venn3 <- read.xlsx("")

venn.diagram(list(CBSsh2vsctrl =venn1$gene_id,CBSsh23vsctrl=venn2$geneid),
             filename = "1_1.tif",
             col = "transparent",
             fill = c("cornflowerblue","darkorchid1"),
             label.col = "black",
             cat.col = c("blue","red"),
             cat.dist = c(0.03,0.03),
             cat.pos = c(-8,8),
             cat.cex = 1.5,
             cex = 2,
             main = "CBSsh DE genes cross",
             main.cex = 2)




















































