##Barplot for alignment statistics##
#load library#
library(ggplot2)
library(optparse)


#set the input and output#
option_list <- list(
  make_option(c("-i","--input"),type="character", default="input.txt",
              help ="Input table is [default %default]"),
  make_option(c("-o","--output"),type="character", default="output.pdf",
              help ="Output table is [default %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))
print(paste("The input file is ", opts$input,  sep = ""))
print(paste("The output file is ", opts$output, sep = ""))

##load read count from bamqc##
count <- read.table(opts$input,sep = "\t",header = F,col.names = c('Read_Type','Read_Number'))
count$Read_Type <- factor(count$Read_Type,levels = count$Read_Type)

##output pdf##
pdf(opts$output, width = 9, height = 7)
ggplot(count)+geom_bar(aes(x=Read_Type,y=Read_Number),stat = "identity",fill=c("red","red","red","red","blue","blue",'blue'))+theme_bw() +theme(axis.text.x =element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()
