##Barplot for alignment statistics##
count <- read.table('./testforbarplot.txt',sep = "\t",header = F,col.names = c('Read_Type','Read_Number'))
count$Read_Type <- factor(count$Read_Type,levels = count$Read_Type)
ggplot(count)+geom_bar(aes(x=Read_Type,y=Read_Number),stat = "identity",fill=c("red","red","red","red","blue","blue",'blue'))+theme_bw() +theme(axis.text.x =element_text(angle = 45,hjust = 1,vjust = 1))
