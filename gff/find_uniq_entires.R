library(tidyverse) 
#Erbay Yigit 2021.08.12
# USAGE make df from GFF file and do some opereations.
# gff3 files columnas are organized as below:
#**********
# gff structure
#seq_id    source  type   start end score   strand  phase   attributes
#**********
df0<-readr::read_delim('BestRefSeq_mRNA.gff',
  col_names=c("seq_id","source","type","start","end","score","strand","phase", "attributes"),
                      comment = '#', delim = '\t', na='.')

df<-as_tibble(df0)

#these are form dplyr. It will find uniqs based on two columns
x<- df %>% distinct(df$"start", df$"end", .keep_all = TRUE )
