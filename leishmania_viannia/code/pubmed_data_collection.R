##### Collecting pubmed data (*code fragment*)
##### title: Collecting pubmed data
##### author: Aisling Murran, Caroline Glidden
##### version: 06/01/2021

######### READ_ME: 
# This code was used after the total_data dataframe was created in order to download
# to collect information about the number of times each species showed up in pubmed. 

# *WARNING* This code does not work in its current form because the total_data 
# dataframe does not exist in this file. This code is a fragment which was used to 
# to create the 'neotropical species pubmed citations.csv' file used in 
# 'clean and merge data_4.1.R'


library(RISmed) #downloading pubmed citations
total_data <- readRDS('../cleaned data/analysis data/total data.rds')

##Collect the names of all the species in our total_data set
names <- paste(total_data$MSW05_Binomial, "AND (parasite or pathogen)")
names <- unique(names)

#Count the number of times they appear in pubmed
pubmed.count <-  sapply(names, 
                        function(i) {
                          QueryCount(EUtilsSummary(i,
                                                   type = 'esearch',db = 'pubmed'))
                        }
)

pubmed.count <- as.data.frame(pubmed.count)
pubmed.count.backup <- pubmed.count

###split column to get species name again
library(dplyr)
library(tidyr)
library(stringr)

pubmed.count$names <- row.names(pubmed.count)

final <- pubmed.count %>% separate(names, c('MSW05_Binomial', 'Extra'), sep=' AND ')

final <- final[,c("pubmed.count","MSW05_Binomial")]

write.csv(final, "pubmed_count_disease_related.csv")






