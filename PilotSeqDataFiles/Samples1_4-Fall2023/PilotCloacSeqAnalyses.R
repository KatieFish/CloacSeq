

#######Pilot data analyses - 10032024 KJF CH ZH

## required packages to load
require(reshape2)
require(iNEXT)
require(vegan)
require(ggplot2)

### Download data from GitHub
### Read in data from your folder. 
#this is my path, it won't be yours!

feature.table<-read.delim("~/qiim2_datafiles/pilotcloacseqfall23/CloacSeq/PilotSeqDataFiles/Samples1_4-Fall2023/pilot_seq_feature_table.txt")
taxonomy_assigmments<-read.delim("~/qiim2_datafiles/pilotcloacseqfall23/CloacSeq/PilotSeqDataFiles/Samples1_4-Fall2023/UNITE_euk_feature_taxonomy.tsv")

#Let's start by looking at total reads in each sample. 
sample.1.reads<-sum(feature.table$sample.1)
#sample 1 60K reads
sample.2.reads<-sum(feature.table$sample.2)
#sample 2 10K reads
sample.3.reads<-sum(feature.table$sample.3)
#sample 3 55K reads
sample.4.reads<-sum(feature.table$sample.4)
#sample 4 9K reads

#consistently about 5X less reads on qtip than swab

#Now lets generate an OTU frequency table by merging our two tables.

OTU_freq<- merge(taxonomy_assigmments, feature.table, by="Feature.ID", all=TRUE)

#let's get rid of concensus values for now. 
OTU_freq<-OTU_freq[c(1, 2, 4:7)]


#first, let's split off a table to make a figure summarizing kingdom level taxonomic assignments.
#import expanded taxonomy file.
taxon_expanded<-read.delim("~/qiim2_datafiles/pilotcloacseqfall23/CloacSeq/PilotSeqDataFiles/Samples1_4-Fall2023/Pilot_taxonomy_expanded.txt")
taxaplot_df<-merge(OTU_freq[c(1, 3:6)], taxon_expanded, by="Feature.ID", all=TRUE)

phylaplot<-melt(taxaplot_df[c(2:6)])

ggplot(phylaplot, aes(x=variable, y=value, fill=Kingdom))+
  ylab("#reads")+
  xlab("sample")+
  geom_col()

##This confirms what we saw before. From here on out we 
#want to look at only the fungal data, so let's generate a list of 
#Fungal Feature IDs to keep in our data. 
Fungal.IDs<- taxonomy_assigmments$Feature.ID[which(grepl("Fungi", taxonomy_assigmments$Taxon))]
##Keep only these in our data tables. 
OTU_freq<-OTU_freq[which(OTU_freq$Feature.ID %in% Fungal.IDs),]


#Ok - let's rarefy!

abundance_data<-list(OTU_freq$sample.1, OTU_freq$sample.2,
                     OTU_freq$sample.3, OTU_freq$sample.4)

names(abundance_data)<-c("Sample1", "Sample2", "Sample3", "Sample4")

inxt_results<-iNEXT(abundance_data, datatype = "abundance", q=0)

ggiNEXT(inxt_results, type=1, color.var="Assemblage")+
  xlab("Number of reads")+
  ylab("Number of OTUs")+
  theme_bw()


### What can we use these data to determine?


#Some diversity stats

stats<-inxt_results$AsyEst

ggplot(stats)+
  geom_errorbar(aes(x=Assemblage, ymin=(Estimator+s.e.), ymax=(Estimator-s.e.)), width=0.2, size=.5, color="darkgrey")+
  geom_point(aes(x=Assemblage, y=Estimator, col=Assemblage))+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.97))+
  facet_wrap(~Diversity, scales = "free")

