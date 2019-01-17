#!/usr/local/bin/Rscript

library(ggplot2)
library(readr)
#setwd('~/Downloads/AHDB/ampseq')
#---
# Read in data
#---
library(optparse)
opt_list = list(
    make_option("--OTU_table", type="character", help="OTUs table"),
    make_option("--prefix", type="character", help="Output prefix")
  )
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$OTU_table
prefix = opt$prefix

df1 <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)

# Create a varible that represents the split first column
x <- strsplit(df1$'#OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE)
# Use this to populate new columns for each element in this variable
df1$OTUs <- sapply( x, "[", 1 )
df1$Genus <- sapply( x, "[", 2 )
df1$Species <- sapply( x, "[", 3 )

# Edit species names:
df1$Species <- gsub("_SH\\d+","",df1$Species,ignore.case=F)
df1$Species <- gsub("_"," ",df1$Species,ignore.case=F)
df1$Species <- gsub("Alternaria oregonensis","Alternaria infectoria",df1$Species,ignore.case=F)
df1$Species <- gsub("Botrytis caroliniana","Botrytis cinerea",df1$Species,ignore.case=F)
df1$Species <- gsub("Ilyonectria mors-panacis","Cylindrocarpon destructans",df1$Species,ignore.case=F)
df1$Species <- gsub("Gibberella zeae","Fusarium graminearum",df1$Species,ignore.case=F)
df1$Species <- gsub("Gibberella tricincta","Fusarium redolens",df1$Species,ignore.case=F)
df1$Species <- gsub("Fusarium neocosmosporiellum","Fusarium solani",df1$Species,ignore.case=F)
df1$Species <- gsub("Monographella nivalis","Microdochium nivale",df1$Species,ignore.case=F)
df1$Species <- gsub("Didymella arachidicola","Phoma arachidicola",df1$Species,ignore.case=F)
df1$Species <- gsub("Thanatephorus cucumeris","Rhizoctonia solani",df1$Species,ignore.case=F)
df1$Species <- gsub("Sclerotium cepivorum","Stromatinia rapulum",df1$Species,ignore.case=F)
# df1$Species <- gsub("","",df1$Species,ignore.case=F)




df2 <- df1
# summary(df2)
# colnames(df2)[1] <- "OTUs"
# summary(df2)
rownames(df2) <- df2$OTUs
df2$'#OTU ID' <- NULL
df2$OTUs <- NULL
df2$Genus <- NULL
df2$Species <- NULL

df2 <- data.frame(t(df2))
df3 <- df2
df3$group <- gsub("_.$","",rownames(df3))

#install.packages('ggfortify')
library(ggfortify)

#---
# Function to summarise data to show SE and means by experimental treatment
#---

# --- Summary SE function ---
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
# --- /Summary SE function ends ---

#---
# Create Facetplots summarising data by experiment and pooled / technical reps
#---
#---
# Create plots combining data of OTUs with duplicate taxa designations
#---

library(reshape2)
df4 <- melt(df1, id.vars=c('#OTU ID',"OTUs", 'Genus', 'Species'),value.name = "Counts", variable.name='Run')
df4$'#OTU ID' <- NULL
df4$group <- gsub("_.$","",df4$Run)
df4$experiment <- gsub("_.$|_rep.$|_rep._.*$|_pooled_reps|_pool","",df4$Run)
df4$experiment <- gsub("^.*dilution","dilution",df4$experiment)
df4$experiment <- gsub("^.*equimolar","equimolar",df4$experiment)
df4$experiment <- gsub("equimolar.*$","equimolar",df4$experiment)
df4$experiment <- gsub("_"," ",df4$experiment)
df4$reps <- ifelse(grepl("pool",df4$Run), "pooled","tech. reps.")
df4$reps <- as.factor(df4$reps)
df4$reps <- factor(df4$reps,levels(df4$reps)[c(2,1)])

df6 <- aggregate(df4$Counts, by=list(df4$Run, df4$Species), FUN=sum)
colnames(df6) <- c('Run', 'Species', 'Counts')
df6$experiment <- gsub("_.$|_rep.$|_rep._.*$|_pooled_reps|_pool","",df6$Run)
df6$experiment <- gsub("^.*dilution","dilution",df6$experiment)
df6$experiment <- gsub("^.*equimolar","equimolar",df6$experiment)
df6$experiment <- gsub("equimolar.*$","equimolar",df6$experiment)
df6$experiment <- gsub("_"," ",df6$experiment)
# print(df6$experiment)
df6$reps <- ifelse(grepl("pool",df6$Run), "pooled","tech. reps.")
df6$reps <- as.factor(df6$reps)
df6$reps <- factor(df6$reps,levels(df6$reps)[c(2,1)])

df7 <- subset(df6, reps != "pooled")

# summary(df6)

# Normalise the data to 1000 reads


df8 <- aggregate(df7$Counts, by=list(df7$Run, df7$experiment, df7$reps), FUN=sum)

colnames(df8) <- c("Run", "experiment", "reps", "Total")

df9 <- merge(df7,df8,by=c("Run", "experiment", "reps"))

df9$norm <- ((df9$Counts / df9$Total) * 1000)
# summary(df9)
df10 <- summarySE(df9[ which(df9$Total > 1000),], measurevar="norm", groupvars=c(c("experiment", "Species")))
# summary(df10)
df10$id = numeric(nrow(df10))
for (i in 1:nrow(df10)){
   dat_temp <- df10[1:i,]
   df10[i,]$id <- nrow(dat_temp[dat_temp$Species == df10[i,]$Species,])
 }
dfx <- dcast(df10, id~Species, value = 'value', value.var = 'norm')
dfx$id <- NULL
dfx <- data.frame(apply(dfx, 2, sort, decreasing=T))
y <- as.logical(dfx[1,] > 10)
dfx <- data.frame(t(dfx))
dfx$keep <- y
dfx <- dfx['keep']
dfx$Species <- rownames(dfx)
dfx$Species <- gsub("\\."," ",dfx$Species,ignore.case=F)
dfx$Species <- gsub("Verticillium albo atrum","Verticillium albo-atrum",dfx$Species,ignore.case=F)

# dfx$Species
# df10$Species
dfy <- merge(dfx, df10, by="Species", all = TRUE)
df11 <- dfy[ which(dfy$keep == TRUE),]


# Create plot
facet_species<-ggplot(data=subset(df11), aes(x=Species, y=norm))
facet_species <- facet_species + geom_bar(stat="identity")
facet_species <- facet_species + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_species <- facet_species + ylab('Read count (per 1000 mapped reads)') + xlab('')
facet_species <- facet_species + geom_errorbar(aes(ymin=norm-se, ymax=norm+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_species <- facet_species + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
facet_species <- facet_species + geom_text(aes(label=round(norm)),  position = position_stack(vjust = 0.5))
# facet_species <- facet_species + facet_grid(experiment ~ reps)
facet_species <- facet_species + facet_grid(experiment ~ .)

filename <- paste(prefix, 'facet_species.pdf', sep='_')
ggsave(filename, plot = facet_species, width =20, height = 20, units = "cm", limitsize = FALSE)
#
#
# #---
# #
# #---
# # Write a table summarising mean reads by treatment
#
# summary_df <- aggregate(df6$Counts, by=list(df6$Species, df6$experiment, df6$reps), FUN=mean)
# colnames(summary_df) <- c('Species', 'experiment', 'reps', 'Counts')
# summary_df$Counts <- round(summary_df$Counts)
# # print(summary_df)
# filename <- paste(prefix, 'summarised_counts.tsv', sep='_')
# write.table(summary_df,file=filename,sep="\t", row.names=FALSE)
