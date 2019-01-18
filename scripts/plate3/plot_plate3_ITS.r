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



df1 <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
# df1 <- read.delim("~/Downloads/AHDB_new/plate3/quantified/mix-A/daffodil/mix-A_ITS_zOTUs_zOTUs_table.txt")
colnames(df1)[1] <- "OTU ID"
# Create a varible that represents the split first column
df1$'OTU ID' <- as.character(df1$'OTU ID')
x <- strsplit(df1$'OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE)
df1$Genus <- sapply( x, "[", 2)
df1$Species <- sapply( x, "[", 3)
df1$Species <- gsub("_SH\\d+","",df1$Species,ignore.case=F)
df1$Species <- gsub("_"," ",df1$Species,ignore.case=F)
df1$Species = ifelse(df1$Species == 'unknown', df1$Genus, df1$Species)
df1$Species <- gsub("Epicoccum pimprinum","Phoma pimprina",df1$Species,ignore.case=F)
df1$Genus = ifelse(df1$Species == 'Phoma pimprina', "Phoma", df1$Genus)
df1$Genus <- gsub("Gibberella","Fusarium",df1$Genus,ignore.case=F)
df1$Genus <- gsub("/"," ",df1$Genus,ignore.case=F)
df1$Genus <- gsub("\\-","",df1$Genus,ignore.case=F)
df1$Species <- gsub("/"," ",df1$Species,ignore.case=F)
df1$Species <- gsub("\\-","",df1$Species,ignore.case=F)
# df1$Species[df1$Species == 'unknown'] <- df1$Genus[df1$Species == 'unkown']

# summary(df1)
df1$Genus <- as.factor(df1$Genus)
# df1$Species <- as.factor(df1$Species)

#--
# By Genus
#--

library(reshape2)
df2 <- melt(df1, id.vars=c('OTU ID', 'Genus', 'Species'),value.name = "Counts", variable.name='Run')
df2$Run <- as.character(df2$Run)
y <- strsplit(df2$Run, '_', fixed = FALSE, perl = FALSE, useBytes = FALSE)
df2$Mix <- sapply( y, "[", 2)
df2$Field <- sapply( y, "[", 3)
df2$Bed <- sapply( y, "[", 4)
df2$Sample <- sapply( y, "[", 5)
df2$Locus <- sapply( y, "[", 6)

df3 <- aggregate(df2$Counts, by=list(df2$Genus, df2$Mix, df2$Field, df2$Bed, df2$Sample, df2$Locus), FUN=sum)
colnames(df3) <- c("Genus", "Mix", "Field", "Bed", "Sample", "Locus", "Counts")

df4 <- aggregate(df2$Counts, by=list(df2$Mix, df2$Field, df2$Bed, df2$Sample, df2$Locus), FUN=sum)
colnames(df4) <- c("Mix", "Field", "Bed", "Sample", "Locus", "Total")

df5 <- merge(df3,df4,by=c("Mix", "Field", "Bed", "Sample", "Locus"))

df5$norm <- ((df5$Counts / df5$Total) * 1000)
df6 <- summarySE(df5[ which((df5$Total > 1000) & (df5$Bed != 'pool')),], measurevar="norm", groupvars=c(c("Mix", "Field", "Bed", "Locus", "Genus")))
# df7 <- df6[ which(df6$norm > 10),]
df6$id = numeric(nrow(df6))
for (i in 1:nrow(df6)){
   dat_temp <- df6[1:i,]
   df6[i,]$id <- nrow(dat_temp[dat_temp$Genus == df6[i,]$Genus,])
 }
dfx <- dcast(df6, id~Genus, value = 'value', value.var = 'norm')
dfx$id <- NULL
dfx <- data.frame(apply(dfx, 2, sort, decreasing=T))
y <- as.logical(dfx[1,] > 10)
dfx <- data.frame(t(dfx))
dfx$keep <- y
dfx <- dfx['keep']
dfx$Genus <- rownames(dfx)
dfx$Genus <- gsub("\\."," ",dfx$Genus,ignore.case=F)
dfx$Genus <- gsub(" incertae sedis ","(incertae sedis)",dfx$Genus,ignore.case=F)

# dfx$Genus
# df6$Genus
dfy <- merge(dfx, df6, by="Genus", all = TRUE)
df7 <- dfy[ which(dfy$keep == TRUE),]



facet_Genus<-ggplot(data=subset(df7), aes(x=Genus, y=norm))
# facet_Genus <- facet_Genus + scale_y_continuous(limits = c(0, 1000))
facet_Genus <- facet_Genus + geom_bar(stat="identity")
facet_Genus <- facet_Genus + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_Genus <- facet_Genus + ylab('Read count (per 1000 mapped reads)') + xlab('')
facet_Genus <- facet_Genus + geom_errorbar(aes(ymin=norm-se, ymax=norm+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_Genus <- facet_Genus + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
facet_Genus <- facet_Genus + facet_grid(Bed ~ .)
# facet_Genus <- facet_Genus + geom_text(aes(label=round(norm)), vjust=-2.5)
facet_Genus <- facet_Genus + geom_text(aes(label=round(norm)),  position = position_stack(vjust = 0.5))
facet_Genus

# prefix <- '/Users/armita/Downloads/AHDB_new/plate3/plots/TEF_FOC'
fig_height <- (length(unique(df2$Bed))*5)+5
filename <- paste(prefix, 'facet_genus.pdf', sep='_')
ggsave(filename, plot = facet_Genus, width =20, height = fig_height, units = "cm", limitsize = FALSE)

# --
# By species
# --

library(reshape2)
df2 <- melt(df1, id.vars=c('OTU ID', 'Genus', 'Species'),value.name = "Counts", variable.name='Run')
df2$Run <- as.character(df2$Run)
y <- strsplit(df2$Run, '_', fixed = FALSE, perl = FALSE, useBytes = FALSE)
df2$Mix <- sapply( y, "[", 2)
df2$Field <- sapply( y, "[", 3)
df2$Bed <- sapply( y, "[", 4)
df2$Sample <- sapply( y, "[", 5)
df2$Locus <- sapply( y, "[", 6)

df3 <- aggregate(df2$Counts, by=list(df2$Species, df2$Mix, df2$Field, df2$Bed, df2$Sample, df2$Locus), FUN=sum)
colnames(df3) <- c("Species", "Mix", "Field", "Bed", "Sample", "Locus", "Counts")

df4 <- aggregate(df2$Counts, by=list(df2$Mix, df2$Field, df2$Bed, df2$Sample, df2$Locus), FUN=sum)
colnames(df4) <- c("Mix", "Field", "Bed", "Sample", "Locus", "Total")

df5 <- merge(df3,df4,by=c("Mix", "Field", "Bed", "Sample", "Locus"))

df5$norm <- ((df5$Counts / df5$Total) * 1000)
df6 <- summarySE(df5[ which((df5$Total > 1000) & (df5$Bed != 'pool')),], measurevar="norm", groupvars=c(c("Mix", "Field", "Bed", "Locus", "Species")))
df6$id = numeric(nrow(df6))
for (i in 1:nrow(df6)){
   dat_temp <- df6[1:i,]
   df6[i,]$id <- nrow(dat_temp[dat_temp$Species == df6[i,]$Species,])
 }
dfx <- dcast(df6, id~Species, value = 'value', value.var = 'norm')
dfx$id <- NULL
dfx <- data.frame(apply(dfx, 2, sort, decreasing=T))
y <- as.logical(dfx[1,] > 10)
dfx <- data.frame(t(dfx))
dfx$keep <- y
dfx <- dfx['keep']
dfx$Species <- rownames(dfx)
dfx$Species <- gsub("\\."," ",dfx$Species,ignore.case=F)

# dfx$Species
# df6$Species
dfy <- merge(dfx, df6, by="Species", all = TRUE)
df7 <- dfy[ which(dfy$keep == TRUE),]


facet_species<-ggplot(data=subset(df7), aes(x=Species, y=norm))
# facet_species <- facet_species + scale_y_continuous(limits = c(0, 1000))
facet_species <- facet_species + geom_bar(stat="identity")
facet_species <- facet_species + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_species <- facet_species + ylab('Read count (per 1000 mapped reads)') + xlab('')
facet_species <- facet_species + geom_errorbar(aes(ymin=norm-se, ymax=norm+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_species <- facet_species + theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
facet_species <- facet_species + facet_grid(Bed ~ .)
# facet_species <- facet_species + geom_text(aes(label=round(norm)), vjust=-2.5)
facet_species <- facet_species + geom_text(aes(label=round(norm)),  position = position_stack(vjust = 0.5))
facet_species

# prefix <- '/Users/armita/Downloads/AHDB_new/plate3/plots/TEF_FOC'
fig_height <- (length(unique(df2$Bed))*10)+5
filename <- paste(prefix, 'facet_species.pdf', sep='_')
ggsave(filename, plot = facet_species, width =30, height = fig_height, units = "cm", limitsize = FALSE)
