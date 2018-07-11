#!/usr/bin/Rscript

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

# ITS
# OTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/ITS/soil_pathogens/ITS_OTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/ITS/soil_pathogens/ITS_OTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# zOTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/ITS/soil_pathogens/ITS_zOTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# Normalised data
#df1 <- read_delim("~/Downloads/AHDB/ampseq/ITS/soil_pathogens/ITS_zOTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#---
# TEF
#---
# -
# Soil pathogens
# -
#   OTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/soil_pathogens/TEF_OTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_soil_pathogens_OTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/soil_pathogens/TEF_OTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_soil_pathogens_OTUs_norm'
#   zOTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/soil_pathogens/TEF_zOTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_soil_pathogens_zOTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/soil_pathogens/TEF_zOTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_soil_pathogens_zOTUs_norm'
# -
# Fusarium spp.
# -
#   OTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/Fusarium_spp/TEF_OTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_Fusarium_spp_OTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/Fusarium_spp/TEF_OTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_Fusarium_spp_OTUs_norm'
#   zOTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/Fusarium_spp/TEF_zOTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_Fusarium_spp_zOTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/TEF/Fusarium_spp/TEF_zOTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'TEF_Fusarium_spp_zOTUs_norm'
#---
# SIX13
#---
# -
# Soil pathogens
# -
#
# -
# Fusarium spp.
# -
#   OTUs
# raw data
#df1 <- read_delim("~/Downloads/AHDB/ampseq/SIX13/Fusarium_spp/SIX13_OTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#prefix <- 'SIX13_Fusarium_spp_OTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/SIX13/Fusarium_spp/SIX13_OTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'SIX13_Fusarium_spp_OTUs_norm'
#   zOTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/SIX13/Fusarium_spp/SIX13_zOTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'SIX13_Fusarium_spp_zOTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/SIX13/Fusarium_spp/SIX13_zOTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'SIX13_Fusarium_spp_zOTUs_norm'
#---
# OG12981
#---
# -
# Soil pathogens
# -
#
# -
# Fusarium spp.
# -
#   OTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG12981/Fusarium_spp/OG12981_OTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG12981_Fusarium_spp_OTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG12981/Fusarium_spp/OG12981_OTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG12981_Fusarium_spp_OTUs_norm'
#   zOTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG12981/Fusarium_spp/OG12981_zOTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG12981_Fusarium_spp_zOTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG12981/Fusarium_spp/OG12981_zOTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG12981_Fusarium_spp_zOTUs_norm
#
#---
# OG13890
#---
# -
# Soil pathogens
# -
#
# -
# Fusarium spp.
# -
#   OTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG13890/Fusarium_spp/OG13890_OTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG13890_Fusarium_spp_OTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG13890/Fusarium_spp/OG13890_OTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG13890_Fusarium_spp_OTUs_norm'
#   zOTUs
# raw data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG13890/Fusarium_spp/OG13890_zOTUs_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG13890_Fusarium_spp_zOTUs'
# Normalised data
# df1 <- read_delim("~/Downloads/AHDB/ampseq/OG13890/Fusarium_spp/OG13890_zOTUs_table_norm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# prefix <- 'OG13890_Fusarium_spp_zOTUs_norm'


# Create a varible that represents the split first column
x <- strsplit(df1$'#OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE)
# Use this to populate new columns for each element in this variable
df1$OTUs <- sapply( x, "[", 1 )
df1$Genus <- sapply( x, "[", 2 )
df1$Species <- sapply( x, "[", 3 )


#---
# Make a PCA plot of all data
#---
df2 <- df1
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
pca_plot <- autoplot(prcomp(df2), data = df3, colour = 'group')
#pca_plot <- pca_plot + theme(legend.title = element_blank())
pca_plot <- pca_plot + theme(legend.position="bottom", legend.title = element_blank())
pca_plot <- pca_plot + theme(legend.text = element_text(colour="black", size=7.5))
pca_plot <- pca_plot + theme(plot.margin=unit(c(1,2.5,0.5,0.5),"cm"))
filename <- paste(prefix, 'pca.pdf', sep='_')
ggsave(filename, plot = pca_plot, width =20, height = 15, units = "cm", limitsize = FALSE)

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

library(reshape2)
df4 <- melt(df1, id.vars=c('#OTU ID',"OTUs", 'Genus', 'Species'),value.name = "Counts", variable.name='Run')
df4$'#OTU ID' <- NULL
df4$group <- gsub("_.$","",df4$Run)
df4$experiment <- gsub("_.$|_pooled_reps|_pool","",df4$Run)
df4$reps <- ifelse(grepl("pool",df4$Run), "pooled","tech. reps.")
df4$reps <- as.factor(df4$reps)
df4$reps <- factor(df4$reps,levels(df4$reps)[c(2,1)])
# Generate SE for experimental reps
df5 <- summarySE(df4, measurevar="Counts", groupvars=c("OTUs", 'Genus', 'Species', 'group', 'experiment', 'reps'))
df5$OTU_label <- paste(df5$Species, df5$OTUs, sep = ' ')

# Create plot
facet_OTUs<-ggplot(data=subset(df5), aes(x=OTU_label, y=Counts))
facet_OTUs <- facet_OTUs + geom_bar(stat="identity")
facet_OTUs <- facet_OTUs + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_OTUs <- facet_OTUs + ylab('Mapped reads') + xlab('')
facet_OTUs <- facet_OTUs + geom_errorbar(aes(ymin=Counts-se, ymax=Counts+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_OTUs <- facet_OTUs + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
facet_OTUs <- facet_OTUs + facet_grid(experiment ~ reps)
facet_OTUs
filename <- paste(prefix, 'facet_OTUs.pdf', sep='_')
ggsave(filename, plot = facet_OTUs, width =25, height = 25, units = "cm", limitsize = FALSE)

#---
# Create plots combining data of OTUs with duplicate taxa designations
#---

df6 <- aggregate(df4$Counts, by=list(df4$Run, df4$Species), FUN=sum)
colnames(df6) <- c('Run', 'Species', 'Counts')
df6$experiment <- gsub("_.$|_pooled_reps|_pool","",df6$Run)
df6$reps <- ifelse(grepl("pool",df6$Run), "pooled","tech. reps.")
df6$reps <- as.factor(df6$reps)
df6$reps <- factor(df6$reps,levels(df6$reps)[c(2,1)])
# Generate SE for experimental reps
df7 <- summarySE(df6, measurevar="Counts", groupvars=c('Species', 'experiment', 'reps'))
df7$OTU_label <- paste(df7$Species, df7$OTUs, sep = ' ')
# Create plot
facet_species<-ggplot(data=subset(df7), aes(x=OTU_label, y=Counts))
facet_species <- facet_species + geom_bar(stat="identity")
facet_species <- facet_species + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_species <- facet_species + ylab('Mapped reads') + xlab('')
facet_species <- facet_species + geom_errorbar(aes(ymin=Counts-se, ymax=Counts+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_species <- facet_species + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
facet_species <- facet_species + facet_grid(experiment ~ reps)
facet_species
filename <- paste(prefix, 'facet_species.pdf', sep='_')
ggsave(filename, plot = facet_species, width =25, height = 25, units = "cm", limitsize = FALSE)
