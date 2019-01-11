#!/usr/local/bin/Rscript

library(ggplot2)
library(readr)
# install.packages('Rfast')
library(Rfast)
#setwd('~/Downloads/AHDB/ampseq')
#---
# Read in data
#---
library(optparse)
opt_list = list(
    make_option("--OTU_table_onion", type="character", help="OTUs table for onion field"),
    make_option("--OTU_table_daffodil", type="character", help="OTUs table for daffodil field"),
    make_option("--OTU_table_stocks", type="character", help="OTUs table for stocks field"),
    make_option("--prefix", type="character", help="Output prefix")
)
opt = parse_args(OptionParser(option_list=opt_list))
f1 = opt$OTU_table_onion
f2 = opt$OTU_table_daffodil
f3 = opt$OTU_table_stocks
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
# --- Summary SE function end ---



# Read in onion table
df1_onion <- read_delim(f1, "\t", escape_double = FALSE, trim_ws = TRUE)
# df1_onion <- read.delim("~/Downloads/AHDB_new/plate3/quantified/mix-B/onion/mix-B_OG4952_zOTUs_zOTUs_table.txt")
colnames(df1_onion)[1] <- "OTU ID"
df1_onion$'OTU ID' <- as.character(df1_onion$'OTU ID')
# x <- strsplit(df1_onion$'OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE) # Create a varible that represents the split first column
x <- strsplit(df1_onion$'OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE) # Create a varible that represents the split first column
df1_onion$Species <- sapply( x, "[", 1)


# Read in daffodil table
df1_daffodil <- read_delim(f2, "\t", escape_double = FALSE, trim_ws = TRUE)
# df1_daffodil <- read.delim("~/Downloads/AHDB_new/plate3/quantified/mix-B/daffodil/mix-B_OG4952_zOTUs_zOTUs_table.txt")
colnames(df1_daffodil)[1] <- "OTU ID"
df1_daffodil$'OTU ID' <- as.character(df1_daffodil$'OTU ID')
x <- strsplit(df1_daffodil$'OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE) # Create a varible that represents the split first column
df1_daffodil$Species <- sapply( x, "[", 1)

# Read in stocks table
df1_stocks <- read_delim(f3, "\t", escape_double = FALSE, trim_ws = TRUE)
# df1_stocks <- read.delim("~/Downloads/AHDB_new/plate3/quantified/mix-B/stocks/mix-B_OG4952_zOTUs_zOTUs_table.txt")
colnames(df1_stocks)[1] <- "OTU ID"
df1_stocks$'OTU ID' <- as.character(df1_stocks$'OTU ID')
x <- strsplit(df1_stocks$'OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE) # Create a varible that represents the split first column
df1_stocks$Species <- sapply( x, "[", 1)

df1_total <- merge(x=df1_onion, y=df1_daffodil, by=c("OTU ID", "Species"), all=TRUE)
df1_total <- merge(x=df1_total, y=df1_stocks, by=c("OTU ID", "Species"), all=TRUE)
df1_total$Species <- gsub("Gibberella","Fusarium",df1_total$Species,ignore.case=F)
# df1_total$Species <- gsub("_[a-zA-Z0-9-]+?$","",df1_total$Species,ignore.case=F, perl=T)

# df1_total$Species
df1_total$Species <- gsub('F\\.',"Fusarium ",df1_total$Species,ignore.case=F)
# df1_total$Species <- gsub("oxysporum.*","oxysporum",df1_total$Species,ignore.case=F, perl=T)
df1_total$Species <- gsub("_"," ",df1_total$Species,ignore.case=F)
df1_total$Species <- gsub("fsp","f.sp.",df1_total$Species,ignore.case=F)
df1_total$Species <- gsub("/"," ",df1_total$Species,ignore.case=F)

# df1_total$Species <- as.factor(df1_total$Species)

#--
# By Species
#--

library(reshape2)
df2 <- melt(df1_total, id.vars=c('OTU ID', 'Species'),value.name = "Counts", variable.name='Run')
df2$Counts[is.na(df2$Counts)] <- 0
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
df5 <- df5[!grepl("pool", df5$Bed),]

# df5$norm <- ((df5$Counts / df5$Total) * 1000)
df6 <- summarySE(df5, measurevar="Counts", groupvars=c(c("Mix", "Field", "Locus", "Species")))

df6$id = numeric(nrow(df6))
for (i in 1:nrow(df6)){
   dat_temp <- df6[1:i,]
   df6[i,]$id <- nrow(dat_temp[dat_temp$Species == df6[i,]$Species,])
 }
dfx <- dcast(df6, id~Species, value = 'value', value.var = 'Counts')
dfx$id <- NULL
dfx <- data.frame(apply(dfx, 2, sort, decreasing=T))
y <- as.logical(dfx[1,] > 1)
dfx <- data.frame(t(dfx))
dfx$keep <- y
dfx$X1 <- NULL
dfx$X2 <- NULL
dfx$X3 <- NULL
dfx$Species <- rownames(dfx)
dfx$Species <- gsub("\\."," ",dfx$Species,ignore.case=F)
dfx$Species <- gsub("f sp ","f.sp.",dfx$Species,ignore.case=F)
# dfx$Species
# df6$Species

dfy <- merge(dfx, df6, by="Species", all = TRUE)

df7 <- dfy[ which(dfy$keep == TRUE),]

df7$Field <- factor(df7$Field, levels = c("onion", "daffodil", "stocks"))

facet_Species<-ggplot(data=subset(df7), aes(x=Species, y=Counts))
# facet_Species <- facet_Species + scale_y_continuous(limits = c(0, 1000))
facet_Species <- facet_Species + geom_bar(stat="identity")
facet_Species <- facet_Species + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_Species <- facet_Species + ylab('Total reads') + xlab('')
facet_Species <- facet_Species + geom_errorbar(aes(ymin=Counts-se, ymax=Counts+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_Species <- facet_Species + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
facet_Species <- facet_Species + facet_grid(Field ~ .)
# facet_Species <- facet_Species + geom_text(aes(label=round(Counts)), vjust=-2.5)
facet_Species <- facet_Species + geom_text(aes(label=round(Counts)),  position = position_stack(vjust = 0.5))
facet_Species

# prefix <- '/Users/armita/Downloads/AHDB_new/plate3/plots/16S_Fields'
fig_height <- (length(unique(df2$Field))*5)+5
filename <- paste(prefix, 'facet_Species.pdf', sep='_')
ggsave(filename, plot = facet_Species, width =15, height = fig_height, units = "cm", limitsize = FALSE)
