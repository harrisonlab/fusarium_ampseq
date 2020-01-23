#!/usr/local/bin/Rscript

library(ggplot2)
library(readr)
# install.packages('Rfast')
library(Rfast)
library("dplyr")
library(reshape2)

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



# Read data:
read.delim2("/Volumes/GGB/Harrison/Projects/AHDB-FOC-2019-C300106/Science/Obj4-Establish-Potential-AmpliconSeq/metabarcoding_data/3 taxon abundance tables/FOC FON FOM soils/16S/16S_2019_samples_zOTUs_table.txt")

colnames(df1)[1] <- "OTU ID"
# Create a varible that represents the split first column
df1$'OTU ID' <- as.character(df1$'OTU ID')
x <- strsplit(df1$'OTU ID', ',', fixed = FALSE, perl = FALSE, useBytes = FALSE)
df1$Genus <- sapply( x, "[", 2)
df1$Species <- sapply( x, "[", 3)
# df1$Species <- gsub("_SH\\d+","",df1$Species,ignore.case=F)
# df1$Species <- gsub("_"," ",df1$Species,ignore.case=F)
# df1$Species = ifelse(df1$Species == 'unknown', df1$Genus, df1$Species)
# df1$Species <- gsub("Epicoccum pimprinum","Phoma pimprina",df1$Species,ignore.case=F)
# df1$Genus = ifelse(df1$Species == 'Phoma pimprina', "Phoma", df1$Genus)
# df1$Genus <- gsub("Gibberella","Fusarium",df1$Genus,ignore.case=F)
# df1$Genus <- gsub("/"," ",df1$Genus,ignore.case=F)
# df1$Genus <- gsub("\\-","",df1$Genus,ignore.case=F)
# df1$Species <- gsub("/"," ",df1$Species,ignore.case=F)
# df1$Species <- gsub("\\-","",df1$Species,ignore.case=F)
# 16S modifications:
df1$Genus <- gsub("Gp","Unkown group ",df1$Genus,ignore.case=F)
df1$Genus <- gsub("_"," ",df1$Genus,ignore.case=F)
df1$Genus <- gsub("genera incertae sedis","(incertae sedis)",df1$Genus,ignore.case=F)


df2 <- melt(df1, id.vars=c('OTU ID', 'Genus', 'Species'),value.name = "Counts", variable.name='Run')



df2$Counts[is.na(df2$Counts)] <- 0
df2$Run <- as.character(df2$Run)
y <- strsplit(df2$Run, '_', fixed = FALSE, perl = FALSE, useBytes = FALSE)
df2$Mix <- sapply( y, "[", 1)
df2$Field <- sapply( y, "[", 3)
df2$Dilution <- sapply( y, "[", 4)
df2$Rep <- sapply( y, "[", 5)
df2$Plate <- sapply( y, "[", 6)
df2$Locus <- sapply( y, "[", 7)

df3 <- aggregate(df2$Counts, by=list(df2$Genus, df2$Mix, df2$Field, df2$Dilution, df2$Rep, df2$Plate, df2$Locus), FUN=sum)
colnames(df3) <- c("Genus", "Mix", "Field", "Dilution", "Rep", "Plate", "Locus", "Counts")

# This step makes an additional dataframe with total reads per run
# This can then be used to normalise the data
df4 <- aggregate(df2$Counts, by=list(df2$Mix, df2$Field, df2$Dilution, df2$Rep, df2$Plate, df2$Locus), FUN=sum)
colnames(df4) <- c("Mix", "Field", "Dilution", "Rep", "Plate", "Locus", "Total")

df5 <- merge(df3,df4,by=c("Mix", "Field", "Dilution", "Rep", "Plate", "Locus"))
# df5 <- df5[!grepl("pool", df5$Dilution),]

df5$norm <- ((df5$Counts / df5$Total) * 1000)
df6 <- summarySE(df5[ which(df5$Total > 1000),], measurevar="norm", groupvars=c(c("Mix", "Field", "Locus", "Genus", "Plate")))

df6$id = numeric(nrow(df6))
for (i in 1:nrow(df6)){
   dat_temp <- df6[1:i,]
   df6[i,]$id <- nrow(dat_temp[dat_temp$Genus == df6[i,]$Genus,])
 }
# Create a data frame that can be used tofilter genera based on abundance
dfx <- dcast(df6, id~Genus, value = 'value', value.var = 'norm')
dfx$id <- NULL
dfx <- data.frame(apply(dfx, 2, sort, decreasing=T))
y <- as.logical(dfx[1,] > 10)
dfx <- data.frame(t(dfx))
dfx$keep <- y
dfx$X1 <- NULL
dfx$X2 <- NULL
dfx$X3 <- NULL
dfx$Genus <- rownames(dfx)

dfy <- merge(dfx, df6, by="Genus")

df7 <- dfy[ which(dfy$keep == TRUE),]

df7$Field <- gsub("FOC","onion",df7$Field,ignore.case=F)
df7$Field <- gsub("FON","daffodil",df7$Field,ignore.case=F)
df7$Field <- gsub("FOM","stocks",df7$Field,ignore.case=F)
df7$Field <- factor(df7$Field, levels = c("onion", "daffodil", "stocks"))

df7$Plate <- gsub("1","plate 1",df7$Plate,ignore.case=F)
df7$Plate <- gsub("2","plate 2",df7$Plate,ignore.case=F)

facet_Genus<-ggplot(data=subset(df7), aes(x=Genus, y=norm))
# facet_Genus <- facet_Genus + scale_y_continuous(limits = c(0, 1000))
facet_Genus <- facet_Genus + geom_bar(stat="identity")
facet_Genus <- facet_Genus + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_Genus <- facet_Genus + ylab('Read count (per 1000 mapped reads)') + xlab('')
facet_Genus <- facet_Genus + geom_errorbar(aes(ymin=norm-se, ymax=norm+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_Genus <- facet_Genus + theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
facet_Genus <- facet_Genus + facet_grid(Field ~ Plate)
# facet_Genus <- facet_Genus + geom_text(aes(label=round(norm)), vjust=-2.5)
facet_Genus <- facet_Genus + geom_text(aes(label=round(norm)),  position = position_stack(vjust = 0.5))
facet_Genus

# prefix <- '/Users/armita/Downloads/AHDB_new/plate3/plots/16S_Fields'
fig_height <- (length(unique(df2$Field))*5)+10
prefix <- "FOC_FOM_FON_16S"
filename <- paste(prefix, 'facet_genus.pdf', sep='_')
ggsave(filename, plot = facet_Genus, width =40, height = fig_height, units = "cm", limitsize = FALSE)




#
# # onion data by plate
# df1_onion <- select(df1,contains("FOC"))
# df1_onion_P1 <- select(df1_onion,contains("_1_ITS"))
# df1_onion_P2 <- select(df1_onion,contains("_2_ITS"))
# # daffodil data by plate
# df1_daffodil <- select(df1,contains("FON"))
# df1_daffodil_P1 <- select(df1_daffodil,contains("_1_ITS"))
# df1_daffodil_P2 <- select(df1_daffodil,contains("_2_ITS"))
# # stocks data by plate
# df1_stocks <- select(df1,contains("FOM"))
# df1_stocks_P1 <- select(df1_stocks,contains("_1_ITS"))
# df1_stocks_P2 <- select(df1_stocks,contains("_2_ITS"))
#
#
# df1_total <- merge(x=df1_onion, y=df1_daffodil, by=c("OTU ID", "Genus"), all=TRUE)
# df1_total <- merge(x=df1_total, y=df1_stocks, by=c("OTU ID", "Genus"), all=TRUE)
# df1_total$Genus <- gsub("Gibberella","Fusarium",df1_total$Genus,ignore.case=F)
