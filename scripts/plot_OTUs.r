#!/usr/bin/Rscript


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

library(ggplot2)
library(readr)
library(optparse)

#---
# Read in data
#---
opt_list = list(
    make_option("--OTU_table", type="character", help="OTUs table"),
    make_option("--prefix", type="character", help="Output prefix")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$OTU_table
prefix = opt$prefix

df1 <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)

#---
# Create Facetplots summarising data by experiment and pooled / technical reps
#---

library(reshape2)
df4 <- melt(df1, id.vars=c('Species'),value.name = "Counts", variable.name='Run')
df4$Run <- gsub("_ITS$|_TEF$|_SIX13$|_OG12981$|_OG13890$","",df4$Run)
df4$group <- gsub("_.$","",df4$Run)
df4$experiment <- gsub("_.$|_pooled_reps|_pool","",df4$Run)
#df4$reps <- ifelse(grepl("pool",df4$Run), "pooled","tech. reps.")
df4$reps <- ifelse(grepl("pooled",df4$Run), "pooled","tech. reps.")
df4$reps <- as.factor(df4$reps)
df4$reps <- factor(df4$reps,levels(df4$reps)[c(2,1)])
# Generate SE for experimental reps
df5 <- summarySE(df4, measurevar="Counts", groupvars=c('Species', 'group', 'experiment', 'reps'))
df5$OTU_label <- paste(df5$Species, df5$OTUs, sep = ' ')

#Our transformation function
# labelFUN <- function(x) { if (sprintf("%.0f", x) < 206) {x <- sprintf("%.0f", x)} else { x <- sprintf("")}}
# df5$count_label <- ifelse(sprintf("%.0f", df5$Counts) < 206, NA, df5$Counts)
# df5$count_label[ df5$count_label<206 ] <- NULL
# LabelFun <- function(x) sprintf("%.0f", x)

# Create plot
facet_OTUs<-ggplot(data=subset(df5), aes(x=OTU_label, y=Counts))
facet_OTUs <- facet_OTUs + geom_bar(stat="identity")
facet_OTUs <- facet_OTUs + theme(axis.text.x=element_text(angle = -45, hjust = 0))
facet_OTUs <- facet_OTUs + ylab('Mapped reads') + xlab('')
# facet_OTUs <- facet_OTUs + geom_text(aes(label=df5$count_label, vjust=0))
# facet_OTUs <- facet_OTUs + geom_text(aes(label=LabelFun(Counts)), vjust=0, size=1)
facet_OTUs <- facet_OTUs + geom_errorbar(aes(ymin=Counts-se, ymax=Counts+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
facet_OTUs <- facet_OTUs + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
facet_OTUs <- facet_OTUs + facet_grid(experiment ~ reps)
facet_OTUs
filename <- paste(prefix, 'facet_species.pdf', sep='_')
ggsave(filename, plot = facet_OTUs, width =25, height = 25, units = "cm", limitsize = FALSE)
