library(ggplot2)
library(scales)

#############################################################################################
# @author:  Alexander Graf
# 
# Parameters:
#   outputPDF:  name of the output PDF file
#   threshold:  minimal coverage for a gene to be expressed [e.g. 10]
#   files:      tuples of file name of the HTSeq-count file and corresponding sample name
#               e.g. FILE_NAME,SAMPLE_NAME FILE_NAME2,SAMPLE_NAME2
#
#############################################################################################

args <- commandArgs(trailingOnly = TRUE)

outputPDF <- args[1]
threshold <- as.integer(args[2])
files <- c()
for(i in 3:length(args)){
        files <- c(files,args[i])
}

getTranscriptsOverThres <- function(data, i, threshold, sum, million){
        transcripts <- data[,i+1]
        transcripts <- transcripts*million/sum
        numberOfTrans <- length(transcripts[transcripts > threshold])
        return(numberOfTrans)
}

data <- data.frame()
summedReads <- data.frame()
plotTable <- data.frame(counts=integer(0), transcripts=integer(0), sample=character(0))

### read and sum #####
for (i in 1:length(files)){
        dataFile <- strsplit( files[i], "\\,")[[1]]
        sample <- dataFile[2]
        htseq <- read.table(dataFile[1], sep="\t", header=F)
        htseq <- htseq[1:(nrow(htseq)-3) ,]
        if(nrow(summedReads) > 0){
                summedReads <- rbind(summedReads, data.frame(Sample=sample, sum=sum(htseq$V2)))
        }else{
                summedReads <- data.frame(Sample=sample, sum=sum(htseq$V2))
        }
        htseq <- htseq[1:(nrow(htseq)-2) ,]

        if(nrow(data) > 0){
                coln <- colnames(data)
                data <- cbind(data, counts=htseq$V2)
                colnames(data) <- c(coln,sample)
        }else{
                data <- data.frame(Gene=htseq$V1, counts=htseq$V2)
                colnames(data) <- c("Gene",sample)
        }
}

summedReads$sum <- summedReads$sum/1000000
xlabValue <- 'Million reads (Mb)'
for(i in 1:nrow(summedReads)){
        sum <- summedReads[i,]$sum
        sample <- as.character(summedReads[i,]$Sample)
        for(j in seq(0,sum,length.out=20)){
                aboveThres <- getTranscriptsOverThres(data, i,threshold, sum, j)
                plotTable <- rbind(plotTable, data.frame(counts=j , transcripts=aboveThres , sample=summedReads[i,]$Sample))
        }

}

p <- ggplot(plotTable, aes(x=counts,y=transcripts,colour=sample)) + geom_point(size=0.75) +  geom_line(size=.7, shape=16) +xlab(xlabValue) +ylab(paste("Transcripts above ", threshold, "x")) +
theme(  axis.ticks.y = element_line(colour="black"),
                axis.ticks.x = element_line(colour="black"),
                axis.text.x = element_text(size=10, colour="black"),
                #axis.ticks.y = element_line(colour="black"),
                axis.text.y = element_text(size=10, colour="black"))

pdf(outputPDF, width=10, height=7)
plot(p)
dev.off()
                       
