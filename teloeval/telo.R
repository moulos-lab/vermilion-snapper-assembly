library(Biostrings)
library(rmarkdown)
library(rmdformats)
library(ggplot2)
library(scales)

# n: how many chromosomes to plot? 0 for all
teloEval <- function(fa,pat=c("TTAGGG","CCCTAA"),refpos=c("start","end"),n=0,
    path=NULL,template="telo.Rmd") {
    if (!file.exists(fa))
        stop("A FASTA file must be provided!")
    if (!all(is.character(pat)))
        stop("pat must be a character vector!")
    refpos <- refpos[1]
    
    message("Reading ",fa)
    S <- readDNAStringSet(fa)
    matchList <- vector("list",length(pat))
    names(matchList) <- pat
    for (i in seq_along(names(matchList))) {
        message("Checking pattern ",pat[i])
        tmp <- vmatchPattern(pat[i],S)
        if (refpos == "start")
            matchList[[i]] <- start(tmp)
        else
            matchList[[i]] <- end(tmp)
    }
    
    if (n>0) {
        for (i in seq_along(names(matchList)))
            matchList[[i]] <- matchList[[i]][seq_len(n)]
    }
    
    if (is.null(path))
        path <- paste0("teloeval_",format(Sys.time(),"%Y%m%d%H%M%S"))
    if (!dir.exists(path))
        dir.create(path,showWarnings=FALSE,recursive=TRUE)
    
    REP_ENV <- new.env(parent=globalenv())
    REP_ENV$matchList <- matchList
    
    #file.copy(#from="/media/raid/software/prisma/inst/eval_report.Rmd",
    #    from=system.file(package="prisma","eval_report.Rmd"),
    #    to=file.path(path,"eval_report.Rmd"),overwrite=TRUE)
    render(
        input=template,
        output_file="index.html",
        output_dir=path,
        envir=REP_ENV
    )
}



#~ g <- ggplot(x,aes(x=Position)) + 
#~  geom_histogram(binwidth=500000,color="darkblue",fill="lightblue") +
#~  xlab("\nPosition (bp)") +
#~  ylab(paste(pat[1]," occurences\n")) +
#~  ggtitle(nam) + 
#~  scale_x_continuous(labels=scales::label_number(
#~      scale_cut=cut_short_scale()),n.breaks=10) +
#~  theme(
#~      plot.title=element_text(hjust=0.5),
#~      axis.title.x=element_text(size=14,face="bold"),
#~      axis.title.y=element_text(size=14,face="bold"),
#~      axis.text.x=element_text(size=10,face="bold"),
#~      axis.text.y=element_text(size=12,face="bold")
#~  )

#~ g <- ggplot(x,aes(x=Position)) + 
#~  geom_density(color="red4",fill="red3",alpha = 0.2) +
#~  xlab("\nPosition (bp)") +
#~  ylab(paste(pat[1]," density\n")) +
#~  ggtitle(nam) + 
#~  scale_x_continuous(labels=scales::label_number(
#~      scale_cut=cut_short_scale()),n.breaks=10) +
#~  scale_y_continuous(label=scientific_10) +
#~  theme(
#~      plot.title=element_text(hjust=0.5),
#~      axis.title.x=element_text(size=14,face="bold"),
#~      axis.title.y=element_text(size=14,face="bold"),
#~      axis.text.x=element_text(size=10,face="bold"),
#~      axis.text.y=element_text(size=12,face="bold")
#~  )
