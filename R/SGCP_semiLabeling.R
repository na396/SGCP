
################################################################################# defineCutoff
defineCutoff <- function(df, percent){

    # defines the quantile based on given percent
    #
    # df: a dataframe returned by GOstats
    # percent: a number in (0,1) for percentile
    #
    # q: corresponding quantile

    pvalue <- df$Pvalue
    q <- quantile(pvalue, percent)
    return(q)
}

################################################################################# cvCutoff
cvCutoff <- function(df, percentile = 0.10, stp = 0.01){

    # find the cutoff, must hold one condition =>
    #                 number clusters for remaning genes must be greater equal to 2
    #
    # percentile: a number in (0,1)
    #          percentile of interest
    # stp: step for decreasing percentile a number of (0,1)
    #
    # cutoff: a number in (0,1) indicating the cutoff value for gene significancy

    flag <- TRUE
    while(flag){

        cutoff <- defineCutoff(df = df, percent = percentile)
        tdf <- df[df$Pvalue < cutoff, ]
        numClus <- length(unique(tdf$clusterNum))

        if(numClus <= 1){
            flag <- TRUE
            percentile <- percentile + stp
        }else{
            flag <- FALSE
        }

    }


    return(cutoff)
}

################################################################################# semiLabeling
semiLabeling <-
    function(geneID, df_GO, GOgenes, cutoff = NULL,
    percent = 0.10, stp = 0.01){


    if(is.null(geneID)){ stop("geneID is NULLL")}

    if(!is.data.frame(df_GO)){ stop("df_GO must be a dataframe")}

    if(!is.list(GOgenes)){ stop("GOgenes must be a list")}

    if(percent >= 1 || percent <= 0){
        warning("percent must be in (0,1) \n making percent to default",
                call. = FALSE)
        percent <- 0.10 }
    if(stp >= 1 || stp <= 0){
        warning("stp must be in (0,1) \n making percent to default",
                call. = FALSE)
        stp <- 0.01 }


    geneLabel <- data.frame(geneID = geneID, label = NA)
    totGenes <- nrow(geneLabel)

    if(is.null(cutoff)){
        cutoff <- cvCutoff(df_GO, percentile = percent, stp = stp) }


    message("cutoff value is ", cutoff)
    df_GO <- df_GO[df_GO$Pvalue < cutoff, ]

    clusterNums <- unique(df_GO$clusterNum)

    ## perform the semilabeling
    for(lab in clusterNums){

        sigGenes <- -1
        geneInClus <- c()
        GOIDs <- df_GO$GOID[df_GO$clusterNum == lab]

        caption <- paste0("Cluster", lab, "_GOTermGenes")
        GOTerms <- GOgenes[[caption]]
        for(go in GOIDs){
            sigGenes <- c(sigGenes, GOTerms[[go]])
            sigGenes <- unique(sigGenes)
        }

        geneLabel$label[geneLabel$geneID %in% sigGenes] <- lab

    }

    geneID <- data.frame(geneID = geneID)
    geneLabel <- inner_join(geneLabel, geneID, by = "geneID")

    newList <- list("cutoff" = cutoff, "geneLabel" = geneLabel)

    message("semiLabeling done!..\n")
    return(newList)

}
