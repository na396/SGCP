
################################################################################# ezSGCP
ezSGCP <- function(expData, geneID, annotation_db, semilabel = TRUE,
            calib = FALSE, norm = TRUE, tom = TRUE,
            saveAdja = FALSE, adjaNameFile = "adjacency.Rdata",
            hm = "adjaHeatMap.png",
            kopt = NULL, method_k = NULL, f.GO = sum, f.conduct = min,
            maxIteration = 1e8, numberStart = 1000, eff.egs = TRUE,
            saveOrig = TRUE, n_egvec = 100, sil = FALSE,
            dir = c("over", "under"), onto = c("BP", "CC", "MF"),
            hgCut = NULL, condTest = TRUE,
            cutoff = NULL, percent = 0.10, stp = 0.01,
            model = "knn", kn = NULL){

    # performs A to Z

    if(is.null(geneID) && semilabel == TRUE){
        warning("geneIDs are required for semilabeling", call. = FALSE)
        message("semilabeling will not performed")
        semilabel <- FALSE}

    if(is.null(geneID)){
        geneID <- paste0(rep("gene", nrow(expData)), seq(1,nrow(expData)))}


    if(length(geneID) != nrow(expData)){
        stop("number of rows in expData must be the same as the number of genes in geneID")}

    if(!is.matrix(expData) && !is.data.frame(expData)){
        stop("expData must be either a dataframe or a matrix")}

    if(ncol(expData) > nrow(expData)){
        warning("number of samples larger than the genes!!! \n Do rows correspond to genes?!"
                , call. = FALSE)}

    if(!is.null(kopt) && kopt != round(kopt) ){
        warning("kopt must be either null or an integer")
        message("making k null")
        kopt <- NULL}

    if(length(setdiff(method_k, c("relativeGap", "secondOrderGap", "additiveGap"))) != 0){
        warning("method_k can be either relativeGap, secondOrderGap, or additiveGap",
                call. = FALSE)
        message("making method to NULL")
        method_k <- NULL }

    if(!is.numeric(maxIteration) || !is.numeric(numberStart)){
        warning("maxIteration and numStart must be numeric and integer")
        message("making maxIteration and numberStart to default")
        maxIteration <- 1e+8
        numberStart <- 1000}

    if(all(dir %!in% c("under", "over"))){
        warning("dir must be in c(under or over) \n making to default",
                call. = FALSE)
        dir <- c("over", "under")}

    if(length(dir) > 2 ){
        warning("dir must be in c(under or over) \n making to default",
                call. = FALSE)
        dir <- c("over", "under")}

    if(all(onto %!in% c("BP", "CC", "MF"))){
        warning(" onto must be in BP CC MF \n making to default", call. = FALSE)
        onto <- c("BP", "CC", "MF")}

    if(length(onto) > 3 ){
        warning(" onto must be in BP CC MF \n making to default", call. = FALSE)
        onto <- c("BP", "CC", "MF")}

    if(!is.null(hgCut) && (hgCut >=1 || hgCut <= 0)){
        warning(" not correct hgCutoff value \n making to default", call. = FALSE) }

    if(condTest != TRUE && condTest != FALSE){
        warning(" condTest must be boolean! \n making to deafult", call. = FALSE)
        condTest <- TRUE}

    if(percent >= 1 || percent <= 0){
        warning("percent must be in (0,1) \n making percent to default",
                call. = FALSE)
        percent <- 0.10 }

    if(stp >= 1 || stp <= 0){
        warning("step must be in (0,1) \n making stp to default",
                call. = FALSE)
        stp <- 0.01 }

    if(!is.null(model) && model != "knn" & model != "lr"){
        warning("model must be either NULL, knn, or lr \n setting to knn", call. = FALSE)
        model <- "knn" }

    if(semilabel == TRUE && is.null(annotation_db)){
        warning("semilabeling need annotation_db \n semilabeling will not be performed",
                call. = FALSE)
        message("making semilabel to FALSE")
        semilabel <- FALSE}

    if(!is.null(semilabel) && semilabel != TRUE && semilabel != FALSE){
        warning("semilabel must be either TRUE or FALSE", call. = FALSE)
        message("making semilabel to FALSE")
        semilabel <- FALSE }

    newList <- list()
    ############################################################################### adjacencyMatrix
    message("starting network construction step...")
    resAdja <- adjacencyMatrix(expData = expData,
                        calibration = calib, norm = norm, tom = tom,
                        saveAdja = saveAdja, adjaNameFile = adjaNameFile,
                        hm = hm)


    ############################################################################### clustering
    message("starting network clustering step...")
    resClus <- clustering(adjaMat = resAdja, geneID = geneID,
                    annotation_db = annotation_db,
                    kopt = kopt, method = method_k,
                    func.GO = f.GO, func.conduct = f.conduct,
                    maxIter = maxIteration, numStart = numberStart,
                    eff.egs = eff.egs,
                    saveOrig = saveOrig, n_egvec = n_egvec, sil = sil)

    geneID <- resClus$geneID

    if(length(resClus$dropped.indices) > 0){
        expData <- expData[-resClus$dropped.indices, ] }


    rm(resAdja, expData)
    newList <- c(newList, setNames(list(resClus), "clustering"))
    ############################################################################### geneOntology

    if(!semilabel){

        clusterLabels <- data.frame(geneID = resClus$geneID,
                                initialClusters = resClus$clusterLabels)
        newList <- c(setNames(list(clusterLabels), "clusterLabels"), newList)


    }else{

        message("starting initial gene ontology enrichment step...")
        resGO <- geneOntology(geneUniv = geneID, clusLab = resClus$clusterLabels,
                        annotation_db = annotation_db,
                        direction = dir, ontology = onto, hgCutoff = hgCut,
                        cond = condTest)

        newList <- c(newList, setNames(list(resGO), "initial.GO" ))
        ############################################################################### semiLabeling

        message("starting semi-labeling stage...")
        resSL <- semiLabeling(geneID = geneID, df_GO = resGO$GOresults,
                        GOgenes = resGO$FinalGOTermGenes ,
                        cutoff = cutoff, percent = percent, stp = stp )

        newList <- c(newList, setNames(list(resSL), "semiLabeling" ))
        ################################################################################ semiSupervised

        message("starting semi-supervised step...")
        resSup <- semiSupervised(specExp = resClus$Y, geneLab = resSL$geneLabel,
                            model = model, kn = kn)

        newList <- c(newList, setNames(list(resSup), "semiSupervised"))
        ################################################################################ geneOntology
        geneLabel <- resSup$FinalLabeling

        message("starting final gene ontology enrichment step...")
        resFinal <- geneOntology(geneUniv = geneLabel$geneID,
                            clusLab = geneLabel$FinalLabel,
                            annotation_db = annotation_db,
                            direction = dir, ontology = onto, cond = condTest)

        newList <- c(newList, setNames(list(resFinal), "final.GO"))

        initial <- data.frame(geneID = resClus$geneID,
                        initialClusters = resClus$clusterLabels)
        final <- data.frame(geneID = geneID,
                        finalClusters = resSup$FinalLabeling$FinalLabel)
        clusterLabels <- inner_join(initial, final, by= "geneID")

        newList <- c(setNames(list(clusterLabels), "clusterLabels"), newList)


        remain <- setdiff(unique(clusterLabels$initialClusters),
                        unique(clusterLabels$finalClusters))

        if(length(remain)== 1){
            temp <- paste0("cluster ", remain, " is wiped out")
            message(temp)
        }else if(length(remain)){
            temp <- paste0("clusters", remain, " are wiped out")
            message(temp)}




    } # end of semiLabel

    message("ezSGCP done!..\n")
    newList <- c(setNames(list(semilabel), "semilabel"), newList)
    return(newList)
}








