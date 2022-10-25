

################################################################################# silhouetteIndex
silhouette <- function(dis.y, clus.labels){

    # dis.y: a dist object
    #       as.dist(y)
    # clus.labels: a vector
    #             contains the clusterlabels for the point

    #checkSym(dis.y, stp = "silhouette index")

    if(!is.matrix(dis.y)){ dis.y <- as.matrix(dis.y)}

    if(length(clus.labels) != nrow(dis.y)){
        stop("number of elements in clus.labels and number of rows in dis.y are not equal")}

    sil <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(sil) <- c("geneIndices", "silIndex", "clusterLabel")

    clusters <- unique(clus.labels)


    for (inclus in clusters) {

        indw <- which(clus.labels == inclus)
        if(length(indw) > 1){

            a <- dis.y[indw, indw, drop = FALSE]

            clus.size.in <- max(2, ncol(a))
            A <- rowSums(a)/(clus.size.in - 1)
            A <- as.matrix(A)

            B <-  matrix(Inf, nrow = length(indw), ncol = 1)
            for(outclus in clusters){

                if(outclus != inclus){
                    indb <- which(clus.labels == outclus)
                    b <- dis.y[indw, indb, drop = FALSE]
                    clus.size.out <- max(1, ncol(b))
                    cur <- as.matrix(rowSums(b)/clus.size.out)
                    B <- pmin(B, cur)
                    B <- as.matrix(B)} # end of if(outclus != inclus)


            } # end of outclus for loop
            S <- (B - A)/pmax(A, B)
            tdf <- data.frame(geneIndices = indw, silIndex = S, clusterLabel = inclus)

        }else{

            tdf <- data.frame(geneIndices = indw, silIndex = 0, clusterLabel = inclus)

        } # end of index length
        sil <- rbind(sil, tdf)

    } # end of inclus for loop

    sil <- sil %>%  arrange(clusterLabel,-silIndex)
    sil$clusterLabel <- as.factor(sil$clusterLabel)

    return(sil)

} # end of function

################################################################################# divideNorm
divideNorm <- function(M, rowWise = TRUE){

    # Divide by Norm
    # Argument:
    #         M: a matrix
    #         rowWise: either TRUE or FALSE
    #               if TRUE it divides each row by its norm
    #               if FALSE it divides each column by its norm
    #
    # Value:
    #       M: the scaled matrix


    if(rowWise == TRUE){
        M <- t(M)
        M <- as.data.frame(lapply(as.data.frame(M), normalization))
        M <- as.matrix(M)
        M <- t(M)

    }else if(rowWise == FALSE){
        #M <- t(M)
        M <- as.data.frame(lapply(as.data.frame(M), normalization))
        M <- as.matrix(M)
        #M <- t(M)

    }else{ stop( " rowWise can be either TRUE or FALSE")}

    return(M)

}

################################################################################# dropNoise
findNoise <- function(evec, thresh = 20){

    # drops noise components using eigenvectors
    #
    # Arguments:
    #         evec: a matrix of eigenvectors, corresponds to columns
    #         thresh: an integer
    #                   minimum cluster size to be considered as noise
    #
    # Value:
    #         final_ind: index of noisy genes

    v2 <- evec[, 2]
    indp <- which(v2 >= 0)
    p <- length(indp)
    indn <- which(v2 < 0)
    n <- length(n)

    if(length(v2) < 2*thresh){
        temp <- paste0("not enough genes, total number of genes is ", length(v2))
        stop(temp)
    }

    if(p > thresh && n > thresh){
        message("network has at least to components and each is large")
        message("this package does not support more than one components")
        message("you need to run SGCP on each components by yourself")

    }else{
        if(p < thresh){
            final_ind <- indp

        }else if(n < thresh) {
            final_ind <- indn
        }

    }

    return(final_ind)
}

################################################################################# clusterPlots
clusterPlots <- function(plt, tit, xname, yname){

    plt <-  plt +
        theme_classic()  +
        theme(axis.text.x = element_text(angle= 45, hjust=1, size = 5, face = 'bold',lineheight = 0.9)) +
        theme(axis.text.y = element_text(size = 15, face = 'bold',lineheight = 0.9)) +
        theme(axis.title.y = element_text(size = 10, face = 'bold',lineheight = 0.9)) +
        theme(plot.title = element_text(size = 20, face = 'bold',lineheight = 0.9, hjust = 0.5)) +
        theme(legend.text = element_text(size = 7))+
        theme(legend.position = c(.9, .9)) +
        labs(x = xname, y = yname, title = tit)
    return(plt)

}

################################################################################# clusterNumber
clusterNumber <- function(egvals, maxNum = 102){

    # Finds the number of clusters using relative eigenvalues ratios
    #
    # Argument:
    #       egvals: a vector of eigenvalues, sorted non increasing
    #       maxNum: an integer, upper bound on the maegvalsimum number of cluster
    #
    # Value:
    #       numOpt: a list
    #             relativeGap: k for method relativeGap
    #             secondOrderGap: k for method secondOrderGap
    #             additiveGap: k for method additiveGap
    #             plt: plots for each methods


    if(!is.numeric(egvals)){ stop("eigenvalues egvals are not numeric") }
    if(is.unsorted(rev(egvals))){
        warning("egvals are not sorted", call. = FALSE)
        message("sorting the eigenva;ues...")
        egvals <- sort(egvals, decreasing = TRUE)}

    if(round(egvals[1], 7) == 1){egvals <- egvals[-1]}

    # method: eigen gap
    #k <- seq(2, maxNum)
    k <- seq_len(maxNum)
    k <- k[-1]
    pairDiff <- abs(diff(egvals[seq_len(maxNum)]))

    dfgap <- as.data.frame(cbind(k, pairDiff))
    #dfgap <- dfgap[-1, ]
    optg <- dfgap[which.max(pairDiff),]$k

    # method: first and second order
    egvals <- 1 - egvals
    egvals <- egvals[seq_len(min(maxNum, length(egvals)))]

    df <- data.frame(eigenVal = egvals, firstOrder = NA, kfirst = NA,
                    secondOrder = NA, ksecond = NA)

    i <- 2
    while(i <= nrow(df)){
        df$firstOrder[i] <- df$eigenVal[i]/df$eigenVal[i-1]
        df$kfirst[i] <- i
        i <- i + 1}

    i <- 3
    while (i <= nrow(df)) {
        df$secondOrder[i] <- df$firstOrder[i-1] - df$firstOrder[i]
        df$ksecond[i] <- i- 1
        i <- i+1}

    df$secondOrder <- c(NA, diff(df$firstOrder))
    #df$indices <- seq(1:nrow(df))
    #print(head(df))

    # method: first
    resf <- df$firstOrder
    resf <- resf[-1]
    optf <- which.max(resf) + 1

    # method: second
    ress <- df$secondOrder
    ress <- ress[-c(1,2)]
    opts <- which.max(ress) + 1

    numOpt <- list("relativeGap" = optf, "secondOrderGap" = opts, "additiveGap" = optg)
    temp <- paste0("number of clusters for relativeGap method is ")
    message(temp, optf)
    temp <- paste0("number of clusters for secondOrderGap method is ")
    message(temp, opts)
    temp <- paste0("number of clusters for additiveGap method is ")
    message(temp, optg)

    #if(plt){

    pltgap <- ggplot(data = dfgap, aes(x = k, y = pairDiff)) +
        geom_point(size = .5) +
        geom_line()

    pltgap <- clusterPlots(pltgap, tit = "Additive Gap Method",
                        xname = "cluster number", yname = "Additive Gap")

    pltfirst <- ggplot(data = df, aes(x = kfirst, y = firstOrder)) +
                    geom_point(size = .5) +
                    geom_line()

    pltfirst <- clusterPlots(pltfirst, tit = "Relative Gap Method",
                            xname = "cluster number", yname = "Relative Gap")

    pltsecond <- ggplot(data = df, aes(x = ksecond, y = secondOrder)) +
        geom_point(size = .5) +
        geom_line()

    pltsecond <- clusterPlots(pltsecond, tit = "Second-Order Gap Method",
                            xname = "cluster number", yname = "Second-Order Gap")

    plt <- list("additiveGap" = pltgap, "relativeGap" = pltfirst,
                "secondOrderGap" = pltsecond)
    numOpt <- c(numOpt, "plots" = plt)
    return(numOpt)
}


################################################################################# conductance
conductance <- function(adja, clusLab){

    # calculates the conductance of a clustering
    #
    # Arguments:
    #           adja: a matrix,
    #               symmetric sqaure adjacency matrix with values in (0,1)
    #           clusLab: a vector of cluster label
    #
    # Values:
    #       a list: mean, median, clusterConductance
    #               mean: mean over all conduct value per cluster
    #               median: median over all conduct value per cluster
    #               conductance: conductance value per cluster

    conduct <- c()
    deg <- matrix(rowSums(adja), ncol = 1)

    for(clus in unique(clusLab)){
        indin <- which(clusLab == clus)
        indout <- which(clusLab != clus)

        cur <- 0
        M <- adja[indin, indout]
        cur <- cur + sum(M)

        conduct <- c(conduct, cur/sum(deg[indin]))

    }
    names(conduct) <- unique(clusLab)

    newList <- list("mean" = mean(conduct), "median" = median(conduct),
                    "clusterConductance" = conduct)
}

################################################################################# cvConductance
cvConductance <- function(adja, k, Y, X, func = min,
                        maxIter = maxIter, numStart = numStart){

    # perform cross validation
    #
    # Arguments:
    #         adja: matrix
    #             symmetric squared with values in (0,1) matrix
    #         k: a list containing
    #           relativeGap: scalar containing k relativeGap method
    #           secondOrderGap: scalar containing k for secondOrderGap method
    #           additiveGap: scalar containing k for additiveGap method
    #           plots: plots of the methods
    #         Y: matrix of eigenvectors
    #         X: vector of eigenvalues
    #         method: either mean, median, or mean
    #                 for comparing conductance indices
    #         maxIter: a number,
    #                 maximum number of iteration in kmeans
    #         numStart: a number
    #                   number of start for kmeans
    #
    # Values:
    #       a list of
    #       method: string of method names
    #       k: scalar for number of cluster
    #       Y: a rectangular matrix of the transformation
    #       X: a vector of eigenvalues correspond to Y
    #       cluster: output of kmeans
    #       clusterLabels: a vector of datapoints labels

    message("Conductance Validation...")

    # first order method
    krelativeGap <- k$relativeGap
    Yf <- Y[, seq_len(min(2*krelativeGap, ncol(Y)))]
    Xf <- X[seq_len(min(2*krelativeGap, length(X)))]
    Yf <- divideNorm(Yf, rowWise = TRUE)
    clusf <- kmeans(Yf, krelativeGap, iter.max = maxIter, nstart = numStart)
    conf <- conductance(adja = adja, clusLab = clusf$cluster)

    # second order method
    ksecondOrderGap <- k$secondOrderGap
    Ys <- Y[, seq_len(min(2*ksecondOrderGap, ncol(Y)))]
    Xs <- X[seq_len(min(2*ksecondOrderGap, length(X)))]
    Ys <- divideNorm(Ys, rowWise = TRUE)
    cluss <- kmeans(Ys, ksecondOrderGap, iter.max = maxIter, nstart = numStart)
    cons <- conductance(adja = adja, clusLab = cluss$cluster)

    # eigengap method
    kadditiveGap <- k$additiveGap
    Yg <- Y[, seq_len(min(2*kadditiveGap, ncol(Y)))]
    Xg <- X[seq_len(min(2*kadditiveGap, length(X)))]
    Yg <- divideNorm(Yg, rowWise = TRUE)
    clusg <- kmeans(Yg, kadditiveGap, iter.max = maxIter, nstart = numStart)
    cong <- conductance(adja = adja, clusLab = clusg$cluster)

    res <- c(func(conf$clusterConductance),
            func(cons$clusterConductance),
            func(cong$clusterConductance))
    names(res) <- c("relativeGap", "secondOrderGap", "additiveGap")

    resmethod <- names(which.min(res))
    temp <- paste0("selected method in cvConductance is ", resmethod )
    message(temp)

    methodCon <- list("relativeGap.conductance" = conf,
                    "secondOrderGap.conductance" = cons,
                    "additiveGap.conductance" = cong)

    if(resmethod == "relativeGap"){
        method <- "relativeGap"
        k <- krelativeGap
        Y <- Yf
        X <- Xf
        clus <- clusf
        clusterLabels <- clus$cluster

    }else if(resmethod == "secondOrderGap"){
        method <- "secondOrderGap"
        k <- ksecondOrderGap
        Y <- Ys
        X <- Xs
        clus <- cluss
        clusterLabels <- clus$cluster


    }else{
        method <- "additiveGap"
        k <- kadditiveGap
        Y <- Yg
        X <- Xg
        clus <- clusg
        clusterLabels <- clus$cluster }

    newList <- list("method" = method, "k" = k, "Y" = Y, "X" = X,
                    "clus" = clus, "clusterLabels" = clusterLabels,
                    "conductance" = methodCon)

    return(newList)
}

################################################################################# sigClusGO
sigClusGO <- function(adja, Y, X, k, annotation, geneID,
                    maxIter = 1e+8, numStart = 1000){

    # performs geneontology for a single cluster produced by
    #           either relativeGap, secondOrderGap, and additiveGap
    #
    # Arguments:
    #         adja: a matrix
    #               symmetric square matrix with values in (0,1)
    #         Y: the eigenvectors
    #         X: eigenvalues
    #         k: number of clusters
    #         annotation: genome wide annotation
    #         geneID: a vector of gene IDs
    #         maxIter: a number, indicates maximum number of iteration for kmeans
    #         numStart: a number, indicates the number of start for kmeans
    #
    # Value
    #     Y: normalized Y
    #     X: the eigenvalues
    #     k: number of clusters
    #     clus: clustering object returns by kmeans function
    #     conductance: the result returns by conductance function
    #     pvalues: the pvalues associated to the cluster GO
    #     df_clus: the cluster df for GO
    #


    Yt <- Y[, seq_len(min(2*k, ncol(Y)))]
    Xt <- X[seq_len(min(2*k, length(X)))]
    Yt <- divideNorm(Yt, rowWise = TRUE)
    clust <- kmeans(Yt, k, iter.max = maxIter, nstart = numStart)
    cont <- conductance(adja = adja, clusLab = clust$cluster)

    lab.t <- names(which.min(cont$clusterConductance))
    ind.t <- which(clust$cluster == lab.t)
    genes.t <- geneID[ind.t]

    nl.t <- GOenrichment(universe_geneIDs = geneID,
                        cluster_geneIDs = genes.t,
                        annotation_db = annotation,
                        lab = lab.t,
                        direction = c("over", "under"),
                        ontology = c("BP", "CC", "MF"),
                        hgCutoff = NULL,
                        condition = TRUE)
    #df = summary(nl.t)

    newList <- list("Y" = Yt, "X" = Xt, "k" = k, "clus" = clust,
                    "conductance" = cont, "pvalues" = nl.t$GO_summary$Pvalue,
                    "df_clus" = nl.t$GO_summary)
    return(newList)
}

################################################################################# cvGO
cvGO <- function(adja, k, Y, X, annotation, geneID,
                maxIter = maxIter, numStart = numStart, func = sum){

    # perform cross validation
    #
    # Arguments:
    #         adja: a matrix
    #               symmetric square matrix with values in (0,1)
    #         Y: the eigenvectors
    #         X: eigenvalues
    #         k: number of clusters
    #         annotation: genome wide annotation
    #         geneID: a vector of gene IDs
    #         maxIter: a number, indicates maximum number of iteration for kmeans
    #         numStart: a number, indicates the number of start for kmeans
    #         func: function to perform cv for the GOs
    #
    # Values:
    #       a list of containing
    #       method: string of method names
    #       k: scalar for number of cluster
    #       Y: a rectangular matrix of the transformation
    #       X: a vector of eigenvalues correspond to Y
    #       cluster: output of kmeans
    #       clusterLabels: a vector of datapoints labels
    #       conductance for the best selected method

    message("\n Gene Ontology Validation...")

    # first order method
    message("method relativeGap....")
    nl.f <- sigClusGO(adja, Y, X, k$relativeGap, annotation, geneID,
                    maxIter = maxIter, numStart = numStart)
    f <- -log10(nl.f$pvalues)

    if(nrow(nl.f$df_clus) >0 ){
        nl.f$df_clus$clusterNum <- "relativeGap"}

    df.f <- nl.f$df_clus

    # second order method

    message("\n method secondOrderGap...")
    if(k$secondOrderGap != k$relativeGap){

        nl.s <- sigClusGO(adja, Y, X, k$secondOrderGap, annotation, geneID,
                        maxIter = maxIter, numStart = numStart)
        s <- -log10(nl.s$pvalues)

        if(nrow(nl.s$df_clus) >0 ){
            nl.s$df_clus$clusterNum <- "secondOrderGap" }

        df.s <- nl.s$df_clus

    }else {

        message("GO enrichment is the same is relativeGap!")
        nl.s <- nl.f
        s <- nl.s$pvalues

        if(nrow(nl.s$df_clus) >0 ){
            nl.s$df_clus$clusterNum <- "secondOrderGap" }

        df.s <- nl.s$df_clus}


    message("\n method additiveGap....")
    if(k$additiveGap == k$relativeGap){
        message("GO enrichment is the same is relativeGap!")

        nl.g <- nl.f
        g <- -log10(nl.g$pvalues)

        if(nrow(nl.g$df_clus) >0 ){
            nl.g$df_clus$clusterNum <- "additiveGap"}

        df.g <- nl.g$df_clus

    }else if(k$additiveGap == k$secondOrderGap){
        message("GO enrichment is the same is secondOrderGap!")

        nl.g <- nl.s
        g <- -log10(nl.g$pvalues)

        if(nrow(nl.g$df_clus) >0 ){
            nl.g$df_clus$clusterNum <- "additiveGap"}

        df.g <- nl.g$df_clus

    }else{
        nl.g <- sigClusGO(adja, Y, X, k$additiveGap, annotation, geneID,
                        maxIter = maxIter, numStart = numStart)
        g <- -log10(nl.g$pvalues)

        if(nrow(nl.g$df_clus) >0 ){
            nl.g$df_clus$clusterNum <- "additiveGap"}

        df.g <- nl.g$df_clus}


    methodCon <- list("relativeGap.conductance" = nl.f$conductance,
                    "secondOrderGap.conductance" = nl.s$conductance,
                    "additiveGap.conductance" = nl.g$conductance)


    # taking cv function over p-values
    res <- c(func(f), func(s), func(g))
    names(res) <- c("relativeGap", "secondOrderGap", "additiveGap")
    resmethod <- names(which.max(res))


    if(resmethod == "relativeGap"){
        method <- "relativeGap"
        k <- nl.f$k
        Y <- nl.f$Y
        X <- nl.f$X
        clus <- nl.f$clus
        clusterLabels <- clus$cluster
        conduct <- nl.f$conductance
        df <- nl.f$df_clus

    }else if(resmethod == "secondOrderGap"){
        method <- "secondOrderGap"
        k <- nl.s$k
        Y <- nl.s$Y
        X <- nl.s$X
        clus <- nl.s$clus
        clusterLabels <- clus$cluster
        conduct <- nl.s$conductance
        df <- nl.s$df_clus

    }else{
        method <- "additiveGap"
        k <- nl.g$k
        Y <- nl.g$Y
        X <- nl.g$X
        clus <- nl.g$clus
        clusterLabels <- clus$cluster
        conduct <- nl.g$conductance
        df <- nl.g$df_clus

    }

    df <- rbind(df.f, df.s, df.g)
    newList <- list("method" = method, "k" = k, "Y" = Y, "X" = X,
                    "clus" = clus, "clusterLabels" = clusterLabels,
                    "conductance" = methodCon, "df" = df)

    return(newList)
}


################################################################################# transformation
clustering <- function(adjaMat, geneID , annotation_db ,
                    kopt = NULL, method = NULL,
                    func.GO = sum, func.conduct = min,
                    maxIter = 1e8, numStart = 1000,
                    saveOrig = TRUE, n_egvec = 100, sil = FALSE){


    if(is.data.frame(adjaMat)){
        warning("adjaMat is a dataframe", call. = FALSE)
        message("converting adjancecy to matrix")
        df2mat(adjaMat)}


    if(!is.matrix(adjaMat)){ stop("input must be a matrix") }
    checkNumeric(adjaMat, " transformation input")
    checkSym(adjaMat, " transformation input")

    ng <- nrow(adjaMat)

    if(!is.null(geneID)){
        if(length(geneID) != ng){
            warning("length of geneID must be equal to number of rows (or columns) in adjaMat",
                    call. = FALSE)
            message("setting geneIDs to NULL")
            geneID <- NULL} }

    if(is.null(geneID)){
        geneID <- paste0(rep("gene", ng), seq_len(ng))}

    if(!is.null(annotation_db) & !is.character(annotation_db)){
        stop("annotation_db must be character or NULL", call. = FALSE)}

    if(!is.character(geneID)){
        warning("type of geneID is not character")
        message("makig geneIDs type to character")
        geneID <- as.character(geneID) }

    if(!is.null(kopt) && kopt != round(kopt)){
        warning("kopt must be either null or an integer", call. = FALSE)
        message("making kopt null")
        kopt <- NULL}

    if(length(setdiff(method, c("relativeGap",
                                "secondOrderGap", "additiveGap"))) != 0){
        warning("method can be either relativeGap, secondOrderGag, or additiveGap",
                call. = FALSE)
        message("making method to NULL")
        method <- NULL }

    if(!is.numeric(maxIter) || !is.numeric(numStart)){
        warning("maxIter and numStart must be numeric and integer", call. = FALSE)
        message("making maxIter and numStart to default")
        maxIter <- 1e+8
        numStart <- 1000}

    if(is.numeric(n_egvec) && n_egvec != round(n_egvec)){
        warning("n_egvec must be an integer or all", call. = FALSE) }


    # Calculate the diagonal degree matrix
    D <- diag(rowSums(adjaMat))
    d <- as.matrix(as.vector(rowSums(adjaMat)))

    # Drop the zero values for D and drop that gene from the network
    ind <- which(d == 0)
    if(length(ind) > 0){
        adjaMat <- adjaMat[-ind, ]
        adjaMat <- adjaMat[, -ind]
        D <- D[-ind, ]
        D <- D[, -ind]
        if(!is.null(geneID)){geneID <- geneID[-ind] } }

    # Calculate the diagonal degree matrix
    D <- diag(rowSums(adjaMat))
    d <- as.matrix(as.vector(rowSums(adjaMat)))

    # Normalized Laplacian
    message("calculating normalized Laplacian \n it may take time...")
    diag(D) <- 1/sqrt(diag(D))
    L <- D %*% adjaMat %*% D

    # Calculate the eigenvalues and eigenvectors
    message("calculating eigenvalues/vectors \n it may take time...")
    egvv <- eigen(L)
    Y <- egvv$vectors
    X <- egvv$values

    if(!is.unsorted(X)){
        warning(" The eigenvalues are not sorted!", call. = FALSE)
        message("sorting the eigenvalues/vectors...")
        eg <- as.matrix(rbind(X, Y))
        eg <- as.data.frame(t(eg))
        colnames(eg)[1] <- "eigenvalues"
        eg <- eg[order(eg$eigenvalues, decreasing = TRUE), ]
        eg <- t(as.matrix(eg))
        X <- eg[1, ]
        temp <- seq_len(nrow(eg))
        temp <- temp[-1]
        Y <- eg[temp, ]}

    # preprocessing to check for noisy genes
    if(sum(round(X, 7) == 1) > 1 ){
        message("dropping noisy genes...")
        nois_ind <- findNoise(evec = Y, thresh = 100)
        adjaMat <- adjaMat[-nois_ind, ]
        adjaMat <- adjaMat[, -nois_ind]
        D <- D[-nois_ind, ]
        D <- D[, -nois_ind]
        ind <- c(ind, nois_ind)

        if(!is.null(geneID)){geneID <- geneID[-nois_ind]
        temp <- paste0("number of noisy genes are... ", length(nois_ind))
        message(temp)
        message("calculation from begining...")}

        # Calculate the diagonal degree matrix
        D <- diag(rowSums(adjaMat))
        d <- as.matrix(as.vector(rowSums(adjaMat)))

        # Normalized Laplacian
        message("calculating normalized Laplacian \n it may take time...")
        diag(D) <- 1/sqrt(diag(D))
        L <- D %*% adjaMat %*% D

        # Calculate the eigenvalues and eigenvectors
        message("calculating eigenvalues/vectors \n it may take time...")
        egvv <- eigen(L)
        Y <- egvv$vectors
        X <- egvv$values }

    #
    Y <- D %*% Y
    Y <- sweep(Y, 2, (t(d) %*% Y)/sum(d))

    #
    D <- diag(rowSums(adjaMat))
    scalar <- function(x, D){return(x/sqrt(as.vector(t(x) %*% D %*% x)))}
    Y <- apply(Y, 2, scalar, D)

    constant <- t(Y[,1]) %*% D %*% Y[, 1]
    if(round(constant, 7) != 1 ){
        warning("eigenvector 1 condition does not hold! ", call. = FALSE)
        temp <- paste0("the condition value is ",
                    as.double(t(Y[,1]) %*% D %*% Y[, 1] != 1))
        message(temp) }


    # Drop first eigenvalue and eigenvectors
    X <- X[-1]
    Y <- Y[, -1]

    # saving original

    if(saveOrig == TRUE){
        if(n_egvec == "all"){
            n_egvec <- ncol(Y)

        }else if(is.character(n_egvec)){
            message("n_egvec is not neigher all nor integer")
            message("setting n_egvect to 110")
            n_egvec <- 110

        }else if(n_egvec != round(n_egvec)){
            warning("n_egvec must be an integer", call. = FALSE)
            message("making n_egvec to 110")
            n_egvec <- 110
        }

        temp <- paste0("n_egvec is ", n_egvec)
        message(temp)
        Yorig <- Y[, seq_len(min(n_egvec, ncol(Y)))]
        Xorig <- X[seq_len(min(n_egvec, length(X)))]
        Yorig <- divideNorm(Yorig, rowWise = TRUE)

    }

    plt <- list()

    #################################################################################

    if(is.null(kopt) && is.null(method) && !is.null(annotation_db)){


        k <- clusterNumber(egvals = X, maxNum = 102)
        plt <- c(plt, setNames(list(k$plots.relativeGap), "relativeGap"),
                    setNames(list(k$plots.secondOrderGap), "secondOrderGap" ),
                    setNames(list(k$plots.additiveGap), "additiveGap") )


        cvopt <- cvGO(adja = adjaMat, k = k, Y = Y, X = X,
                    annotation = annotation_db, geneID = geneID,
                    maxIter = maxIter, numStart = numStart, func = func.GO)

        method <- cvopt$method
        k <- cvopt$k
        Y <- cvopt$Y
        X <- cvopt$X
        clus <- cvopt$clus
        clusterLabels <- cvopt$clusterLabels
        conduct <- cvopt$conductance
        df <- cvopt$df
        cv <- "cvGO"

        temp <- paste0("\n method ", method, " is selected using GO validation and k is ", k)
        message(temp)

        newList <- list("dropped.indices" = ind, "geneID" = geneID, "method" = method,
                        "k" = k, "Y" = Y, "X" = X, "cluster" = clus,
                        "clusterLabels" = as.character(clusterLabels),
                        "conductance" = conduct, "cvGOdf" = df, "cv" = cv,
                        "clusterNumberPlot" = plt)

    } else if(is.null(kopt) && is.null(method) && is.null(annotation_db)){

        k <- clusterNumber(egvals = X, maxNum = 102)
        plt <- c(plt, setNames(list(k$plots.relativeGap), "reltiveGap"),
                    setNames(list(k$plots.secondOrderGap), "secondOrderGap" ),
                    setNames(list(k$plots.additiveGap), "additiveGap") )

        cvopt <- cvConductance(adja = adjaMat, k = k, Y = Y, X = X,
                            func = func.conduct,
                            maxIter = maxIter, numStart = numStart)

        method <- cvopt$method
        k <- cvopt$k
        Y <- cvopt$Y
        X <- cvopt$X
        clus <- cvopt$clus
        clusterLabels <- cvopt$clusterLabels
        conduct <- cvopt$conductance
        cv <- "cvConductance"

        temp <- paste0("\n method ", method, " is selected using conductance index validation and k is ", k)
        message(temp)

        newList <- list("dropped.indices" = ind, "geneID" = geneID, "method" = method,
                        "k" = k, "Y" = Y, "X" = X, "cluster" = clus,
                        "clusterLabels" = as.character(clusterLabels),
                        "conductance" = conduct, "cv" = cv,
                        "clusterNumberPlot" = plt)

    }else if(is.null(kopt) && !is.null(method)){
        k <- clusterNumber(egvals = X, maxNum = 102)

        if(method == "relativeGap"){

            # first order method
            krelativeGap <- k$relativeGap
            Yf <- Y[, seq_len(min(2*krelativeGap, ncol(Y)))]
            Xf <- X[seq_len(min(2*krelativeGap, length(X)))]
            Yf <- divideNorm(Yf, rowWise = TRUE)

            clusf <- kmeans(Yf, krelativeGap, iter.max = maxIter, nstart = numStart)
            conf <- conductance(adja = adjaMat, clusLab = clusf$cluster)

            method <- "relativeGap"
            k <- krelativeGap
            Y <- Yf
            X <- Xf
            clus <- clusf
            clusterLabels <- clus$cluster
            conduct <- conf
            cv <- "userMethod"

            temp <- paste0("\n method ", method, " is selected using user and k is ", k)
            message(temp)

            newList <- list("dropped.indices" = ind, "geneID" = geneID, "method" = method,
                            "k" = k, "Y" = Y, "X" = X, "cluster" = clus,
                            "clusterLabels" = as.character(clusterLabels),
                            "conductance" = conduct, "cv" = cv)

        } else if(method == "secondOrderGap"){

            # second order method
            ksecondOrderGap <- k$secondOrderGap
            Ys <- Y[, seq_len(min(2*ksecondOrderGap, ncol(Y)))]
            Xs <- X[seq_len(min(2*ksecondOrderGap, length(X)))]
            Ys <- divideNorm(Ys, rowWise = TRUE)

            cluss <- kmeans(Ys, ksecondOrderGap, iter.max = maxIter, nstart = numStart)
            cons <- conductance(adja = adjaMat, clusLab = cluss$cluster)

            method <- "secondOrderGap"
            k <- ksecondOrderGap
            Y <- Ys
            X <- Xs
            clus <- cluss
            clusterLabels <- clus$cluster
            conduct <- cons
            cv <- "userMethod"

            temp <- paste0("\n method ", method, " is selected using user and k is ", k)
            message(temp)

            newList <- list("dropped.indices" = ind, "geneID" = geneID, "method" = method,
                            "k" = k, "Y" = Y, "X" = X, "cluster" = clus,
                            "clusterLabels" = as.character(clusterLabels),
                            "conductance" = conduct, "cv" = cv)

        } else if(method == "additiveGap"){

            # gap method
            kadditiveGap <- k$additiveGap
            Yg <- Y[, seq_len(min(2*kadditiveGap, ncol(Y)))]
            Xg <- X[seq_len(min(2*kadditiveGap, length(X)))]
            Yg <- divideNorm(Yg, rowWise = TRUE)

            clusg <- kmeans(Yg, kadditiveGap, iter.max = maxIter, nstart = numStart)
            cong <- conductance(adja = adjaMat, clusLab = clusg$cluster)

            method <- "additiveGap"
            k <- "additiveGap"
            Y <- Yg
            X <- Xg
            clus <- clusg
            clusterLabels <- clus$cluster
            conduct <- cong
            cv <- "userMethod"

            temp <- paste0("\n method ", method, " is selected using user and k is ", k)
            message(temp)

            newList <- list("dropped.indices" = ind, "geneID" = geneID, "method" = method,
                            "k" = k, "Y" = Y, "X" = X, "cluster" = clus,
                            "clusterLabels" = as.character(clusterLabels),
                            "conductance" = conduct, "cv" = cv)

        } else  {
            stop(" method for number of clusters can be either relativeGap,
                secondOrderGap, or additiveGap")}

    }else if(!is.null(kopt)){

        # optimal by user
        Yopt <- Y[, seq_len(min(2*kopt, ncol(Y)))]
        Xopt <- X[seq_len(min(2*kopt, length(X)))]
        Yopt <- divideNorm(Yopt, rowWise = TRUE)

        clusopt <- kmeans(Yopt, kopt, iter.max = maxIter, nstart = numStart)
        conopt <- conductance(adja = adjaMat, clusLab = clusopt$cluster)

        method <- "userkopt"
        k <- kopt
        Y <- Yopt
        X <- Xopt
        clus <- clusopt
        clusterLabels <- clus$cluster
        conduct <- conopt
        cv <- "userkopt"

        temp <- paste0("\n selected k is ", k, " by user")
        message(temp)

        newList <- list("dropped.indices" = ind, "geneID" = geneID, "method" = method,
                        "k" = k, "Y" = Y, "X" = X, "cluster" = clus,
                        "clusterLabels" = as.character(clusterLabels),
                        "conductance" = conduct, "cv" = cv)
    }



    ################### silhouette
    if(sil){

        message("calculating the Silhouette index \n it may take time...")

        sil <- silhouette(dis.y = dist(Y), clus.labels = clusterLabels)
        sil <- sil %>%  arrange(clusterLabel ,-silIndex)
        sil$clusterLabel <- as.factor(sil$clusterLabel)
        genes <- data.frame(geneID = geneID, geneIndices = seq(1, length(geneID)))
        sil <- inner_join(sil, genes, by = "geneIndices")
        sil <-  sil[ , !(names(sil) %in% "geneIndices")]
        newList <- c(newList, setNames(list(sil), "silhouette"))

    }


    if(saveOrig == TRUE){

        tlist <- list("Yorig" = Yorig, "Xorig" = Xorig, "n_egvec" = ncol(Yorig))
        newList <- c(newList, setNames(list(tlist), "original"))}

    message("network clustering is done...\n")

    return(newList)

}



