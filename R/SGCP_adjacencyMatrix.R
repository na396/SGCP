################################################################################# Package Dependencies

################################################################################ %!in% operation
'%!in%' <- function(x,y)!('%in%'(x,y))


################################################################################# df2mat
df2mat <- function(df){
    # convert df into matrix

    if(is.data.frame(df)){

        colnames(df) <- NULL
        rownames(df) <- NULL
        df <- as.matrix(df)}
    return(df)
}


################################################################################# checkSym
checkSym <- function(mat, stp){
    # checks if mat is symmetric
    # checks if values in mat are between (0, 1)

    caption_sym <- paste0(" output of ", stp, " , is not symmetric")
    if(!isSymmetric(mat))
        stop(caption)

    caption_01 <- paste0(" output of ", stp, "are not in (0,1)")
    if(length(table(between(mat, 0, 1)["FALSE"])) != 0)
        stop(caption_01)

    rm(caption_sym, caption_01)

}


################################################################################# checkNumeric
checkNumeric <- function(x, stp){
    # checks if all the values in x are numeric

    #res <- which(vapply(x, class, numeric(1)) != "numeric")
    caption <- paste0(" at least one non-numeric value at ", stp)
    if(!is.numeric(x))
        stop(caption)

}


################################################################################# sigmoid
sigmoid <- function(x, p=1, q=0){
    # sigmoid function

    res <- 1/(1 + exp(-p *(x - q) ))
    return(res)
}


################################################################################# calibration
calibration <- function(v){
    # perform calibraion on vector v

    average <- mean(v)
    variance <- var(v)
    res <- sigmoid(v, p = (1/variance), q = average)
    return(res)
}


################################################################################# norm_vec
normalization <- function(x){
    # divide the vector x by its norm2

    res <- x/sqrt(sum(x^2))
    return(res)
}


################################################################################# GaussianKernel
GaussianKernel <- function(x, sigma){
    # calculates Gaussian kernel for matrix x

    res <- exp(-1 * as.matrix(dist(x)^2)/sigma)
    return(res)
}


################################################################################# DOM: deepOverlapMeasure
DOM <- function(mat){

    # mat: matrix, squared, symmetric, values are in (0,1)
    #
    # value: Deep Overlap Measure (DOM) on mat

    diag(mat) <- 0

    degreeRow <- replicate(dim(mat)[1], rowSums(mat))
    degreeCol <- t(replicate(dim(mat)[1], colSums(mat)))
    degreeMin <- pmin(degreeRow, degreeCol)
    rm(degreeRow, degreeCol)

    degreeRow <- replicate(dim(mat)[1], rowSums(mat)^2)
    degreeCol <- t(replicate(dim(mat)[1], colSums(mat)^2))
    degreeMin2 <- pmin(degreeRow, degreeCol)

    numerator <- mat + (mat %^% 2) + (mat %^% 3)
    denominator <- degreeMin2 + degreeMin + (1 - mat)

    res <- numerator/denominator
    diag(res) <- 1

    rm(degreeCol, degreeRow, degreeMin, degreeMin2)
    return(as.matrix(res))

}

################################################################################# TOM: Topological Overlap Measure
TOM <- function(mat){
    # mat: matrix, squared, symmetric, values are in (0,1)
    # Niloofar Aghaieabiane
    # October 2021
    # value: Topology Overlap Matrix (TOM) on mat

    diag(mat) <- 0

    degreeRow <- replicate(dim(mat)[1], rowSums(mat))
    degreeCol <- t(replicate(dim(mat)[1], colSums(mat)))
    degreeMin <- pmin(degreeRow, degreeCol)

    numerator <- (mat %^% 2) + mat
    denominator <- degreeMin + (1 - mat)

    res <- numerator/denominator
    diag(res) <- 1

    rm(degreeCol, degreeRow, degreeMin)
    return(as.matrix(res))
}

################################################################################# adjacencyMatrix
adjacencyMatrix <- function(expData, calibration = FALSE, norm = TRUE,
                            tom = TRUE, saveAdja = FALSE,
                            adjaNameFile = "adjacency.RData",
                            hm = "adjaHeatMap.png"){


    if(!is.data.frame(expData) && !is.matrix(expData)){
        stop(" the expressoion input must be eighter a data frame or a matrix")}

    if(ncol(expData) >= nrow(expData)){
        warning("number of genes is smaller than the samples,
            are you sure that rows denote genes??", call. = FALSE)}
    checkNumeric(expData, "expression input")


    if(is.data.frame(expData)){
        expData <- df2mat(expData)
        colnames(expData) <- NULL
        rownames(expData) <- NULL }

    expData <- t(expData)

    if(calibration == TRUE){
        message("calibration...")
        expData <- as.data.frame(lapply(as.data.frame(expData), calibration))
        expData <- df2mat(expData)}


    if(norm == TRUE){
        message("normalization...")
        expData <- as.data.frame(lapply(as.data.frame(expData), normalization))
        expData <- df2mat(expData)}

    message("Gaussian kernel...")
    message("it may take time...")
    totalDis <- dist(t(as.matrix(expData)), method = "euclidean")
    adja <- GaussianKernel( x = t(as.matrix(expData)), sigma = var(totalDis))
    rm(totalDis, expData)
    checkSym(adja, "Gussian kernel")

    adja <- df2mat(adja)


    if(tom == TRUE){
        message("TOM...\n it may take time...")
        adja <- TOM(adja)
        checkSym(adja, "TOM")

    }

    diag(adja) <- 0

    if(saveAdja){
        if(!is.character(adjaNameFile)){
            warning("adjaNamFile is not string", call. = FALSE)
            temp <- paste0("using following name \n ",  " adjacency.RData")
            message(temp)
            adjaNameFile <- "adjacency.RData"}

        save(adja, file = adjaNameFile)
        message("adjancency matrix is stored") }


    if(!is.null(hm)){
        hm_plt <- SGCP_plot_heatMap(adja, tit = "Adjacency Heatmap",
                                    xname = "genes", yname = "genes")
        jpeg(hm)
        show(hm_plt)
        dev.off()
        rm(hm_plt)}


    message("network is created, done!...\n")
    return(adja)


}

