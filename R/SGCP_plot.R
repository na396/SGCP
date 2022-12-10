
################################################################################# plot_extension
SGCP_plot_extension <- function(pl, xname, yname, tit,
                                xtextsize = 7, xtextface = "bold",
                                xangle = 360, yangle = 360,
                                xtitlesize = 10, xtitleface = "bold",
                                ytextsize = 15, ytextface = "bold",
                                ytitlesize = 10, ytitleface = "bold",
                                titlesize = 15, titleface = "bold",
                                legendsize = 7, legendface = "bold",
                                legendtitlesize = 7, legendtitleface = "bold",
                                legendpos = "bottom"){

    res <- pl +
        theme_classic()  +

        theme(axis.text.x = element_text(size = xtextsize, angle = xangle,
                                    face = xtextface, lineheight = 0.9)) +
        theme(axis.title.x = element_text(size = xtitlesize,
                                    face = xtitleface,lineheight = 0.9)) +
        theme(axis.text.y = element_text(size = ytextsize, angle = yangle,
                                    face = ytextface, lineheight = 0.9)) +
        theme(axis.title.y = element_text(size = ytitlesize,
                                    face = ytitleface,lineheight = 0.9)) +
        theme(plot.title = element_text(size = titlesize,
                                    face = titleface, lineheight = 0.9,
                                    hjust = 0.5)) +

        theme(legend.text = element_text(size = legendsize, face = legendface)) +
        theme(legend.title = element_text(size = legendtitlesize,
                                        face = legendtitleface))+
        theme(legend.position = legendpos) +
        labs(x = xname, y = yname, title = tit)

    return(res)
}

################################################################################# plot_pca
SGCP_plot_pca <- function(m, clusLabs, tit = "PCA plot", ps = .5){

    # plot the pca of m
    # clusLabs: cluster labels
    # m: dataset for PCA

    if(!is.null(clusLabs)){
        clusLabs <- as.numeric(clusLabs) }

    if(!is.null(clusLabs)){
        mat_pca <- prcomp(m)
        mat_pca <- as.data.frame(mat_pca$x)
        mat_pca <- mat_pca[, c(1,2)]
        mat_pca <- data.frame(mat_pca, label = clusLabs)
        mat_pca$label <- as.factor(mat_pca$label)

        plt_trans <- ggplot(mat_pca,
                            aes(x = mat_pca[,1],y = mat_pca[,2], color = label )) +
            geom_point(size = ps)

    } else{

        mat_pca <- prcomp(m)
        mat_pca <- as.data.frame(mat_pca$x)
        plt_trans <- ggplot(mat_pca,aes(x = mat_pca[,1],y = mat_pca[,2] )) +
            geom_point(size = ps) }


    plt_trans <- SGCP_plot_extension(plt_trans,
                                xname = "component1", yname = "component2",
                                tit = tit,
                                xtextsize = 10, xtextface = "bold", xangle = 360,
                                xtitlesize = 12, xtitleface = "bold", yangle = 360,
                                ytextsize = 10, ytextface = "bold",
                                ytitlesize = 12, ytitleface = "bold",
                                titlesize = 15, titleface = "bold",
                                legendsize = 7, legendface = "bold",
                                legendtitlesize = 7, legendtitleface = "bold",
                                legendpos = "bottom")


    return(plt_trans)
}

################################################################################# plot_silhouette
SGCP_plot_silhouette <- function(df, tit = "Gene Silhouette Index",
                                xname = "genes", yname = "silhouette index"){
    # df: df of gene index with their corresponding silhouette index

    mm <- mean(df$silIndex)
    df$index <- seq(1, nrow(df))
    df$clusterLabel <- as.factor(as.numeric(df$clusterLabel))
    plt_sil <- ggplot(df, aes(x = index, y = silIndex, fill = clusterLabel)) +
        geom_bar(stat = 'identity') +
        geom_hline(yintercept = mm, linetype="dashed")



    plt_sil <- SGCP_plot_extension(plt_sil,
                                xname = xname, yname = yname,
                                tit = tit,
                                xtextsize = 0, xtextface = "bold", xangle = 360,
                                xtitlesize = 12, xtitleface = "bold", yangle = 360,
                                ytextsize = 10, ytextface = "bold",
                                ytitlesize = 12, ytitleface = "bold",
                                titlesize = 15, titleface = "bold",
                                legendsize = 7, legendface = "bold",
                                legendtitlesize = 7, legendtitleface = "bold",
                                legendpos = "right")

    return(plt_sil)


}

################################################################################# plot_heatMap
SGCP_plot_heatMap <- function(m, tit = "Adjacency Heatmap",
                            xname = "genes", yname = "genes"){

    # m: the adjacency matrix, square symmetric matrix

    melted_m <- melt(m)
    n <- min(nrow(melted_m), 1000)

    geo_heatmap <- ggplot(melted_m, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_distiller(palette = 'Spectral') +
        theme_classic() +

        theme(axis.title.x = element_text(size = 10, face = 'bold',
                                        lineheight = 0.9)) +
        theme(axis.title.y = element_text(size = 10, face = 'bold',
                                        lineheight = 0.9)) +
        theme(plot.title = element_text(size = 15,face = 'bold',
                                        lineheight = 0.9, hjust = 0.5)) +
        theme(legend.title = element_text(size = 8, face="bold")) +
        theme(legend.text = element_text(size = 8)) +
        labs(x = xname, y = yname, title = tit)

    rm(melted_m)

    return(geo_heatmap)

}

################################################################################# plot_conductance
SGCP_plot_conductance <- function(conduct, tit = "Clustering Conductance Index",
                                xname = "cluster", yname = "conductance"){

    l.f <- length(conduct$relativeGap.conductance$clusterConductance)
    con.f <- conduct$relativeGap.conductance$clusterConductance
    methFirst <- data.frame(conductVal = con.f, Method = rep("relativeGap",
                                                        length(con.f)),
                            lb = paste0("rg", as.numeric(names(con.f))),
                            lbs = as.numeric(names(con.f)))
    methFirst <- methFirst[order(methFirst$lbs), ]
    methFirst$lb <- factor(methFirst$lb, levels = methFirst$lb)

    l.s <- length(conduct$secondOrderGap.conductance$clusterConductance)
    con.s <- conduct$secondOrderGap.conductance$clusterConductance
    methSecond <- data.frame(conductVal = con.s, Method = rep("secondOrderGap",
                                                        length(con.s)),
                            lb = paste0("sg", as.numeric(names(con.s))),
                            lbs = as.numeric(names(con.s)))
    methSecond <- methSecond[order(methSecond$lbs), ]
    methSecond$lb <- factor(methSecond$lb, levels = methSecond$lb)

    l.g <- length(conduct$additiveGap.conductance$clusterConductance)
    con.g <- conduct$additiveGap.conductance$clusterConductance
    methGap <- data.frame(conductVal = con.g, Method = rep("additiveGap",
                                                        length(con.g)),
                            lb = paste0("ag", as.numeric(names(con.g))),
                            lbs = as.numeric(names(con.g)))
    methGap <- methGap[order(methGap$lbs), ]
    methGap$lb <- factor(methGap$lb, levels = methGap$lb)

    fconduct <- rbind(methGap, methFirst, methSecond)
    fconduct$index <- seq(1, nrow(fconduct))
    fconduct$Method <- factor(fconduct$Method,
                            levels = c("additiveGap", "relativeGap",
                                        "secondOrderGap"))

    maxClus <- max(l.f, l.s, l.g)

    plt_conduct <- ggplot(fconduct, aes(x = lb, y = conductVal, fill = Method)) +
        geom_bar(stat = 'identity', position=position_dodge(), colour="black") +
        scale_fill_discrete(
            breaks = c( "additiveGap","relativeGap", "secondOrderGap"))

    xangle <- 360
    if(maxClus > 4){xangle <- 90}

    plt_conduct <- SGCP_plot_extension(plt_conduct,
                                    xname = xname, yname = yname,
                                    tit = tit,
                                    xtextsize = 10, xtextface = "bold", xangle = xangle,
                                    xtitlesize = 12, xtitleface = "bold", yangle = 360,
                                    ytextsize = 10, ytextface = "bold",
                                    ytitlesize = 12, ytitleface = "bold",
                                    titlesize = 15, titleface = "bold",
                                    legendsize = 7, legendface = "bold",
                                    legendtitlesize = 7, legendtitleface = "bold",
                                    legendpos = "right")


    return(plt_conduct)
}

################################################################################# plot_density
SGCP_plot_density <- function(df, tit = "p-values Density",
                            xname = "module", yname = "-log10 p-value"){

    # df returned by gene ontology
    names(df)[names(df) == "clusterNum"] <- "clusterLabel"
    if(!is.numeric(df$clusterLabel)){
        df$clusterLabel <- as.numeric(df$clusterLabel) }

    df$clusterLabel <- factor(df$clusterLabel, levels = sort(unique(df$clusterLabel)))
    df$Temp <- "GO p-value Densitys"

    plt_density <- ggplot(df, aes(x = -log10(Pvalue), y = clusterLabel,
                                fill = clusterLabel)) +
        geom_density_ridges() +
        coord_flip() +
        facet_wrap(~Temp, scales = "free", ncol = 1) +


        theme_classic()  +

        theme(axis.text.x = element_blank()) +
        theme(axis.title.x = element_blank() ) +

        theme(axis.text.y = element_text(size = 13, face = 'bold', lineheight = 0.9)) +
        theme(axis.title.y = element_text(size = 14, face = 'bold', lineheight =.9)) +

        theme(plot.title =  element_text(size = 16, face = 'bold', lineheight = 0.5,
                                        hjust = 0.5)) +

        theme(legend.text = element_text(face = "bold") ) +

        theme(strip.background = element_blank()) +
        theme(strip.text = element_blank()) +

        theme(legend.position = "right") +
        theme(legend.title = element_text(size = 8, face = "bold", lineheight = 0.5)) +

        labs(y = "", x = "-log10 p-value",  title = tit )

    plt_density <- plt_density + scale_fill_discrete(name = xname)

    return(plt_density)

}

################################################################################# plot_jitter
SGCP_plot_jitter <- function(df, tit = "p-values Distribution",
                            xname = "module", yname = "-log10 p-value", ps = 3){

    df$clusterNum <- as.numeric(df$clusterNum)
    df$clusterNum <- factor(df$clusterNum, levels = sort(unique(df$clusterNum)) )


    names(df)[names(df) == "clusterNum"] <- "clusterLabel"
    plt_jitter <- ggplot(df, aes(x = clusterLabel, y = -log10(Pvalue),
                                colour = clusterLabel)) +

        geom_jitter(size = ps, height = 0)

    plt_jitter <- SGCP_plot_extension(plt_jitter,
                                    xname = xname, yname = yname, tit = tit,
                                    xtextsize = 0, xtextface = "bold",
                                    xangle = 360, yangle = 360,
                                    xtitlesize = 0, xtitleface = "bold",
                                    ytextsize = 15, ytextface = "bold",
                                    ytitlesize = 10, ytitleface = "bold",
                                    titlesize = 15, titleface = "bold",
                                    legendsize = 7, legendface = "bold",
                                    legendtitlesize = 7, legendtitleface = "bold",
                                    legendpos = "right")

    plt_jitter <- plt_jitter + scale_colour_discrete(name = xname)

    return(plt_jitter)

}

################################################################################# plot_mean
SGCP_plot_bar <- function(df, tit = "mean -log10 p-values",
                        xname = "module", yname = "-log10 p-value"){

    #df: returned by GOstats

    if(!is.numeric(df$clusterNum)){df$clusterNum <- as.numeric(df$clusterNum)}
    df$index <- unclass(df$clusterNum)
    df$clusterNum <- factor(df$clusterNum, levels = sort(unique(df$clusterNum)))
    df$logPvalue <- -log10(df$Pvalue)
    names(df)[names(df) == "clusterNum"] <- "clusterLabel"
    clusNums <- sort(unique(df$clusterLabel))

    df.sd <- data.frame(matrix(nrow = nlevels(df$clusterLabel), ncol = 1))
    colnames(df.sd) <- "sd"
    rownames(df.sd) <- levels(df$clusterLabel)

    for(clus in levels(df$clusterLabel)){
        temp <- df[df$clusterLabel == clus, ]
        ff <- MeanCI(temp$logPvalue, sides = "two.sided")
        df.sd[clus, "sd"] <- ff[1] - ff[2]
        rm(temp)}
    df.sd$clusterLabel <- rownames(df.sd)

    temp <- aggregate(logPvalue~clusterLabel, data = df, mean)
    temp <- inner_join(temp, df.sd, by = "clusterLabel")
    temp$clusterLabel <- factor(temp$clusterLabel, levels = clusNums)

    plt_mean <- ggplot(temp, aes(clusterLabel, logPvalue, fill = clusterLabel)) +
        geom_bar(stat = "identity", position = "dodge", colour = "black") +
        geom_errorbar(aes(ymin = logPvalue - sd, ymax= logPvalue + sd),
                    width = 0.5)


    plt_mean <- SGCP_plot_extension(plt_mean,
                                xname = xname, yname = yname,
                                tit = tit,
                                xtextsize = 0, xtextface = "bold", xangle = 360,
                                xtitlesize = 0, xtitleface = "bold",
                                ytextsize = 10, ytextface = "bold",
                                ytitlesize = 12, ytitleface = "bold",
                                titlesize = 15, titleface = "bold",
                                legendsize = 7, legendface = "bold",
                                legendtitlesize = 7, legendtitleface = "bold",
                                legendpos = "right")


    plt_mean <- plt_mean + scale_fill_discrete(name = xname)

    return(plt_mean)

}

################################################################################# piechart
SGCP_plot_pie <- function(df, tit = "GO Analysis",
                        xname = "module", yname = "count", posx = 1.9){

    if(!is.numeric(df$clusterNum)){ df$clusterNum <- as.numeric(df$clusterNum)}
    df <- df[order(df$clusterNum), ]
    clusNums <- as.character(unique(df$clusterNum))

    temp <- df[, c("clusterNum", "GOtype", "Pvalue")]
    temp$logPvalue <- -log10(df$Pvalue)

    ppl <- temp %>%
            group_by(clusterNum, GOtype) %>%
            summarise(max = max(logPvalue), count = n())

    cluslabs <- unique(ppl$clusterNum)
    onto <- c("overBP", "overCC", "overMF", "underBP", "underCC", "underMF")

    for(c in cluslabs){

        t <- ppl[ppl$clusterNum == c, ]
        difEle <- setdiff(onto, t$GOtype)
        if(length(difEle) > 1){
            addEle <- data.frame(clusterNum = c, "GOtype" = difEle,
                                "max" = 0, "count" = 0)
            ppl <- rbind(ppl, addEle) }

    }

    ppl$max <- round(ppl$max, 2)
    ppl$clusterNum <- paste0(rep(xname, nrow(ppl)), ppl$clusterNum)
    ppl$clusterNum <- factor(ppl$clusterNum, levels = paste0(xname, clusNums))


    plt_pie <- ggplot(ppl, aes(x = "", y = count, fill = GOtype)) +
        geom_bar(stat = "identity", position = position_fill()) +
        geom_text(aes(x = posx, label = max),
                position = position_fill(vjust = 0.5),
                size = 2, fontface = "bold") +
        coord_polar("y", start=0) +
        facet_wrap(~clusterNum) +

        theme_classic()  +
        theme(strip.background = element_blank(),
            strip.text.x = element_text(size = 15,
                                        face = 'bold',lineheight = 0.9)) +

        theme(legend.text = element_text(size = 10, face="bold")) +

        theme(axis.text.x = element_blank()) +
        theme(axis.text.y = element_text(size = 15,
                                        face = 'bold',lineheight = 0.9)) +

        theme(axis.title.y = element_text(size = 10,
                                        face = 'bold',lineheight = 0.9)) +
        theme(axis.title.x = element_text(size = 10,
                                        face = 'bold',lineheight = 0.9)) +

        theme(plot.title = element_text(size = 20,
                                        face = 'bold',lineheight = 0.9,
                                        hjust = 0.5)) +
        theme(legend.title = element_text(size = 10, face = "bold"))+
        scale_fill_discrete(name = "GeneOntology") +
        theme(legend.position = "bottom") +
        labs(x = "", y = "", title = tit)

    return(plt_pie)
}

################################################################################# plot2excel
SGCP_plot2excel <- function(pl, wb, shname,
                            ind, wid = 6, heigh = 6, ftype = "png",
                            uni = "in", sr = 2, sc = 2){

    # saves the plots in excel file

    # pl: the plots
    # the saving file
    # index of the sheet


    addWorksheet(wb, shname, gridLines = FALSE)
    show(pl)

    insertPlot(wb, ind, width = wid, height = heigh, startRow = sr, startCol = sc,
            fileType = ftype, units = uni)

    return(wb)

}

#################################################################################
SGCP_saveResult <- function(wb, xlsxName){

    saveWorkbook(wb, xlsxName, overwrite = TRUE)
}





