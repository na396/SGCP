
################################################################################# SGCP_ezPlot

SGCP_ezPLOT <- function(sgcp, expreData, keep = FALSE,
                        pdf.file = TRUE, pdfname = "ezSGCP.pdf",
                        excel.file = TRUE,
                        xlsxname = "ezSGCP.xlsx", w = 6, h = 6, sr = 2, sc = 2,
                        ftype = "png",  uni = "in",
                        expressionPCA = TRUE, pointSize1 = .5,
                        exprePCATitle0 = "Expression Data PCA Without Labels",
                        exprePCATitle1 = "Expression Data PCA With Initial Labels",
                        exprePCATitle2 = "Expression Data PCA With Final Labels",
                        transformedPCA = TRUE, pointSize2 = 0.5,
                        transformedTitle0 = "Transformed Data PCA Without Labels",
                        transformedTitle1 = "Transformed Data PCA Initial Labels",
                        transformedTitle2 = "Transformed Data PCA Final Labels",
                        conduct = TRUE,
                        conductanceTitle = "Cluster Conductance Index",
                        conductx = "clusterLabel", conducty = "conductance index",
                        clus_num = TRUE,
                        silhouette_index = FALSE,
                        silTitle = "Gene Silhouette Index",
                        silx = "genes", sily = "silhouette index",
                        jitt1 = TRUE,
                        jittTitle1 = "Initial GO p-values", jps1 = 3,
                        jittx1 = "cluster", jitty1 = "-log10 p-value",
                        jitt2 = TRUE,
                        jittTitle2 = "Final GO p-values", jps2 = 3,
                        jittx2 = "module", jitty2 = "-log10 p-value",
                        density1 = TRUE,
                        densTitle1 = "Initial GO p-values Density",
                        densx1 = "cluster", densy1 = "-log10 p-value",
                        density2 = TRUE,
                        densTitle2 = "Final GO p-values Density",
                        densx2 = "module", densy2 = "-log10 p-value",
                        mean1 = TRUE,
                        meanTitle1 = "Cluster Performance",
                        meanx1 = "cluster", meany1 = "mean -log10 p-value",
                        mean2 = TRUE,
                        meanTitle2 = "Module Performance",
                        meanx2 = "module", meany2 = "mean -log10 p-value",
                        pie1 = TRUE, pieTitle1 = "Initial GO Analysis",
                        piex1 = "cluster", piey1 = "count", posx1 = 1.8,
                        pie2 = TRUE, pieTitle2 = "Final GO Analysis",
                        piex2 = "module", piey2 = "count", posx2 = 1.8

){


    pl.final <- list()
    wb <- createWorkbook()
    pdf.out <- vector(15, mode = "list")

    ind <- 1
    counter <- 1
    ############################################# expression PCA plot
    if(expressionPCA){
        message("plotting PCA of expression data...")
        ############################################################
        expPCA.lab0 <- SGCP_plot_pca(m = expreData, clusLabs = NULL,
                                    tit = exprePCATitle0, ps = pointSize1)

        if(excel.file){
            wb <- SGCP_plot2excel(expPCA.lab0, wb, "expressionPCAwithNoLabel",
                                ind , w, h, ftype, uni, sr, sc)
            ind <- ind + 1}

        if(pdf.file){
            expPCA.lab0
            pdf.out[[counter]] <- recordPlot()
            counter <- counter + 1}

        if(keep){
            pl.final <- c(pl.final, setNames(list(expPCA.lab0) ,
                                            "expressionPCAwithNoLabel"))
        }else{ rm(expPCA.lab0) }


        ############################################################
        if(length(sgcp$clustering$dropped.indices)> 0){
            expreData <- expreData[-sgcp$clustering$dropped.indices]
        }
        expPCA.lab1 <- SGCP_plot_pca(m = expreData,
                                    clusLabs = sgcp$clusterLabels$initialClusters,
                                    tit = exprePCATitle1, ps = pointSize1)

        if(excel.file){
            wb <- SGCP_plot2excel(expPCA.lab1, wb, "expressionPCAwithInitialLabel",
                                ind , w, h, ftype, uni, sr, sc)
            ind <- ind +1}

        if(pdf.file){
            expPCA.lab1
            pdf.out[[counter]] <- recordPlot()
            counter <- counter + 1}

        if(keep){
            pl.final <- c(pl.final, setNames(list(expPCA.lab1) ,
                                            "expressionPCAwithInitialLabel"))
        }else{ rm(expPCA.lab1) }

        ############################################################
        if(sgcp$semilabel){
            expPCA.lab2<- SGCP_plot_pca(m = expreData,
                                        clusLabs = sgcp$clusterLabels$finalClusters,
                                        tit = exprePCATitle2, ps = pointSize1)
            if(excel.file){
                wb <- SGCP_plot2excel(expPCA.lab2, wb, "expressionPCAwithFinalLabel",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind +1}

            if(pdf.file){
                expPCA.lab2
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(expPCA.lab2) ,
                                                "expressionPCAwithFinalLabel")) }
        }else{ rm(expPCA.lab2)}

    } # end of expressionPCA

    ######################################### transformed PCA plot
    if(transformedPCA){
        message("plotting PCA of transformed data...")
        ############################################################
        transPCA.lab0 <- SGCP_plot_pca(m = sgcp$clustering$Y,
                                    clusLabs = NULL,
                                    tit = transformedTitle0, ps = pointSize2)

        if(excel.file){
            wb <- SGCP_plot2excel(transPCA.lab0, wb, "transformedPCAwithNoLabel",
                                ind , w, h, ftype, uni, sr, sc)
            ind <- ind +1}

        if(pdf.file){
            transPCA.lab0
            pdf.out[[counter]] <- recordPlot()
            counter <- counter + 1}

        if(keep){
            pl.final <- c(pl.final, setNames(list(transPCA.lab0) ,
                                            "transformedPCAwithNoLabel"))
        }else{ rm(transPCA.lab0)}

        ############################################################
        transPCA.lab1 <- SGCP_plot_pca(m = sgcp$clustering$Y,
                                    clusLabs = sgcp$clusterLabels$initialClusters,
                                    tit = exprePCATitle1, ps = pointSize2)

        if(excel.file){
            wb <- SGCP_plot2excel(transPCA.lab1, wb, "transformedPCAwithInitialLabel",
                                ind , w, h, ftype, uni, sr, sc)
            ind <- ind +1}

        if(pdf.file){
            transPCA.lab1
            pdf.out[[counter]] <- recordPlot()
            counter <- counter + 1}

        if(keep){
            pl.final <- c(pl.final, setNames(list(transPCA.lab0) ,
                                            "transformedPCAwithInitialLabel"))
        }else{ rm(transPCA.lab1) }

        ############################################################
        if(sgcp$semilabel){
            transPCA.lab2 <- SGCP_plot_pca(m = sgcp$clustering$Y,
                                        clusLabs = sgcp$clusterLabels$finalClusters,
                                        tit = exprePCATitle2, ps = pointSize2)

            if(excel.file){
                wb <- SGCP_plot2excel(transPCA.lab2, wb, "transformedPCAwithFinalLabel",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind +1}

            if(pdf.file){
                transPCA.lab2
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}


            if(keep){
                pl.final <- c(pl.final, setNames(list(transPCA.lab0) ,
                                                "transformedPCAwithFinalLabel"))
            }else{ rm(transPCA.lab2) }

        } # end of semilabel

    } # end of transformedPCA

    ################################################ cluster number
    if(clus_num){
        if(sgcp$clustering$cv == "cvGO" || sgcp$clustering$cv == "cvConductance"){
            message("plotting relativeGap, secondOrderGap, additiveGap methods ...")

            rg <- sgcp$clustering$clusterNumberPlot$relativeGap
            sg <- sgcp$clustering$clusterNumberPlot$secondOrderGap
            ag <- sgcp$clustering$clusterNumberPlot$additiveGap

            if(excel.file){
                wb <- SGCP_plot2excel(ag, wb, "clusterNuumberAdditiveGap",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind +1}

            if(pdf.file){
                ag
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}


            if(excel.file){
                wb <- SGCP_plot2excel(rg, wb, "clusterNumberRelativeGap",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind +1}

            if(pdf.file){
                rg
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}


            if(excel.file){
                wb <- SGCP_plot2excel(sg, wb, "clusterNumberSecondOrderGap",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind +1 }


            if(pdf.file){
                sg
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}


            if(keep){
                pl.final <- c(pl.final, setNames(list(rg) ,"clusterNuumberRelativeGap"),
                            setNames(list(sg) ,"clusterNumberSecondOrderGap"),
                            setNames(list(ag) ,"clusterNumberAdditiveGap"))
            }else{ rm(rg, sg, ag)}

        }else {
            message("methods (relativeGap, secondOrderGap, additiveGap) are needed for number of clusters" )
        }
    }

    ############################################# cluster Conductance
    if(conduct == TRUE){

        if(sgcp$clustering$cv == "cvGO" || sgcp$clustering$cv == "cvConductance" ){
            message("plotting cluster conductance index...")
            conduct_plt <- SGCP_plot_conductance(conduct = sgcp$clustering$conductance,
                                                tit = conductanceTitle,
                                                xname = conductx, yname = conducty)

            if(excel.file){
                wb <- SGCP_plot2excel(conduct_plt, wb, "clusterConductance",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind +1}

            if(pdf.file){
                conduct_plt
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}


            if(keep){
                pl.final <- c(pl.final, setNames(list(conduct_plt) ,
                                                "clusterConductance"))
            }else{ rm(conduct_plt)}

        }else{ message("cannot plot conductance \n make sure k is NULL in ezSGCP")}
    }

    ######################################### silhouette index
    if(silhouette_index){

        if(!is.null(sgcp$clustering[["silhouette"]])){
            message("plotting cluster silhouette index...")
            sil_plt <- SGCP_plot_silhouette(sgcp$clustering$silhouette,
                                        tit = silTitle, xname = silx, yname = sily)

            if(excel.file){
                wb <- SGCP_plot2excel(sil_plt, wb, "clusterSilhouetteIndex",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                sil_plt
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(sil_plt) ,
                                                "clusterSilhouetteIndex"))
            }else{ rm(sil_plt) }

        }else {
            message("cannot plot silhouette \n, make sure sil is TRUE in ezSGCP")
        }
    } # end of if silhouette_index

    ######################################### pvalue jitter1
    if(jitt1){
        if(sgcp$semilabel == TRUE){
            message("plotting initial GO p-values jitters...")
            jitt_plt1 <- SGCP_plot_jitter(df = sgcp$initial.GO$GOresults, tit = jittTitle1,
                                        xname = jittx1, yname = jitty1, ps = jps1)


            if(excel.file){
                wb <- SGCP_plot2excel(jitt_plt1, wb, "InitialGOpvalueJitter",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                jitt_plt1
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(jitt_plt1) ,
                                                "InitialGOpvalueJitter"))

            }else{rm(jitt_plt1)}
        }else{
            message("cannot plot jitter \n, make sure semilabel is TRUE in ezSGCP")
            }
    }

    ####################################### pvalue density1
    if(density1){
        if(sgcp$semilabel == TRUE){
            message("plotting initial GO p-values density...")
            den_plt1 <- SGCP_plot_density(df = sgcp$initial.GO$GOresults,
                                        tit = densTitle1,
                                        xname = densx1, yname = densy1)


            if(excel.file){
                wb <- SGCP_plot2excel(den_plt1, wb, "InitialGOpvalueDensity",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                den_plt1
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(den_plt1) ,
                                                "InitialGOpvalueDensity"))
            }else{ rm(den_plt1)}

        }else{
            message("cannot plot density \n, make sure semilabel is TRUE in ezSGC") }

    }

    ######################################## pvalue mean1
    if(mean1){
        if(sgcp$semilabel == TRUE){
            message("plotting cluster performance...")
            mean_plt1 <- SGCP_plot_bar(df = sgcp$initial.GO$GOresults, tit = meanTitle1,
                                    xname = meanx1, yname = meany1)

            if(excel.file){
                wb <- SGCP_plot2excel(mean_plt1, wb, "InitialGOpvalueMean",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                mean_plt1
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(mean_plt1) ,
                                                "InitialGOpvalueMean"))
            }else{ rm(mean_plt1)}

        }else{
            message("cannot plot mean \n, make sure semilabel is TRUE in ezSGCP") }

    }

    #############################################pvalue pei1
    if(pie1){
        if(sgcp$semilabel == TRUE){
            message("plotting final GO analysis...")
            pie_plt1 <- SGCP_plot_pie(df = sgcp$initial.GO$GOresults, tit = pieTitle1,
                                    xname = piex1, yname = piey1, posx = posx1)

            if(excel.file){
                wb <- SGCP_plot2excel(pie_plt1, wb, "InitialGOAnalysis",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                pie_plt1
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(pie_plt1) ,
                                                "InitialGOAnalysis"))
            }else{ rm(pie_plt1)}

        }else{
            message("cannot plot pie \n, make sure semilabel is TRUE in ezSGCP") }

    }

    ########################################## pvalue jitter2
    if(jitt2){
        if(sgcp$semilabel == TRUE){
            message("plotting final GO p-values jitters...")
            jitt_plt2 <- SGCP_plot_jitter(sgcp$final.GO$GOresults, tit = jittTitle2,
                                        xname = jittx2, yname = jitty2, ps = jps2)


            if(excel.file){
                wb <- SGCP_plot2excel(jitt_plt2, wb, "FinalGOpvalueJitter",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                jitt_plt2
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}


            if(keep){
                pl.final <- c(pl.final, setNames(list(jitt_plt2) ,
                                                "FinalGOpvalueJitter"))

            }else{rm(jitt_plt2)}

        }else{
        message("cannot plot jitter \n, make sure semilabel is TRUE in ezSGCP")
        }
    }

    ######################################### pvalue density2
    if(density2){
        if(sgcp$semilabel == TRUE){
            message("plotting final GO p-values density...")
            den_plt2 <- SGCP_plot_density(df = sgcp$final.GO$GOresults,
                                        tit = densTitle2,
                                        xname = densx2, yname = densy2)

            if(excel.file){
                wb <- SGCP_plot2excel(den_plt2, wb, "FinalGOpvalueDensity",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                den_plt2
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(den_plt2) ,
                                                "FinalGOpvalueDensity"))
            }else{rm(den_plt2)}

        }else{
            message("cannot plot density \n, make sure semilabel is TRUE in ezSGCP") }

    }

    ##################################### pvalue mean2
    if(mean2){
        if(sgcp$semilabel == TRUE){
            message("plotting module performance...")
            mean_plt2 <- SGCP_plot_bar(df = sgcp$final.GO$GOresults, tit = meanTitle2,
                                    xname = meanx2, yname = meany2)


            if(excel.file){
                wb <- SGCP_plot2excel(mean_plt2, wb, "FinalGOpvalueMean",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                mean_plt2
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(mean_plt2) ,
                                                "FinalGOpvalueMean"))
            }else{ rm(mean_plt2)}

        }else{
            message("cannot plot mean \n, make sure semilabel is TRUE in ezSGCP") }

    }

    ################################################## pvalue pei2
    if(pie2){
        if(sgcp$semilabel == TRUE){
            message("plotting final GO analysis...")
            pie_plt2 <- SGCP_plot_pie(df = sgcp$final.GO$GOresults, tit = pieTitle2,
                                    xname = piex2, yname = piey2, posx = posx2)


            if(excel.file){
                wb <- SGCP_plot2excel(pie_plt2, wb, "FinalGOAnalysis",
                                    ind , w, h, ftype, uni, sr, sc)
                ind <- ind + 1}

            if(pdf.file){
                pie_plt2
                pdf.out[[counter]] <- recordPlot()
                counter <- counter + 1}

            if(keep){
                pl.final <- c(pl.final, setNames(list(pie_plt2) ,
                                                "FinalGOAnalysis"))
            }else{rm(pie_plt2)}

        }else{
            message("cannot plot pie \n, make sure semilabel is TRUE in ezSGCP") }

    }

    ###################################
    if(excel.file){
        SGCP_saveResult(wb, xlsxname)}

    if(pdf.file){

        pdf(pdfname, onefile = TRUE)
        for(plt in pdf.out){
            if(!is.null(plt)){
                replayPlot(plt)
            }
        }
        graphics.off()
    }

    if(keep){ return(pl.final) }

} # end of ezPLOT






