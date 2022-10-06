################################################################################
################################################################################
################################################################################
############################# SGC R package ####################################
############################## Gene Ontology ###################################
################################################################################
################################################################################
################################################################################

################################################################################# Package Dependencies


################################################################################# geneLabeling
GeneOfGOTerm = function(hg, hg_summary, geneID, condition){

  # Identifies the significant geneIDs
  # Input:
  #     hg: an object returned by GOstats
  #     df_summary:  a summary dataframe over hg object
  #     geneID: geneIDs in a given cluster
  #
  # Value:
  #       newList: a list containing labeledGenes, and GOTermgenes
  #       labeledGenes, all significant geneEntrezIDs
  #       GOTermgenes, all the geneIDs assigned to sinificant GO terms
  #                   that fall in the cluster
  #

  GO_Genes = hg@goDag@nodeData@data
  labeledGenes = -1
  GOTermGenes = list()

  for(ind in 1:nrow(hg_summary)){
    #
    go = hg_summary[ind, 1]

    if(condition == TRUE) {
      temp = GO_Genes[[go]]$condGeneIds
    }else{ temp = GO_Genes[[go]]$geneIds}

    GOTermGenes[[go]] = intersect(temp, geneID)
    labeledGenes = c(labeledGenes, temp)
    }

  labeledGenes = labeledGenes[-1]
  newList = list("labeledGenes" = labeledGenes, "GOTermGenes" = GOTermGenes)

  return(newList)
}


################################################################################# GOenrichment
GOenrichment = function(universe_geneIDs, cluster_geneIDs, annotation_db, lab,
                        direction, ontology, hgCutoff = NULL, condition){

  # perform GO enrichment using GOstats package
  #
  # Arguments:
  #         universe_geneIDs: a vector containing the geneIDs of
  #                           the entire epxression data
  #         cluster_geneIDs: a vector containing the geneIDs of
  #                           the nodes within a cluster
  #         annotation_db: annotation db by user
  #         lab: an integer, the cluster number or label
  #         direction: a vector including c("over", "under")
  #         onotlogy: onology of interest including c("BP", "CC", "MF")
  #         hgCutoff: a number between 0 to 1, lower bound on the significant pvalues
  #         condition: TRUE or FALSE
  #                 if TRUE, conditional test will be performed
  # Values:
  #       newList: a list containing labeledGenes, GOTermGenes, GO_summary
  #       labeledGenes, a vector containing the geneIDs of significant genes at
  #                     hgCutoff level
  #       GOTermGenes, a vector contating the geneIDs that appeared in a GOTerm
  #       GO_summary, a dataframe containing the summary of the GO terms at the
  #                   hgCutoff significance level


  if(is.null(condition)){
    warning(" undefined condition!", call. = FALSE)
    message("setting cond TRUE")
    condition = TRUE}

  if(typeof(universe_geneIDs) != typeof(cluster_geneIDs)) {
    warning(" unverse_geneIDs and cluster_geneIDs must be the same type!",
            call. = FALSE)}

  if(length(cluster_geneIDs) >= length(universe_geneIDs)){
    stop(" cluster_geneIDs must be smaller than universe_geneIDs")}

  if(is.null(annotation_db)){ stop("annotation_db is missing!")}

  if(all(direction %!in% c("under", "over"))){
    warning("direction must be in c(under or over)", call. = FALSE)
    message("making direction equal to under or over")
    direction = c("under", "over")}

  if(length(direction) > 2 ){
        warning("direction must be in c(under or over)", call. = FALSE)
        message("making direction equal to under or over")
        direction = c("under", "over")}


  if(all(ontology %!in% c("BP", "CC", "MF"))){
    warning(" ontology must be in BP CC MF", call. = FALSE)
    message("making ontology to BP, CC, MF")
    ontology = c("BP", "CC", "MF") }

    if(length(ontology) > 3){
        warning(" ontology must be in BP CC MF", call. = FALSE)
        message("making ontology to BP, CC, MF")
        ontology = c("BP", "CC", "MF") }

  if(!is.null(hgCutoff) && (hgCutoff>=1 || hgCutoff <= 0)){
    warning(" not correct hgCutoff value", call. = FALSE)
    message("making hgCuroff to default")
    hgCutoff = NULL}

  if(condition != TRUE && condition != FALSE){
    warning("cond must be boolean!", call. = FALSE)
    message("making condition to TRUE")
    condition = TRUE}

  #Creating the summary data frame for all possible ontology outcomes
  GO_summary = data.frame(matrix(nrow = 0, ncol = 9))
  summary_column_names = c("clusterNum", "GOtype", "GOID" ,"Pvalue", "OddsRatio",
                           "ExpCount", "Count", "Size", "Term")
  colnames(GO_summary) = summary_column_names

  labeledGenes = -1
  GOTermGenes = list()

  for(direct in direction){
    #
    for(onto in ontology){
      #

      if(!is.null(hgCutoff)){
        params = new("GOHyperGParams",
                     geneIds = cluster_geneIDs,
                     universeGeneIds = universe_geneIDs,
                     annotation = annotation_db,
                     ontology = onto,
                     pvalueCutoff = hgCutoff,
                     conditional = condition,
                     testDirection = direct)
      }else{
        params = new("GOHyperGParams",
                     geneIds = cluster_geneIDs,
                     universeGeneIds = universe_geneIDs,
                     annotation = annotation_db,
                     ontology = onto,
                     #pvalueCutoff = hgCutoff,
                     conditional = condition,
                     testDirection = direct) }

      caption = paste0("calling GOstats for ", paste0(direct, onto), "...")
      message(caption)
      hg = hyperGTest(params)
      df_hg = summary(hg)

      if(nrow(df_hg) != 0){

        nl = GeneOfGOTerm(hg, df_hg, cluster_geneIDs, condition = condition)
        message("identifying  genes for each GOTerm...")
        labeledGenes = c(labeledGenes, nl[["labeledGenes"]])
        GOTermGenes = c(GOTermGenes, nl[["GOTermGenes"]])

        df_hg$clusterNum = lab
        df_hg$GOtype = paste0(direct, onto)
        df_hg = df_hg[,c(8,9,1,2,3,4,5,6,7)]
        col_temp = paste("GO", onto, "ID", sep = "")
        colnames(df_hg)[colnames(df_hg) == col_temp] = "GOID"

        GO_summary = rbind.fill(GO_summary, df_hg)
        rm(hg, df_hg)
        } # end of if

    } # for loop onto

  } # for loop direction


  #sort the gene ontology result of all combination based on the p-values
  GO_summary = GO_summary[order(GO_summary$Pvalue),]

  if(length(labeledGenes) > 1){
    labeledGenes = labeledGenes[-1]
    labeledGenes = unique(labeledGenes)}

  newList = list("labeledGenes" = labeledGenes,
                 "GOTermGenes" = GOTermGenes,
                 "GO_summary" = GO_summary)
  return(newList)

}


################################################################################# geneOntology
geneOntology = function(geneUniv, clusLab, annotation_db,
                        direction = c("over", "under"),
                        ontology = c("BP", "CC", "MF"), hgCutoff = NULL,
                        cond = TRUE){

  # perform gene ontology using GOstats package
  #
  # Arguments:
  #           geneUniv: a vctor of all the geneIDs in the expression dataset
  #           clusLab: a vector of cluster label for each geneID
  #           annotation_db: annotation db
  #           direction: a vector including c("over", "under")
  #           onotlogy: ontology of interest including c("BP", "CC", "MF")
  #           hgCutoff: a number between 0 to 1,
  #           cond: TRUE or FALSE
  #                   if TRUE, conditional test will be performed
  #


  if(all(direction %!in% c("under", "over"))){
    warning("direction must be in c(under or over) \n making to default",
            call. = FALSE)
    direction = c("over", "under")}

  if(length(direction) > 2){
    warning("direction must be in c(under or over) \n making to default",
            call. = FALSE)
    direction = c("over", "under")}

  if(all(ontology %!in% c("BP", "CC", "MF")) ){
    warning(" ontology must be in BP CC MF \n making to default", call. = FALSE)
    ontology = c("BP", "CC", "MF")}

  if(length(ontology) > 3){
    warning(" ontology must be in BP CC MF \n making to default", call. = FALSE)
    ontology = c("BP", "CC", "MF")}

  if(!is.null(hgCutoff) && (hgCutoff >=1 || hgCutoff <= 0)){
    warning(" not correct hgCutoff value \n making to default", call. = FALSE) }

  if(cond != TRUE && cond != FALSE){
    warning(" cond must be boolean! \n making to deafult", call. = FALSE)
    cond = TRUE}

  if(length(geneUniv) != length(clusLab)){
    stop(" length of the geneUni must be equal to length of clusLab")}

  if(is.null(annotation_db)){ stop("annotation_db must be defined by user!") }
  # if(!is.character(annotation_db)){
  #   stop("type of annotation_db must be character") }

  if(is.factor(clusLab)){
    clusLab = data.frame(clusterLabel = clusLab)
    clusLab = transform(clusLab,
                clusterLabel = as.numeric(levels(clusterLabel))[clusterLabel])
    clusLab = clusLab[["clusterLabel"]] }

  if(typeof(geneUniv) != typeof(clusLab)){
    warning(" type of geneUniv and clusLab not equal", call. = FALSE)}

  geneClus = rbind(geneUniv, clusLab)
  geneClus = as.data.frame(t(geneClus))
  colnames(geneClus) = c("geneID", "clusterLabel")

  # initialization
  GO_res = data.frame(matrix(nrow = 0, ncol = 9))
  colnames(GO_res) = c("clusterNum", "GOtype", "GOID" ,"Pvalue",
                       "OddsRatio", "ExpCount", "Count", "Size", "Term")

  # creating a list for all GO genes
  FinalGOTermGenes = list()

  for(lab in unique(geneClus$clusterLabel)){

    message(paste0("GOenrichment for cluster ", lab))
    cluster = geneClus$geneID[which(geneClus$clusterLabel == lab)]

    nl = GOenrichment(universe_geneIDs = geneUniv,
                      cluster_geneIDs = cluster,
                      annotation_db = annotation_db,
                      lab = lab,
                      direction = direction,
                      ontology = ontology,
                      hgCutoff = hgCutoff,
                      condition = cond)

    GO_res = rbind.fill(GO_res, nl[["GO_summary"]])

    caption = paste0("Cluster", lab, "_GOTermGenes")
    temp = setNames(list(nl[["GOTermGenes"]]), caption)
    FinalGOTermGenes = c(FinalGOTermGenes, temp)

    rm(caption, cluster, temp, nl)

  }

  newList = list("GOresults" = GO_res, "FinalGOTermGenes" = FinalGOTermGenes)

  message("gene ontology done!..\n")
  return(newList)
}


