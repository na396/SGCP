
################################################################################# Package Dependencies

#################################################################################
semiSupervised <- function(specExp, geneLab, model = "knn", kn = NULL){


    # creating train and test

    if(nrow(geneLab) != nrow(specExp)){
        stop(" number of rows in geneLab and specExp are not equal")}

    if(!is.null(model) & model != "knn" & model != "lr"){
        warning("model must be either NULL, knn, or lr \n setting to knn")
        model <- "knn" }

    if(is.factor(geneLab$label)){
        geneLab <- transform(geneLab, label = as.character(levels(label))[label])}

    message("creating train and test sets based on remarkable and unremarkable genes...")

    rownames(specExp) <- geneLab$geneID
    train <- specExp[rownames(specExp) %in% geneLab$geneID & !is.na(geneLab$label),  ]
    test <- specExp[rownames(specExp) %in% geneLab$geneID & is.na(geneLab$label),  ]

    train <- as.data.frame(train)
    test <- as.data.frame(test)

    caption <- paste0(rep("col", ncol(train)),seq(1, ncol(train)))
    colnames(train) <- caption
    colnames(test) <- caption
    rm(caption)

    train <- data.frame(train, geneID = rownames(train))
    train <- inner_join(train, geneLab, by = "geneID" )
    rownames(train) <- train$geneID
    train <- train[, which(names(train) %!in% "geneID" )]

    temp <- paste0("number of remarkable genes is ", nrow(train))
    message(temp)
    temp <- paste0("number of unremarkable genes is ", nrow(test))
    message(temp)

    if(model == "knn"){

        message("performing knn...\n it may take time")

        if(is.null(kn)){
            gg <- geneLab[complete.cases(geneLab), ]
            numClus <- length(unique(gg$label))
            up <- 20 + 2*numClus
            up <- min(up, nrow(specExp))
            kn <- seq(20, min(30, up))  }

        semiSuper <-
        train(label~., method = "knn", tuneGrid = expand.grid(k = kn),
        metric = "Accuracy", data = train)

        semiSuper
        show(summary(semiSuper))
        semiSuper$results
        show(semiSuper$finalModel)

        prediction <- predict(semiSuper, newdata = test)
        message("class assignment for unremarkable genes...")
        message("top class labels, and bottom number of predicted genes")
        show(table(prediction))

    }else if(model == "lr"){

        message("performing losgistic regression...")
        semiSuper <- train(label ~., method = "multinom", data = train)

        semiSuper
        show(summary(semiSuper))
        semiSuper$results
        show(semiSuper$finalModel)

        prediction <- predict(semiSuper, newdata = test)
        message("class assignment for unremarkable genes...")
        message("top class labels, and bottom number of predicted genes ")
        show(table(prediction))

    }

    test <- data.frame(test, label = prediction)

    df <- as.data.frame(rbind(train, test))
    df <- data.frame("geneID" = rownames(df), "label" = df$label)
    FinalLabeling<- inner_join(geneLab, df, by = "geneID")
    colnames(FinalLabeling) <- c("geneID", "semiLabel", "FinalLabel")

    newList <- list("semiSupervised" = semiSuper, "prediction" = prediction,
                    "FinalLabeling" = FinalLabeling)

    message("semi-supervised done!..\n")
    return(newList)

}

