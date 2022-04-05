
# getSymbols() -----------------------------------------------------------

#' Function to retrieve symbols from formula-formatted input
#'
#' Formula like strings are parsed, and returned in form of a list
#'
#' @param char input
#' @param exclude parts of the formula which will be excluded from output
#'
#' @return data frame with columns "name", "time", "value" and other
#' columns describing the measurements.
#'
#' @importFrom utils getParseData
#'
#' @noRd
getSymbols <- function(char, exclude = NULL) {
    ## Parse input data
    char <- char[char != "0"]
    out <- parse(text = char, keep.source = TRUE)
    out <- utils::getParseData(out)

    ## Split the input string at non-caracter symbols
    names <- unique(out$text[out$token == "SYMBOL"])

    ## Remove strings passed as "exclude"
    if (!is.null(exclude)) {
        names <- names[!names %in% exclude]
    }
    return(names)
}




# identifyEffects() ------------------------------------------------------

#' Method to identify the names and values of passed to biological, scaling and
#' error parameters
#'
#' The formula-formatted input will be parsed by \link{getSymbols} and the
#' respective output then passed back
#'
#' @param biological input for biological parameters
#' @param scaling input for scaling parameters
#' @param error input for error parameters
#'
#' @return two lists: values and parameter names
#'
#' @noRd

identifyEffects <- function(biological = NULL, scaling = NULL, error = NULL) {

    ## Get biological scale and error parameters
    biologicaValues <- getSymbols(as.character(biological)[3])
    scalingValues <- getSymbols(as.character(scaling)[3])
    errorValues <- getSymbols(as.character(error)[3])


    # Determine to which class parameters belong
    biologicalPars <- getSymbols(as.character(biological)[2])
    scalingPars <- getSymbols(as.character(scaling)[2])
    errorPars <- getSymbols(as.character(error)[2])

    effectsValues <- list(
        biologicaValues = biologicaValues,
        scalingValues = scalingValues,
        errorValues = errorValues
    )
    effectsPars <- list(
        biologicalPars = biologicalPars,
        scalingPars = scalingPars,
        errorPars = errorPars
    )

    return(list(effectsValues = effectsValues, effectsPars = effectsPars))
}



# replaceSymbols ----------------------------------------------------------

#' Method for replacing parts of a string. Used to paste expresisons for the
#' (error-) model.
#'
#' The formula-formatted input will be parsed by \link{getSymbols} and the
#' respective output then passed back
#'
#' @param what string or vector of strings to be replaced
#' @param by string or vector of strings that will take the replaced place
#' @param x string in which \code{what} will be replaced by \code{by}
#'
#' @return string, \code{x} with the replaced object.
#'
#' @importFrom utils getParseData
#'
#' @noRd
#'
replaceSymbols <- function(what, by, x) {
    xOrig <- x
    notZero <- which(x != "0")
    x <- x[notZero]
    useNames <- names(x)
    xParsed <- parse(text = x, keep.source = TRUE)
    data <- utils::getParseData(xParsed)
    by <- rep(by, length.out = length(what))
    names(by) <- what
    data$text[data$text %in% what] <- by[data$text[data$text %in% what]]
    data <- data[data$token != "expr", ]
    breaks <- c(0, which(diff(data$line1) == 1), length(data$line1))
    out <- lapply(
        seq_len((length(breaks) - 1)),
        function(i) {
            paste(
                data$text[seq((breaks[i] + 1), (breaks[i + 1]))],
                collapse = ""
            )
        }
    )
    names(out) <- useNames
    out <- unlist(out)
    xOrig[notZero] <- out
    return(xOrig)
}



# paste_() ----------------------------------------------------------------

#' Just a shorthand
#' @param ... the thing that should be pasted.
#'
#' @noRd
paste_ <- function(...) paste(..., sep = "_")




# analyzeBlocks() --------------------------------------------------------

#' Method to analyze which elements of the given matrix with <number of unique
#' set of scaling parameters> columns, and <unique set of biological
#' parameters> rows, are on the same scale (same block) and pass a list of lists
#' of the respective indices.
#' Each row describes a unique set of biological conditions, e.g. a
#' specific target measured under a specific condition at one time point.
#' The columns describe the sets of scaling parameters as target name and
#' gel (in case of western blot).
#'
#'
#' @param blockMatrix Input matrix, with entries equal to 1 wherever the set of
#' biological effects (row) is measured und the respective set of scaling
#' effects (column), i.e. each row has a one for each scaling under which it was
#' measured and each column has a one indicating which biological effects
#' where measured at the respective scaling.
#' All entries are either one or zero.
#'
#' @return List with one entry per scale. Each entry contains a list of row
#' indices of the sets of biological effects measured under the respective
#' scale.
#'
#' @noRd

analyzeBlocks <- function(blockMatrix) {
    out <- which(apply(blockMatrix, 1, sum) == 0)
    if (length(out) > 0) {
        blockMatrix <- blockMatrix[-out, ]
        cat("matrix contains zero rows which have been eliminated\n")
    }

    noUniqueBiological <- dim(blockMatrix)[1]
    rComponents <- list()
    cComponents <- list()

    counter <- 0
    while (length(unlist(rComponents)) < noUniqueBiological) {
        counter <- counter + 1

        if (length(unlist(rComponents)) == 0) {
            w <- 1
        } else {
            my_sample <- (1:noUniqueBiological)[-unlist(rComponents)]
            w <- min(my_sample)
        }

        repeat {
            v <- unique(
                rapply(
                    as.list(w),
                    function(i) which(blockMatrix[i, ] == 1)
                )
            )
            wNew <- unique(
                rapply(
                    as.list(v),
                    function(j) which(blockMatrix[, j] == 1)
                )
            )
            if (length(wNew) == length(w)) break
            w <- wNew
        }
        rComponents[[counter]] <- w
        cComponents[[counter]] <- v
    }

    return(rComponents)
}



# splitForScaling()  ----------------------------------------------------


#' split_data
#'
#' Split data in independent blocks according to biological and scaling
#' variables as being defined for \link{alignReplicates}. Each block will be given an
#' individual scaling factor.
#'
#' @param data data frame with columns "name", "time", "value" and others
#' @param effectsValues two-sided formula, see \link{alignReplicates}
#' @param normalizeInput logical, if set to TRUE, the input data will be
#' normalized by dividing all entries belonging to one scaling factor by their
#' respective mean. This prevents convergence failure on some hardware when the
#' data for different scaling effects differ by to many orders of magnitude.
#' @return list of data frames
#'
#' @noRd
splitForScaling <- function(data,
                              effectsValues,
                              normalizeInput) {
    if (!"1" %in% colnames(data)) {
        data["1"] <- 1
        intercept <- FALSE
    } else {
        intercept <- TRUE
    }

    # Construnct strings containing the values of the respective effects
    scalingStrings <- Reduce(paste_, data[effectsValues[[2]]])
    biologicalStrings <- Reduce(paste_, data[effectsValues[[1]]])

    scalingStringsUnique <- unique(scalingStrings)
    biologicalStringsUnique <- unique(biologicalStrings)

    # Initialize matrix
    blockMatrix <- matrix(
        0,
        ncol = length(scalingStringsUnique),
        nrow = length(biologicalStringsUnique)
    )

    ## For every datapoint set the entry in the matrix to 1, corresponding to the
    ## index in the respective unique lists of fixed (row) and specific (col)
    for (i in seq_len(nrow(data))) {
        useRow <- which(biologicalStringsUnique == biologicalStrings[i])
        useCol <- which(scalingStringsUnique == scalingStrings[i])
        blockMatrix[useRow, useCol] <- 1
    }

    ## compile list of scalings
    scalingsList <- analyzeBlocks(blockMatrix)

    ## Remove the "1"-column added above
    if (!intercept) {
        data <- data[, -which(colnames(data) == "1")]
    }

    ## Compile the list
    outList <- lapply(
        scalingsList,
        function(l) {
            data[biologicalStrings %in% biologicalStringsUnique[l], ]
        }
    )

    ## normalize the data
    if (normalizeInput) {
        for (i in seq_len(length(outList))) {
            outList[[i]]$value <- outList[[i]]$value /
                mean(outList[[i]]$value)
            ## give a warning if division by zero is possible
            if (min(outList[[i]]$value) < 0) {
                warning(
                    paste0(
                        "'normalizeInput == TRUE' and negative input: ",
                        "division by zero possible."
                    )
                )
            }
        }
    }

    return(outList)
}



# inputCheck() -----------------------------------------------------------

#' Check input parameters for structural errors
#'
#' the input of \link{alignReplicates} is checked for user mistakes.
#'
#' @param data Input data, usually output of \link{readWide}
#' @param model Model definition as a string
#' @param errorModel Error model definition as a string
#' @param biological Definition of the biological effects
#' @param scaling Definition of the scaling effects
#' @param error Definition of the error effects
#' @param fitLogscale logical, defines the parameter fit scale
#'
#' @noRd

inputCheck <- function(data = NULL,
                        model = NULL,
                        errorModel = NULL,
                        biological = NULL,
                        scaling = NULL,
                        error = NULL,
                        fitLogscale = NULL,
                        outputScale = NULL) {

    ## Check if parameters are present
    if (is.null(model) | is.null(errorModel) | is.null(biological) |
        is.null(scaling) | is.null(error)) {
        stop(
            "All of model, errorModel, biological, scaling, error ",
            "must be set."
        )
    }

    ## check data
    columnNames <- names(data)
    expectedNames <- union(
        getSymbols(as.character(biological)[3]),
        union(
            getSymbols(as.character(scaling)[3]),
            getSymbols(as.character(error)[3])
        )
    )

    if (!all(expectedNames %in% columnNames)) {
        stop(
            "Not all column names set in 'biological', 'scaling' and ",
            "'error' are present in 'data'."
        )
    }

    ## Check scaling
    if (!is.logical(fitLogscale)) {
        stop(
            "'fitLogscale' must be logical"
        )
    }



    ## Stop if formulas have the wrong specification
    if (length(as.character(biological)) == 1 |
        length(as.character(scaling)) == 1 | length(as.character(error)) == 1) {
        stop("Do not pass biological, scaling or error as string.")
    }
    if (length(as.character(biological)) < 3) {
        stop("Left and right-hand side of formula 'biological' is needed")
    }
    if (length(as.character(scaling)) < 3) {
        stop("Left and right-hand side of formula 'scaling' is needed")
    }
    if (length(as.character(error)) < 3) {
        stop("Left and right-hand side of formula 'error' is needed")
    }

    returnMsg <- "All input checks passed."
    return(returnMsg)
}



# scaleTarget() ----------------------------------------------------------

#' Scaling function applied to the current dataset
#'
#' Heart of the package. The passed set is assumed to belong to one scale. The
#' respective scaling factors are calculated by optimizing
#' \link{objFunction}.
#'
#' @param currentData current data set, belongs to one scale.
#' @param passParList from the the \code{passParList} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{effectsValues}}{
#'      Named list of vectors. The names correspond to the effects and the
#'      vectors contain the 'values'; the right hand side of the respective
#'      effects parameters passed to \link{alignReplicates} the names of the columns
#'      associated with the respective effect.
#'  }
#'  \item{\code{parameterData}}{
#'      If the \code{data} parameter passed to \link{alignReplicates} is already an
#'      output of \link{alignReplicates}, it contains the fitted parameters, otherwise
#'      \code{NULL}.
#'  }
#'  \item{\code{averageTechRep}}{
#'      Logical, if \code{TRUE} technical replicates, that can not be separated
#'      by the respective \code{biological} and \code{scaling} values will be
#'      averaged.
#'  }
#'  \item{\code{verbose}}{
#'      Logical, if \code{TRUE} additional output will be generated.
#'  }
#'  \item{\code{covariates}}{
#'      String, the covariates as defined by the (error-) model expressions as
#'      defined in the respective parameters of \link{alignReplicates}
#'  }
#'  \item{\code{parameters}}{
#'      Named vector, contains the variables for the three effects, and the
#'      respective effects as names.
#'  }
#'  \item{\code{fitLogscale}}{
#'      logical, defining if the parameters are fitted on a log scale.
#'  }
#'  \item{\code{effectsPars}}{
#'      Similar to \code{parameters}, but as a named list.
#'  }
#'  \item{\code{modelExpr}}{
#'      What is passed to \code{model} in \link{alignReplicates} but with the entries
#'      of \code{parameters} replaced by the respective name (e.g. 'yi' is
#'      replaced by 'biological'). The result is parsed as an expression.
#'  }
#'  \item{\code{errorModelExpr}}{
#'      Same as \code{modelExpr} but with the value of \code{errorModel} from
#'      \link{alignReplicates}.
#'  }
#'  \item{\code{constrExpr}}{
#'      If \code{normalize} is set to \code{TRUE} in \link{alignReplicates},
#'      \code{constrExpr} contains expression of the normalization
#'      constraint, otherwise it is empty.
#'  }
#'  \item{\code{modelJacobianExpr}}{
#'      Named list with the derivatives of \code{modelExpr} with respect to the
#'      respective effects as entries.
#'  }
#'  \item{\code{errorModelJacobianExpr}}{
#'      Analogue to \code{modelJacobianExpr} but with \code{errorModelExpr}.
#'  }
#'  \item{\code{constStrength}}{
#'      Numerical 1000, if \code{normalize} is set to \code{TRUE} in
#'      \link{alignReplicates}, otherwise unused.
#'  }
#'  \item{\code{normalize}}{
#'      Logical, direct taken from the \code{normalize} parameter of
#'      \link{alignReplicates}.
#'  }
#' }
#'
#' @return A list of \code{data.frame}S with the entries:
#' \describe{
#'  \item{\code{outPrediction}}{
#'      \code{data.frame} with the columns \code{name}, \code{time},
#'      \code{value}, \code{sigma}, \code{biological}, \code{scaling} and
#'      \code{error}.
#'
#'      \code{values} contains the predictions by evaluation of the (error-)
#'      model with the fitted parameters on the original scale (if
#'      \code{normalizeInput} is set to \code{FALSE} in \link{alignReplicates}). The
#'      last three columns contain strings pasted from the respective values of
#'      the current entry.
#'  }
#'  \item{\code{outScaled}}{
#'      \code{data.frame} with the original columns and additionally
#'      \code{sigma} and \code{1}.
#'
#'      The \code{values} are the original ones scaled to common scale by
#'      applying the inverse of the model with the fitted parameters. The
#'      \code{sigma} entries are the results of the evaluated error model scaled
#'      to common scale adhering to Gaussian error propagation.
#'  }
#'  \item{\code{outAligned}}{
#'      \code{data.frame} with the original column names, the column
#'      \code{value} contains the estimated biological parameters.
#'
#'      The \code{values} are estimated biological parameters, while the errors
#'      \code{sigma} are from the Fisher information.
#'  }
#'  \item{\code{currentData}}{
#'      \code{data.frame} with the original data with added columns \code{sigma}
#'      and \code{1}.
#'  }
#'  \item{\code{outOrigParameters}}{
#'      Same \code{data.frame} as \code{currentData} but with added columns for
#'      the \code{parameters} ad the respective fitted values.
#'  }
#'  \item{\code{parTable}}{
#'      \code{data.frame} with the columns:
#'      \itemize{
#'          \item{\code{level}:}{
#'          String pasting all the unique effect entries.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{parameter}:}{
#'          Parameter associated with the current effect (compare with
#'          \code{parameters}).
#'          }
#'      }
#'      \itemize{
#'          \item{\code{value}:}{
#'          Value of the estimated parameters. The ones corresponding to the
#'          biological parameters coincide with the \code{values} column of
#'          \code{outAligned}.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{sigma}:}{
#'          Errors of the estemated parameter values by utilizing the Fisher
#'          information
#'          }
#'      }
#'      \itemize{
#'          \item{\code{nll}:}{
#'          Negative of twice the log likelihood of the fit in which the
#'          paramters are estimated.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{noPars}:}{
#'          Number of fitted parameters (biological, scaling and error), minus
#'          one if \code{normalize = TRUE} in the call of \link{alignReplicates}.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{noData}:}{
#'          Length of the current data file, i.e. number of measurements in the
#'          set that is currently scaled.
#'          }
#'      }
#'  }
#' }
#'
#' @importFrom rootSolve multiroot
#' @importFrom trust trust
#' @importFrom MASS ginv
#'
#' @noRd

scaleTarget <- function(currentData,
                         passParList) {

    ## Retrieve parameters from list
    effectsValues <- passParList$effectsValues
    parameterData <- passParList$parameterData
    averageTechRep <- passParList$averageTechRep
    verbose <- passParList$verbose
    covariates <- passParList$covariates
    parameters <- passParList$parameters
    fitLogscale <- passParList$fitLogscale
    effectsPars <- passParList$effectsPars
    normalize <- passParList$normalize
    outputScale <-  passParList$outputScale
    iterlim <- passParList$iterlim


    currentName <- unique(currentData$name)

    ## Add column for error
    currentData$sigma <- NaN

    ## Add a dummy column filled with 1
    currentData[["1"]] <- "1"

    ## Initialize the parameter for the current target
    currentParameter <- NULL

    if (!is.null(parameterData)) {
        currentParameter <- parameterData[
            parameterData$name %in% currentData$name,
        ]
    }


    ## Build a list of string of biological/scaling values present
    dataFitBiological <- do.call(
        paste_, currentData[, effectsValues[[1]], drop = FALSE]
    )

    dataFitScaling <- do.call(
        paste_, currentData[, effectsValues[[2]], drop = FALSE]
    )

    ## Reducing datapoints
    if (averageTechRep) {
        if (verbose) {
            cat("Analyzing technical replicates ... ")
        }
        groups <- interaction(dataFitBiological, dataFitScaling)
        if (any(duplicated(groups))) {
            currentData <- do.call(rbind, lapply(
                unique(groups),
                function(g) {
                    subdata <- currentData[groups == g, ]
                    outdata <- subdata[1, ]
                    outdata$value <- mean(subdata$value)
                    return(outdata)
                }
            ))
            cat(
                "data points that could not be distinguished by either ",
                "biological\nor scaling variables have been averaged.\n"
            )
        } else {
            if (verbose) {
                cat("none found.\n")
            }
        }
    }

    ## compile dataset fit for fitting
    dataFit <- data.frame(
        currentData[
            ,
            union(c("name", "time", "value", "sigma"), covariates)
        ],
        biological = do.call(
            paste_,
            currentData[, effectsValues[[1]], drop = FALSE]
        ),
        scaling = do.call(
            paste_,
            currentData[, effectsValues[[2]], drop = FALSE]
        ),
        error = do.call(
            paste_,
            currentData[, effectsValues[[3]], drop = FALSE]
        ),
        stringsAsFactors = FALSE
    )

    ## Retrieve (unique) list of biological, scaling and error values
    levelsList <- list(
        biological = unique(as.character(dataFit$biological)),
        scaling = unique(as.character(dataFit$scaling)),
        error = unique(as.character(dataFit$error))
    )

    ## Combine above lists allLevels will therefore have the length of the
    ## sum of all different biological, scaling and error values (in that order)
    allLevels <- unlist(
        lapply(
            seq_along(parameters),
            function(k) {
                switch(names(parameters)[k],
                    biological = levelsList[[1]],
                    scaling = levelsList[[2]],
                    error = levelsList[[3]]
                )
            }
        )
    )
    #
    initPars <- genInitPars(
        parameters,
        fitLogscale,
        levelsList
    )
    #
    mask <- genMask(
        initPars,
        parameters,
        allLevels,
        dataFit
    )

    fitParsBiological <- NULL

    if (verbose) {
        cat("Starting fit\n")
    }

    passParList2 <- list(
        dataFit = dataFit,
        levelsList = levelsList,
        fitParsBiological = fitParsBiological,
        effectsPars = effectsPars,
        mask = mask,
        initPars = initPars
    )

    # * call of trust() function ----------------------------------------------


    fitResult <- trust::trust(
        objfun = objFunction,
        parinit = initPars,
        rinit = 1,
        rmax = 10,
        iterlim = iterlim,
        blather = verbose,
        passParList = passParList,
        passParList2 = passParList2
    )

    if (!fitResult$converged) {
        warning(paste("Non-converged fit for target", currentName))
    } else if (verbose == TRUE) {
        cat("Fit converged.\n")
    }

    residualsFit <- resolveFunction(
        currentPars = fitResult$argument,
        passParList = passParList,
        passParList2 = passParList2,
        calcDeriv = FALSE
    )

    bessel <- sqrt(
        nrow(dataFit) / (nrow(dataFit) - length(initPars) +
            normalize)
    )

    ## Get singular values (roots of non negative eigenvalues of M^* \cdot M)
    ## here it is synonymous with eigenvalue.
    singleValues <- svd(fitResult[["hessian"]])[["d"]]

    ## Define a tollerance threshold, the root of the machine precission is
    ## a usual value for this threshold.
    tol <- sqrt(.Machine$double.eps)

    ## Define nonidentifiable as being "to small to handle", judged by the
    ## above defined threshold
    nonIdentifiable <- which(singleValues < tol * singleValues[1])
    if (length(nonIdentifiable) > 0) {
        warning("Eigenvalue(s) of Hessian below tolerance. Parameter
                  uncertainties might be underestimated.")
    }

    # * parameter table -------------------------------------------------------
    # Generate parameter table

    parTable <- data.frame(
        name = currentName,
        level = c(
            rep(levelsList[[1]], length(effectsPars[[1]])),
            rep(levelsList[[2]], length(effectsPars[[2]])),
            rep(levelsList[[3]], length(effectsPars[[3]]))
        ),
        parameter = names(fitResult$argument),
        value = fitResult$argument,



        nll = fitResult$value,
        noPars = length(fitResult$argument) - normalize,
        noData = nrow(dataFit)
    )

    ## Calculating the error from the inverse of the Fisher information
    ## matrix which is in this case the Hessian, to which the above
    ## calculated Bessel correction is applied.
    parTable$sigma <- as.numeric(
      sqrt(diag(2 * MASS::ginv(fitResult$hessian)))
    ) * bessel




    ## transform parameters back if fitted on log scale
    if (fitLogscale == TRUE) {
        parTable$value <- exp(parTable$value)
        parTable$sigma <- parTable$value * parTable$sigma
        parTable$lower <- parTable$value -
          parTable$sigma
        parTable$upper <- parTable$value +
          parTable$sigma




    }

    if (verbose) {
        cat("Estimated parameters on non-log scale:\n")
        print(parTable)
        cat(
            "converged:", fitResult$converged, ", iterations:",
            fitResult$iterations, "\n"
        )
        cat("-2*LL: ", fitResult$value, "on", nrow(currentData) +
            normalize - length(fitResult$argument), "degrees of freedom\n")
    }

    attr(parTable, "value") <- fitResult$value
    attr(parTable, "df") <- nrow(currentData) + normalize -
        length(fitResult$argument)


    # * outPrediction --------------------------------------------------------
    ## Predicted data
    outPrediction <- currentData
    outPrediction$value <- fitResult$prediction
    outPrediction$sigma <- fitResult$sigma * bessel
    outPrediction$lower <- outPrediction$value - outPrediction$sigma
    outPrediction$upper <- outPrediction$value + outPrediction$sigma





    # * outScaled ------------------------------------------------------------
    ## Initialize list for scaled values
    initialValuesScaled <- rep(0, nrow(currentData))
    if (verbose) {
        cat("Inverting model ... ")
    }


    valuesScaled <- try(
        rootSolve::multiroot(
            f = evaluateModel,
            start = initialValuesScaled,
            jacfunc = evaluateModelJacobian,
            parList = residualsFit[-1],
            verbose = FALSE,
            passParList = passParList,
            passParList2 = passParList2
        )$root,
        silent = TRUE
    )
    if (verbose) {
        cat("done\n")
    }

    if (inherits(valuesScaled, "try-error")) {
        outScaled <- NULL
        warning("Rescaling to common scale not possible.
                  Equations not invertible.")
    } else {
        myDeriv <- abs(
            evaluateModelJacobian(
                valuesScaled,
                residualsFit[-1],
                passParList = passParList,
                passParList2 = passParList2
            )
        )
        sigmasScaled <- fitResult$sigma * bessel / myDeriv

        # transform parameters back if fitted on log scale
        if (fitLogscale == TRUE) {
            valuesScaled <- exp(valuesScaled)
            sigmasScaled <- valuesScaled * sigmasScaled
        }
        upperScaled <- valuesScaled + sigmasScaled
        lowerScaled <- valuesScaled - sigmasScaled

        outScaled <- currentData
        outScaled$value <- valuesScaled
        outScaled$sigma <- sigmasScaled
        outScaled$upper <- upperScaled
        outScaled$lower <- lowerScaled
    }


    # * outAligned -----------------------------------------------------------
    noInitial <- length(levelsList[[1]])

    ## Use one datapoint per unique set of fixed parameters
    outAligned <- currentData[
        !duplicated(dataFit$biological),
        intersect(effectsValues[[1]], colnames(currentData))
    ]

    ## The values are the fitted parameters for the respective fixed
    ## parameter ensembles.
    outAligned$value <- residualsFit[[effectsPars[[1]][1]]]
    outAligned$parameter <- c(
        do.call(
            rbind,
            lapply(
                seq_len(nrow(outAligned)),
                function (i) {
                    out <- parTable[
                        which(parTable$level == do.call(
                            paste_,
                            c(outAligned[i,effectsValues[[1]]])
                        )
                        ),
                    ]$parameter
                }
            )
        )
    )

    outAligned$sigma <- as.numeric(
      sqrt(diag(2 * MASS::ginv(fitResult$hessian)))
    )[seq_len(noInitial)] * bessel


    # transform parameters back if fitted on log scale
    if (fitLogscale == TRUE) {
        outAligned$value <- exp(outAligned$value)
        outAligned$sigma <- outAligned$value * outAligned$sigma
        outAligned$lower <- outAligned$value -
          outAligned$sigma
        outAligned$upper <- outAligned$value +
          outAligned$sigma
    }


    # * outOrigParameters -------------------------------------------------
    outOrigParameters <- currentData
    for (k in seq_along(parameters)) {
        effect <- names(parameters)[k]
        useLevels <- as.character(dataFit[[effect]])
        index0 <- which(
            as.character(
                gsub('[_0-9]+', '', parTable$parameter)
                ) == parameters[k]
        )
        index1 <- match(useLevels, as.character(parTable$level[index0]))
        index <- index0[index1]
        outOrigParameters[[parameters[k]]] <- parTable$value[index]
    }
    out <- list(
        outPrediction = outPrediction,
        outScaled = outScaled,
        outAligned = outAligned,
        original = currentData,
        outOrigParameters = outOrigParameters,
        parTable = parTable
    )




    return(out)
}



# genInitPars() -------------------------------------------------

#' Method to generate a set of initial parameters for \link{scaleTarget}
#'
#' @param parameters Named vector, contains the variables for the three effects,
#' and the respective effects as names.
#' @param fitLogscale logical, identifying if the parameters are
#' fited on log scale
#' @param levelsList Named list with one entry per effect containing a vector
#' of strings composed with the respective pasted entries of the effects columns
#'
#' @return Named vector with one entry per parameter, all are 1 if
#' \code{fitLogscale = TRUE} and 0 if not. The names are the entries
#' of \code{parameters} of the corresponding effect. Each entry is repeated as
#' often as there are parameters of the respective effect.
#'
#' @noRd

genInitPars <- function(parameters,
                                  fitLogscale,
                                  levelsList) {
    ## Create a named vector with initial values for all parameters. The name
    ## is the parameter, and the initial value is 0 for log and 1 for linear.
    ## each parameter name-value entry is repeated as many times as there are
    ## measurements of the i'th target with this parameter.
    ## Example: biological is passed as "ys ~ condition", so the entry with name
    ##   "ys" will be repeated for as many times as there are measurements
    ##   with the same "condition" entry.
    initPars <- do.call(
        "c",
        lapply(
            seq_along(parameters),
            ## Go through all parameters with current parameter n
            function(n) {
                if (fitLogscale == TRUE) {
                    v <- 0
                } else {
                    v <- 1
                }
                ## Let l be the number of measurements of same type (biological,
                ## scaling or error) for the current n
                l <- switch(
                    ## Get "biological", "scaling" or "error" from the current n
                    names(parameters)[n],
                    biological = length(levelsList[[1]]),
                    scaling = length(levelsList[[2]]),
                    error = length(levelsList[[3]])
                )

                ## Get the n'th parameter
                p <- as.character(parameters[n])

                ## The variable "parValues" will have l times the entry
                ## "v" (0 for log, 1 for linear) with the corresponding
                ## parameter as the name.
                parValues <-
                    structure(
                        rep(v, l),
                        names = rep(p, l)
                    )
                return(parValues)
            }
        )
    )

    ## give intividual names
    names(initPars) <- paste(
        names(initPars),
        seq_along(initPars),
        sep = "_"
        )

    return(initPars)
}


# genMask() ---------------------------------------------------------

#' Creates a identifier list
#'
#' A list with one element per entry in initPars is generated. For
#' each of those elements, it is checked which elements of the data set
#' \code{dataFit} coincides with the current element of \code{allLevels}.
#' The entry is a list with 1 indicating where the current element of
#' \code{allLevels} is found in the respective effect-column of
#' \code{dataFit}.
#' Keep in mind, that the \code{initPars} and \code{allLevels} are
#' structured analogue.
#'
#' @param initPars named vector, set of parameters as output of
#' \link{genInitPars}
#' @param parameters named vector, with the 'values' of the three effects
#' @param allLevels character vector, strings describing all levels in the
#' current data set
#' @param dataFit data frame containing the data that will be fitted
#'
#' @return list with one entry per entry of \code{allLevels}. Each of this
#' entries is a numerical list of the length of \code{dataFit}. The entries are
#' either 1, if the corresponding entry of \code{dataFit} is resembled by the
#' current entry of \code{allLevels}.
#'
#' @noRd

genMask <- function(initPars,
                          parameters,
                          allLevels,
                          dataFit) {
    mask <- lapply(
        seq_along(initPars),
        function(k) {
            effect <- getParClass(
                parameters[
                    match(
                        getParClass(initPars)[k],
                        parameters
                    )
                ]
            )
            mask_vector <- as.numeric(
                dataFit[[effect]] == allLevels[k]
            )
            return(mask_vector)
        }
    )

    return(mask)
}



# resolveFunction() -----------------------------------------------------

#' Function to sort the current set of parameters into the three effect
#' categories.
#'
#' Additionally, residuals of model evaluations for the set of
#' \code{currentPars} \code{fitParsBiological} are calculated by
#' evaluation of \link{rssModel}.
#'
#' @param currentPars named vector of vectors to be tested currently
#' @param passParList from the the \code{passParList} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{parameters}}{Named list of the parameters of the three effects,
#'  the names are the respective effects}
#' }
#' @param passParList2 from the the \code{passParList2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{fitParsBiological}}{TODO: Check if it is always \code{NULL}?}
#'  \item{\code{levelsList}}{Named list of vectors. One entry per effect with
#'  the respective name. The entry contains a list of the unique strings
#'  composed from the entries of \code{effect_values} of the respective effect,
#'  i.e. the entries of the columns containing e.g. the biological-effects
#'  (name, time, condition etc.).}
#'  \item{\code{effectsPars}}{Named list of vectors. One entry per effect with
#'  the respective name. The entry then contains the string of the effect
#'  parameter, i.e. the variable name (e.g. "sj" for scaling).}
#' }
#'
#' @return named list of the residuals calculated in \link{rssModel} and the
#' \code{currentPars} sorted in the effects with the respective entries
#' from \code{levelsList} as names.
#'
#' @noRd
#'
resolveFunction <- function(currentPars,
                             passParList,
                             passParList2,
                             calcDeriv) {
    if (FALSE) {
        currentPars <- initPars
    }


    parameters <- passParList$parameters
    fitParsBiological <- passParList2$fitParsBiological
    levelsList <- passParList2$levelsList
    effectsPars <- passParList2$effectsPars

    parsAll <- c(currentPars, fitParsBiological)
    parList <- lapply(
        parameters,
        function(n) {
            parClass <- gsub('[_0-9]+', '', names(parsAll))
            subpar <- parsAll[parClass == n]
            if (n %in% effectsPars[[1]]) {
                ## Rename entries of the first list (fixed)
                names(subpar) <- levelsList[[1]]
            }
            if (n %in% effectsPars[[2]]) {
                ## Rename entries of the second list (latent)
                names(subpar) <- levelsList[[2]]
            }
            if (n %in% effectsPars[[3]]) {
                ## Rename entries of the third list (error)
                names(subpar) <- levelsList[[3]]
            }
            return(subpar)
        }
    )

    ## Rename the thre lists
    names(parList) <- parameters

    ## Create list of residuals and attach the corresponding current
    ## parameters
    c(
        list(
            residuals = rssModel(
                parList,
                passParList,
                passParList2,
                calcDeriv = calcDeriv
            )
        ),
        parList
    )
}



# rssModel() -------------------------------------------------------------

#' Calculate model residuals
#'
#' Residuals are calculated by the difference (prediction - value), where
#' \code{prediction} is the evaluation of the model with the current set of
#' parameters and \code{value} the respective measurements. Optionally, the
#' derivatives are also calculated
#'
#'
#' @param parList Named list with one entry per effect (with the respective
#' name). The entry contains a named vector with the unique stings of the
#' respective effect values i.e. the entries of the respective columns (e.g.
#' name, time, condition etc. for 'biological')
#' @param passParList from the the \code{passParList} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{modelExpr}}{The argument of the \code{model} parameter of
#'  \link{alignReplicates} parsed as an executable expression.}
#'  \item{\code{errorModelExpr}}{The argument of the \code{errorModel}
#'  parameter of \link{alignReplicates} parsed as en executable expression.}
#'  \item{\code{constrExpr}}{If the logical argument \code{normalize} of
#'  \link{alignReplicates} is set to \code{TRUE}, an executable constraint expression
#'  is passed. If \code{normalize = FALSE} it is empty.}
#'  \item{\code{modelJacobianExpr}}{The Jacobian of the model as a named list
#'  with the respective derivatives as entries passed as expression.}
#'  \item{\code{errorModelJacobianExpr}}{The Jacobian of the error model as a
#'  named list with the respective derivatives as entries passed as expression.}
#' }
#'
#' @param passParList2 from the the \code{passParList2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{dataFit}}{The current dataset}
#' }
#'
#' @param calcDeriv logical, if \code{TRUE}, derivatives will also be
#' calculated.
#'
#' @return list of residuals with the derivatives (optional) as attributes.
#'
#' @noRd
rssModel <- function(parList,
                      passParList,
                      passParList2,
                      calcDeriv = TRUE) {
    modelExpr <- passParList$modelExpr
    errorModelExpr <- passParList$errorModelExpr
    constrExpr <- passParList$constrExpr
    modelJacobianExpr <- passParList$modelJacobianExpr
    errorModelJacobianExpr <- passParList$errorModelJacobianExpr

    dataFit <- passParList2$dataFit

    with(
        ## paste list with entries: name, time, value, sigma, the effects and
        ## their paramters
        c(as.list(dataFit), parList), {
          if(F) {
            name <- c(as.list(dataFit), parList)$name
            time <- c(as.list(dataFit), parList)$time
            value <- c(as.list(dataFit), parList)$value
            sigma <- c(as.list(dataFit), parList)$sigma
            biological <- c(as.list(dataFit), parList)$biological
            scaling <- c(as.list(dataFit), parList)$scaling
            error <- c(as.list(dataFit), parList)$error
            yi <- c(as.list(dataFit), parList)$yi
            sj <- c(as.list(dataFit), parList)$sj
            sigmaR <- c(as.list(dataFit), parList)$sigmaR
            sigmaA <- c(as.list(dataFit), parList)$sigmaA
          }
            ## Generate lists var and prediction with <NoOfMeasurements> entries
            ## and initialize the it with 1
            prediction <- var <- rep(1, length(value))

            ## Get the prediction by evaluating the model
            prediction[seq_along(prediction)] <- eval(modelExpr)

            ## Create residuals: differences between prediction and measurements
            res <- prediction - value

            ## Generate variances by squaring the evaluated error model
            var[seq_along(var)] <- eval(errorModelExpr)^2

            ## Evaluate constraints
            constr <- eval(constrExpr)

            ## Create list of residuals
            residuals <- c(res / sqrt(var), var, constr)

            ## Initialize variables for derivatives
            residualDeriv <- varianceDeriv <- NULL

            ## Calculate derivatives if wished
            if (calcDeriv) {
                jac <- lapply(
                    seq_len(length(parList)),
                    function(k) {
                        ## Initialize lists
                        vMod <- vErr <- rep(0, length(value))

                        ## Get derivatives for model and errormodel
                        vMod[seq_len(length(value))] <- eval(
                            modelJacobianExpr[[k]]
                        )
                        vErr[seq_len(length(value))] <- eval(
                            errorModelJacobianExpr[[k]]
                        )

                        residualDeriv <- vMod / sqrt(var) - vErr * res / var
                        varianceDeriv <- vErr * 2 * sqrt(var)
                        list(residualDeriv, varianceDeriv)
                    }
                )

                residualDeriv <- lapply(jac, function(j) j[[1]])
                varianceDeriv <- lapply(jac, function(j) j[[2]])
            }
            attr(residuals, "residualDeriv") <- residualDeriv
            attr(residuals, "varianceDeriv") <- varianceDeriv

            return(residuals)
        }
    )
}


# objFunction() ----------------------------------------------------

#' Objective function for trust region optimizer
#'
#' @param currentPars Current set of parameters to be tested. Will be
#' iteratively changed by the objective function.
#' @param passParList from the the \code{passParList} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{parameters}}{
#'      Named vector, contains the variables for the three effects, and the
#'      respective effects as names.
#'  }
#'  \item{\code{fitLogscale}}{
#'      Logical, defining the parameter fit scale. See
#'      \code{fitLogscale} in \link{alignReplicates}.
#'  }
#'  \item{\code{constStrength}}{
#'      Numerical 1000, if \code{normalize} is set to \code{TRUE} in
#'      \link{alignReplicates}, otherwise unused.
#'  }
#' }
#'
#' @param passParList2 from the the \code{passParList2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{dataFit}}{
#'      The current data set
#'  }
#'  \item{\code{levelsList}}{Named list of vectors. One entry per effect with
#'      the respective name. The entry contains a list of the unique strings
#'      composed from the entries of \code{effect_values} of the respective
#'      effect, i.e. the entries of the columns containing e.g. the biological
#'      effects (name, time, condition etc.).
#'  }
#'  \item{\code{mask}}{
#'      Output of \link{genMask}
#'  }
#' }
#'
#' @param calcDeriv Logical, indicates if derivatives should be also
#' calculated.
#'
#' @return named list with the entries \code{value}, \code{gradient} and
#' \code{hessian}. The \code{value}S are the residual sum of squares with the
#' residuals as calculated by \link{rssModel} (via \link{resolveFunction}).
#'
#' @noRd
objFunction <- function(currentPars,
                               passParList,
                               passParList2,
                               calcDeriv = TRUE) {
    if (FALSE) {
        currentPars <- initPars
        calcDeriv <- TRUE
    }

    parameters <- passParList$parameters
    fitLogscale <- passParList$fitLogscale
    constStrength <- passParList$constStrength

    dataFit <- passParList2$dataFit
    levelsList <- passParList2$levelsList
    mask <- passParList2$mask

    noData <- nrow(dataFit)

    ## Recover residuals from output of res_fn()
    calculatedResiduals <- resolveFunction(
        currentPars = currentPars,
        passParList = passParList,
        passParList2 = passParList2,
        calcDeriv = calcDeriv
    )$residuals

    ## Retrieve residuals of model, errormodel and constraint as well as
    ## the derivatives.
    residuals <- calculatedResiduals[1:noData]
    variances <- calculatedResiduals[seq((noData + 1), (2 * noData))]
    constraint <- calculatedResiduals[2 * noData + 1]
    residualDeriv <- attr(calculatedResiduals, "residualDeriv")
    varianceDeriv <- attr(calculatedResiduals, "varianceDeriv")

    ## Set bessel correction factor to 1
    bessel <- 1

    ## Calculate the current value
    value <- sum(residuals^2) + bessel * sum(log(2*pi * variances)) + constraint^2
    gradient <- NULL
    hessian <- NULL

    ## Get derivatives for the measurements by applying the predefined mask
    if (calcDeriv) {
        calculatedResidualsJacobian <- do.call(cbind, lapply(
            seq_along(currentPars),
            function(k) {

                ## Get the "class" (fixed, latent, or error) of the current k'th
                ## parameter
                whichPar <- match(getParClass(currentPars)[k], parameters)

                ## Apply the mask to the residuals
                residualJacobian <- residualDeriv[[whichPar]] * mask[[k]]
                varianceJacobian <- varianceDeriv[[whichPar]] * mask[[k]]
                constrainJacobian <- as.numeric(
                    getParClass(currentPars)[k] == parameters["biological"]
                ) * constStrength / length(levelsList[[1]])

                ## Convert to log if wanted
                if (fitLogscale == TRUE) {
                    constrainJacobian <- constrainJacobian *
                        exp(currentPars[k])
                }

                ## Stitch the results together
                c(residualJacobian, varianceJacobian, constrainJacobian)
            }
        ))

        ## Split the above list into corresponding parts
        residualJacobian <- calculatedResidualsJacobian[
            1:noData, ,
            drop = FALSE
        ]
        jacVars <- calculatedResidualsJacobian[
            (noData + 1):(2 * noData), ,
            drop = FALSE
        ]
        constrainJacobian <- calculatedResidualsJacobian[
            2 * noData + 1, ,
            drop = FALSE
        ]

        ## Compose to gradient vector and hessian matrix
        gradient <- as.vector(
            2 * residuals %*% residualJacobian + (bessel / variances) %*%
                jacVars + 2 * constraint * constrainJacobian
        )
        hessian <- 2 * t(rbind(residualJacobian, constrainJacobian)) %*%
            (rbind(residualJacobian, constrainJacobian))
    }

    ## Takes residual (ultimative res <- prediction - value and res / sqrt(var)
    ## from rssModel()) to retrive the presdiction
    prediction <- residuals * sqrt(variances) + dataFit$value

    ## Calculate errors
    sigma <- sqrt(variances)

    ## Compose all of the above to one list
    list(
        value = value, gradient = gradient, hessian = hessian,
        prediction = prediction, sigma = sigma
    )
}


# evaluateModel() --------------------------------------------------------

#' Model evaluation method
#'
#' Method to evaluate the model with a given set of parameters.
#'
#' @param initPars Parameters used for model evaluation
#'
#' @param parList Named list with one entry per effect (with the respective
#' name). The entry contains a named vector with the unique stings of the
#' respective effect values i.e. the entries of the respective columns (e.g.
#' name, time, condition etc. for 'biological')
#' @param passParList from the the \code{passParList} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{effectsPars}}{
#'      Named list with one entry per effect. The values are the respective
#'      effect parameters.
#'  }
#'  \item{\code{modelExpr}}{
#'      The argument of the \code{model} parameter of \link{alignReplicates} parsed as
#'      an executable expression.
#'  }
#' }
#'
#' @param passParList2 from the the \code{passParList2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{dataFit}}{The current dataset.}
#' }
#'
#'
#' @noRd

evaluateModel <- function(initPars,
                           parList,
                           passParList,
                           passParList2) {
    ## Create a list with with values from 1 to the number of parameters

    effectsPars <- passParList$effectsPars
    modelExpr <- passParList$modelExpr
    dataFit <- passParList2$dataFit




    biological <- seq_along(initPars)

    ## Generate a list with the entries:
    ##   parameters, the sequence saved generated as "fixed", name, time,
    ##   value, sigma, fixed, latent, error, ys, sj, sigmaR
    useList <- c(
        list(initPars, biological = biological),
        as.list(dataFit),
        parList
    )
    names(useList)[1] <- effectsPars[[1]][1]

    values <- with(useList, eval(modelExpr) - dataFit$value)

    return(values)
}



# evaluateModelJacobian() -----------------------------------------------

#' Model jacobian evaluation method
#'
#' Method to evaluate the model jacobian with a given set of parameters.
#'
#' @param initPars Parameters used for model jacobian evaluation
#'
#' @param parList Named list with one entry per effect (with the respective
#' name). The entry contains a named vector with the unique stings of the
#' respective effect values i.e. the entries of the respective columns (e.g.
#' name, time, condition etc. for 'biological')
#' @param passParList from the the \code{passParList} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{effectsPars}}{
#'      Named list with one entry per effect. The values are the respective
#'      effect parameters.
#'  }
#'  \item{\code{modelDerivExpr}}{
#'      The jacobian of the argument of the \code{model} parameter of
#'      \link{alignReplicates} parsed as an executable expression.
#'  }
#' }
#'
#' @param passParList2 from the the \code{passParList2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{dataFit}}{The current dataset.}
#' }
#'
#' @noRd
#'

evaluateModelJacobian <- function(initPars,
                                    parList,
                                    passParList,
                                    passParList2) {
    effectsPars <- passParList$effectsPars
    modelDerivExpr <- passParList$modelDerivExpr

    dataFit <- passParList2$dataFit


    biological <- seq_along(initPars)
    useList <- c(
        list(initPars, biological = biological), as.list(dataFit),
        parList
    )
    names(useList)[1] <- effectsPars[[1]][1]
    derivative_values <- with(useList, eval(modelDerivExpr))

    return(derivative_values)
}



# getMaxY() --------------------------------------------------------------

#' Determine the maximal y values for plots
#'
#' A method producing a list of y values that will result in a eye pleasing
#' set of y tics when used in ggplot. Part of \link{plotIt}
#'
#' @param x numerical vector containing the y data which will be used to
#' determine the set of y max
#'
#' @return numerical vector with the maximal y values.
#'
#' @noRd
getMaxY <- function(x) {
    rx <- round(x, 1)
    lower <- floor(x)
    upper <- ceiling(x)
    lowDiff <- rx - lower
    uppDiff <- upper - rx

    if (x > 5) {
        if (upper %% 2 == 0) {
            out <- upper
        } else {
            out <- lower + 1
        }
    } else {
        if (x < 2) {
            if (rx > x) {
                if (rx %% 0.2 == 0) {
                    out <- rx
                } else {
                    out <- rx + 0.1
                }
            } else {
                if (rx %% 0.2 != 0) {
                    out <- rx + 0.1
                } else {
                    out <- rx + 0.2
                }
            }
        } else {
            if (lowDiff < uppDiff) {
                out <- lower + 0.5
            } else {
                out <- upper
            }
        }
    }

    return(out)
}


# scaleValues() ----------------------------------------------------------

#' scale the values according to the current parameter scale
#'
#' @param fitLogscale the current scale
#' @param value value to be scaled
#' @param upper upper error bound, i.e. value + sigma
#' @param lower lower error bound, i.e. value - sigma
#' @param sigma sigma to be scaled
#'
#' @noRd
scaleValues <- function(fitLogscale,
                         outputScale,
                         value,
                         upper,
                         lower,
                         sigma) {
    if (fitLogscale == TRUE) {
        value <- exp(value)
        upper <- exp(upper)
        lower <- exp(lower)
        sigma <- value * sigma
    } else {
        value <- value
        upper <- upper
        lower <- lower
        sigma <- sigma
    }

    upper <- value + sigma
    lower <- value - sigma





    return(list(value = value, upper = upper, lower = lower, sigma = sigma))
}



# getParClass() -------------------------------------------------------

#' get parameter class from from individualized name
#'
#' @param paramlist list of parameters with names of form 'p_i' with parameter class p and
#' index i. e.g. 'yi_1'
#'
#' @return vector with classes e.g. yi
#'
#' @noRd
getParClass <- function(paramlist) {
    classes <- gsub('[_0-9]+', '', names(paramlist))
}
