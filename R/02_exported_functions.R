# readWide() -------------------------------------------------------------

#' Import measurement data
#'
#' Data from a given wide style csv-file is imported. While importing, the data
#' is converted into a long table
#'
#' @param file character, the name of the data file.
#' @param description numeric index vector of the columns containing
#' the description.
#' @param time numeric index of length 1: the time column.
#' @param header a logical value indicating whether the file contains
#' the names of the variables as its first line.
#' @param ... further arguments being passed to \link{read.csv}.
#'
#' @return data frame with columns "name", "time", "value" and other
#' columns describing the measurements.
#'
#' @export
#' @importFrom utils read.csv
#'
#' @examples
#' ## Import example data set
#' simDataWideFile <- system.file(
#'     "extdata", "simDataWide.csv",
#'     package = "blotIt"
#' )
#' readWide(simDataWideFile, description = seq_len(3))


readWide <- function(file, description = NULL, time = 1, header = TRUE, ...) {
    useData <- read.csv(file, header = header, ...)
    allNames <- colnames(useData)

    if (is.null(description)) {
        stop("Specify columns containing descriptions.")
    }
    if (is.character(description)) {
        noDescription <- which(colnames(useData) %in% description)
    } else {
        noDescription <- description
    }
    ## Check the index of the "time" column
    if (is.character(time)) {
        noTime <- which(colnames(useData) == time)
    } else {
        noTime <- time
    }

    ## Check availability of description and time
    if (length(noDescription) < length(description)) {
        warning(
            "Not all columns proposed by argument 'description' are available",
            " in file.\nTaking the available ones."
        )
    }

    if (length(noTime) == 0) {
        stop(
            "File did not contain a time column as proposed by 'time' argument."
        )
    }

    ## biological description data from measurement data
    descrEntries <- useData[, noDescription]
    restLong <- unlist(useData[, -noDescription])

    ## Create output data frame
    newData <- data.frame(
        descrEntries,
        name = rep(allNames[-noDescription], each = dim(useData)[1]),
        value = restLong
    )

    ## Remove missing items
    newData <- newData[!is.nan(newData$value), ]
    newData <- newData[!is.na(newData$value), ]


    colnames(newData)[noTime] <- "time"

    return(newData)
}



# plotIt() ---------------------------------------------------------

#' All-in-one plot function for blotIt
#'
#' Takes the output of \link{alignReplicates} and generates graphs. Which data will be
#' plotted can then be specified separately.
#'
#' @param inputList \code{inputList} file as produced by \link{alignReplicates}, a list
#' of data.frames.
#' @param plotPoints String to specify which data set should be plotted in form
#' of points with corresponding error bars. It must be one of
#' \code{c("original", "scaled", "prediction", "aligned")}.
#' @param plotLine Same as above but with a line and error band.
#' @param spline Logical, if set to \code{TRUE}, what is specified as
#' \code{plotLine} will be plotted as a smooth spline instead of straight lines
#' between points.
#' @param scales String passed as \code{scales} argument to \link{facet_wrap}.
#'
#' @param alignZeros Logical, if \code{TRUE}, the zero ticks are aligned
#' between the facets.
#' @param plotCaption Logical, if \code{TRUE}, a caption describing the plotted
#' data is added to the plot.
#' @param ncol Numerical passed as \code{ncol} argument to
#' \link{facet_wrap}.
#'
#' @param useColors vector of custom color values as taken by the \code{values}
#' argument in the \link{scale_color_manual} method for \code{ggplot} objects.
#'
# @param duplicateZeroPoints Logical, if set to \code{TRUE} all zero time
# points are assumed to belong to the first condition. E.g. when the different
# conditions consist of treatments added at time zero. Default is \code{FALSE}.
#'
#' @param useOrder Optional list of target names in the custom order that will
#' be used for faceting.
#'
#' @param plotScaleX character, defining the scale of the x axis
#'
#' @param plotScaleY character, defining the scale of the y axis
#'
#' @param doseResponse Logical, indicates if the plot should be dose response
#'
#' @param labelX Optional, value passed to \code{xlab} parameter of \link{ggplot}
#' for the x-axis. Default is \code{NULL} leading to 'Time' or 'Dose',
#' respectively.
#'
#' @param labelY Optional, value passed to \code{ylab} parameter of \link{ggplot}
#' for the y-axis. Default is \code{NULL} leading to 'Signal'.
#'
#' @param ... Logical expression used for subsetting the data frames, e.g.
#' \code{name == "pERK1" & time < 60}.
#'
#' \describe{
#' To reproduce the known function \code{plot1}, \code{plot2} and \code{plot3},
#' use:
#' \item{plot1}{
#' \code{plotPoints} = 'original', \code{plotLine} = 'prediction'
#' }
#' \item{plot2}{
#' \code{plotPoints} = 'scaled', \code{plotLine} = 'aligned'
#' }
#' \item{plot3}{
#' \code{plotPoints} = 'aligned', \code{plotLine} = 'aligned'
#' }
#' }
#'
#' @import ggplot2 data.table
#'
#' @return ggplot object
#'
#' @export
#' @author Severin Bang and Svenja Kemmer

plotIt <- function(
  inputList,
  ...,
  plotPoints = "aligned",
  plotLine = "aligned",
  spline = FALSE,
  scales = "free",
  alignZeros = TRUE,
  plotCaption = TRUE,
  ncol = NULL,
  useColors = NULL,
  #duplicateZeroPoints = FALSE,
  useOrder = NULL,
  plotScaleY = NULL,
  plotScaleX = NULL,
  doseResponse = FALSE,
  labelX = NULL,
  labelY = NULL
) {
    if (!plotPoints %in% c("original", "scaled", "prediction", "aligned") |
        !plotLine %in% c("original", "scaled", "prediction", "aligned")) {
        stop(
            "\n\t'plotPoints' and 'plotLine' must each be one of
            c('original', 'scaled', 'prediction', 'aligned')\n"
        )
    }

    ## change plotting order from default
    if (!is.null(useOrder)) {
        if (length(setdiff(levels(inputList[[1]]$name), useOrder)) != 0) {
            stop("useOrder doesn't contain all protein names.")
        } else {
            inputList$aligned$name <- factor(
                inputList$aligned$name,
                levels = useOrder
            )
            inputList$scaled$name <- factor(
                inputList$scaled$name,
                levels = useOrder
            )
            inputList$prediction$name <- factor(
                inputList$prediction$name,
                levels = useOrder
            )
            inputList$original$name <- factor(
                inputList$original$name,
                levels = useOrder
            )
        }
    }

    biological <- inputList$biological
    scaling <- inputList$scaling

    plotList <- inputList

    ## duplicate 0 values for all doses # not working at the moment!
    # if (duplicateZeroPoints) {
    #     for (ndat in 1) {
    #         dat <- plotList[[ndat]]
    #         subsetZeros <- copy(subset(dat, time == 0))
    #         myDoses <- setdiff(unique(dat$dose), 0)
    #         myZerosAdd <- NULL
    #         for (d in seq(1, length(myDoses))) {
    #             subsetZerosD <- copy(subsetZeros)
    #             subsetZerosD$dose <- myDoses[d]
    #             myZerosAdd <- rbind(myZerosAdd, subsetZerosD)
    #         }
    #         dat <- rbind(dat, myZerosAdd)
    #         plotList[[ndat]] <- dat
    #     }
    # }

    ## add columns containing the respective scaling and biological effects

    # aligned
    if (doseResponse == TRUE) {
        xValue <- "dose"
    } else {
        xValue <- "time"
    }


    plotList$aligned$biological <- do.call(
        paste0,
        plotList$aligned[
            ,
            biological[!(biological %in% c("name", xValue))],
            drop = FALSE
        ]
    )
    plotList$aligned$scaling <- NA

    ## scaled
    plotList$scaled$biological <- do.call(
        paste0,
        inputList$scaled[
            ,
            biological[!(biological %in% c("name", xValue))],
            drop = FALSE
        ]
    )
    plotList$scaled$scaling <- do.call(
        paste0,
        inputList$scaled[, scaling[scaling != "name"], drop = FALSE]
    )

    ## prediction
    plotList$prediction$biological <- do.call(
        paste0,
        inputList$prediction[
            ,
            biological[!(biological %in% c("name", xValue))],
            drop = FALSE
        ]
    )
    plotList$prediction$scaling <- do.call(
        paste0,
        inputList$prediction[, scaling[scaling != "name"], drop = FALSE]
    )

    ## original
    plotList$original$biological <- do.call(
        paste0,
        inputList$original[
            ,
            biological[!(biological %in% c("name", xValue))],
            drop = FALSE
        ]
    )
    plotList$original$scaling <- do.call(
        paste0,
        inputList$original[, scaling[scaling != "name"], drop = FALSE]
    )

    plotListPoints <- plotList[[plotPoints]]
    plotListLine <- plotList[[plotLine]]

    plotListPoints <- subset(plotListPoints, ...)
    plotListLine <- subset(plotListLine, ...)

    ## build Caption
    usedErrors <- list(
        aligned = "Fisher Information",
        scaled = "Propergated error model to common scale",
        prediction = "Error model",
        original = "None"
    )

    usedData <- list(
        aligned = "Estimated true values",
        scaled = "Original data scaled to common scale",
        prediction = "Predictions from model evaluation on original scale",
        original = "Original data"
    )

    captionText <- paste0(
        "Datapoints: ", usedData[[plotPoints]], "\n",
        "Errorbars: ", usedErrors[[plotPoints]], "\n",
        "Line: ", usedData[[plotLine]], "\n",
        if (plotPoints != plotLine) {
            paste0("Errorband: ", usedErrors[[plotLine]], "\n")
        },
        "\n",
        "Date: ", Sys.Date()
    )

    ## we want to keep the x ticks!
    if (scales == "biological") {
        scales <- "free_x"
    }



# * settings for dose response/time course --------------------------------

    if (doseResponse == TRUE) {
        Xlabel <-  "Dose"
        xVariable <- "dose"
    } else {
        Xlabel <-  "Time"
        xVariable <- "time"
    }

    Ylabel <- "Signal"

    if (!is.null(labelX)){
        Xlabel <- labelX
    }

    if (!is.null(labelY)){
        Ylabel <- labelY
    }

    if (!is.null(plotScaleX)) {
        errWidth <- 0
    } else {
        errWidth <- max(plotListPoints[xVariable])/50
    }


# * plotting --------------------------------------------------------------
    if (plotPoints == "aligned" & plotLine == "aligned") {
        g <- ggplot(
            data = plotListPoints,
            aes_string(
                x = xVariable,
                y = "value",
                group = "biological",
                color = "biological",
                fill = "biological"
            )
        )
        g <- g + facet_wrap(~name, scales = scales, ncol = ncol)
        if (is.null(useColors)) {
            # useColors <- scale_color_brewer()
            # # c(
            # # "#000000",
            # # "#C5000B",
            # # "#0084D1",
            # # "#579D1C",
            # # "#FF950E",
            # # "#4B1F6F",
            # # "#CC79A7",
            # # "#006400",
            # # "#F0E442",
            # # "#8B4513",
            # # rep("gray", 100)
            # # )
            # g <- g + scale_color_manual("Condition", values = useColors) +
            #     scale_fill_manual("Condition", values = useColors)
        } else {
            useColors <- c(useColors, rep("gray", 100))
            g <- g + scale_color_manual("Biological", values = useColors) +
                scale_fill_manual("Biological", values = useColors)
        }
    } else {
        g <- ggplot(
            data = plotListPoints,
            aes_string(
                x = xVariable,
                y = "value",
                group = "scaling",
                color = "scaling",
                fill = "scaling"
            )
        )
        g <- g + facet_wrap(
            ~ name * biological,
            scales = scales,
            ncol = ncol
        )
        if (is.null(useColors)) {
            # useColors <- scale_color_brewer()
            # # c(
            # # "#000000",
            # # "#C5000B",
            # # "#0084D1",
            # # "#579D1C",
            # # "#FF950E",
            # # "#4B1F6F",
            # # "#CC79A7",
            # # "#006400",
            # # "#F0E442",
            # # "#8B4513",
            # # rep("gray", 100)
            # # )
            # g <- g + scale_color_manual("Scaling", values = useColors) +
            #     scale_fill_manual("Scaling", values = useColors)
        } else {
            useColors <- c(useColors, rep("gray", 100))
            g <- g + scale_color_manual("Scaled", values = useColors) +
                scale_fill_manual("Scaled", values = useColors)
        }
    }


    g <- g + geom_point(data = plotListPoints, size = 2.5)
    g <- g + geom_line(data = plotListLine, size = 1)


# * * plotPoints == original ----------------------------------------------
    if (plotPoints == "original") {
        if (plotLine != "original") {
            if (spline == TRUE) {
                g <- g + geom_smooth(
                    data = plotListLine,
                    se = FALSE,
                    method = "lm",
                    formula = y ~ poly(x, 3)
                )
            } else {
                g <- g + geom_ribbon(
                    data = plotListLine,
                    aes(
                        ymin = lower, # value - sigma,
                        ymax = upper, # value + sigma#,
                        # fill = "grey",
                        # color = "grey"
                    ),
                    alpha = 0.3,
                    lty = 0
                )
            }
        }


# * * plotPoints == prediction -------------------------------------------
    } else if (plotPoints == "prediction") {
        if (plotLine != "prediction") {
            g <- g + geom_errorbar(
                data = plotListPoints,
                aes(
                    ymin = lower, # value - sigma,
                    ymax = upper # value + sigma
                ),
                size = 0.5,
                width = errWidth,
                alpha = 0.5
            )
        }
        if (spline == TRUE) {
            g <- g + geom_smooth(
                data = plotListLine,
                se = FALSE,
                method = "lm",
                formula = y ~ poly(x, 3)
            )
        } else {
            g <- g + geom_ribbon(
                data = plotListLine,
                aes(
                    ymin = lower, # value - sigma,
                    ymax = upper, # value + sigma#,
                    # fill = "grey",
                    # color = "grey"
                ),
                alpha = 0.3,
                lty = 0
            )
        }


# * * plotPoints == scaled ------------------------------------------------
    } else if (plotPoints == "scaled") {
        if (plotLine != "scaled") {
            g <- g + geom_errorbar(
                data = plotListPoints,
                aes(
                    ymin = lower, # value - sigma,
                    ymax = upper # value + sigma
                ),
                size = 0.5,
                width = errWidth,
                alpha = 0.5
            )
        }
        if (spline == TRUE) {
            g <- g + geom_smooth(
                data = plotListLine,
                se = FALSE,
                method = "lm",
                formula = y ~ poly(x, 3)
            )
        } else {
            g <- g + geom_ribbon(
                data = plotListLine,
                aes(
                    ymin = lower, # value - sigma,
                    ymax = upper, # value + sigma#,
                    # fill = "grey",
                    # color = "grey"
                ),
                alpha = 0.3,
                lty = 0
            )
        }

#  * * plotPoints == aligned ----------------------------------------------
    } else if (plotPoints == "aligned") {
        if (plotLine != "aligned") {
            g <- g + geom_errorbar(
                data = plotListPoints,
                aes(
                    ymin = lower, # value - sigma,
                    ymax = upper # value + sigma
                ),
                size = 0.5,
                width = errWidth,
                alpha = 0.5
            )
        } else {
            if (spline == TRUE) {
                g <- g + geom_smooth(
                    data = plotListLine,
                    se = FALSE,
                    method = "lm",
                    formula = y ~ poly(x, 3)
                )
            } else {
                g <- g + geom_ribbon(
                    data = plotListLine,
                    aes(
                        ymin = lower, # value - sigma,
                        ymax = upper, # value + sigma#,
                        # fill = "grey",
                        # color = "grey"
                    ),
                    alpha = 0.3,
                    lty = 0
                )
            }
        }
    }





    g <- g + theme_bw(base_size = 20) +
        theme(
            legend.position = "top",
            legend.key = element_blank(),
            strip.background = element_rect(color = NA, fill = NA),
            axis.line.x = element_line(size = 0.3, colour = "black"),
            axis.line.y = element_line(size = 0.3, colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
        )
    g <- g + xlab(paste0("\n",Xlabel)) + ylab(paste0(Ylabel,"\n"))

    if (alignZeros) {
        if (plotPoints != "original") {
            ## scale y-axes (let them start at same minimum determined by
            ## smallest value-sigma and end at individual ymax)
            plotListPoints <- as.data.table(plotListPoints)
            blankData <- plotListPoints[
                ,
                list(ymax = max(upper), ymin = min(lower)),
                by = c("name", "biological", "scaling")
            ]
            blankData[, ":="(ymin = min(ymin))] # same minimum for all proteins
            blankData[
                ,
                ":="(ymax = getMaxY(ymax)),
                by = c("name", "biological", "scaling")
            ] # protein specific maximum
            blankData <- melt(
                blankData,
                id.vars = c("name", "biological", "scaling"),
                measure.vars = c("ymax", "ymin"),
                value.name = "value"
            )
            blankData[, ":="(xVariable = 0, variable = NULL)]
            setnames(blankData, "xVariable", xVariable)
            g <- g + geom_blank(
                data = as.data.frame(blankData),
                aes_string(x = xVariable, y = "value")
            )
        }
    }


    if (plotCaption) {
        g <- g + labs(caption = captionText)
    }

    if (is.null(plotScaleY)) {
        if (inputList$outputScale != "linear") {
            g <- g + coord_trans(y = inputList$outputScale)
        }
    } else if (plotScaleY %in% c("log", "log2", "log10")) {
        g <- g + scale_y_continuous(trans = plotScaleY)
    }

    if (is.null(plotScaleX)) {
        if (inputList$outputScale != "linear") {
            g <- g + coord_trans(x = inputList$outputScale)
        }
    } else if (plotScaleX == "pseudoLog10") {
      g <- g + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))
    } else if (plotScaleX %in% c("log", "log2", "log10")) {
        g <- g + scale_x_continuous(trans = plotScaleX)
    }


    return(g)
}


# llrTest() --------------------------------------------------------------

#' Method for hypothesis testing
#'
#' Two outputs of \link{alignReplicates} can be tested as nested hypothesis. This can
#' be used to test if e.g. buffer material influence can be neglected or a
#' specific measurement point is an outlier.
#'
#' @param H0 output of \link{alignReplicates} obeying the null hypothesis. A special
#' case of \code{H1}.
#' @param H1 output of \link{alignReplicates}, the general case
#'
#' @return list with the log-likelihood ratio, statistical information and the
#' numerical p-value calculated by the evaluating the chi-squared distribution
#' at the present log-likelihood ratio with the current degrees of freedom.
#'
#' @examples
#' ## load provided example data file
#' lrrDataPath <- system.file(
#'     "extdata", "exampleLlrTest.csv",
#'     package = "blotIt"
#' )
#'
#' ## import data
#' llrData <- readWide(
#'     file = lrrDataPath,
#'     description = seq(1, 4),
#'     sep = ",",
#'     dec = "."
#' )
#'
#' ## generate H0: the buffer column is not named as a biological effect e.g.
#' ## not considered as a biological different condition
#' H0 <- alignReplicates(
#'     data = llrData,
#'     model = "yi / sj",
#'     errorModel = "value * sigmaR",
#'     biological = yi ~ name + time + stimmulus,
#'     scaling = sj ~ name + ID,
#'     error = sigmaR ~ name + 1,
#'     fitLogscale = FALSE,
#'     normalize = TRUE,
#'     averageTechRep = FALSE,
#'     verbose = FALSE,
#'     normalizeInput = TRUE
#' )
#'
#' ## generate H1: here the buffer column is named in the biological parameter
#' ## therefore different entries are considered as biologically different
#' H1 <- alignReplicates(
#'     data = llrData,
#'     model = "yi / sj",
#'     errorModel = "value * sigmaR",
#'     biological = yi ~ name + time + stimmulus + buffer,
#'     scaling = sj ~ name + ID,
#'     error = sigmaR ~ name + 1,
#'     fitLogscale = FALSE,
#'     normalize = TRUE,
#'     averageTechRep = FALSE,
#'     verbose = FALSE,
#'     normalizeInput = TRUE
#' )
#'
#' ## perform test
#' llrTest(H0, H1)
#' @export
#'
llrTest <- function(H0, H1, check = TRUE) {
    biological0 <- union(
        H0$biological[!(H0$biological %in% c("name", "time"))],
        "1"
    )
    biological1 <- union(
        H1$biological[!(H1$biological %in% c("name", "time"))],
        "1"
    )

    scaling0 <- union(H0$scaling[H0$scaling != "name"], "1")
    scaling1 <- union(H1$scaling[H1$scaling != "name"], "1")

    error0 <- union(H0$error[H0$error != "name"], "1")
    error1 <- union(H1$scaling[H1$scaling != "name"], "1")

    if (check) {
        if (
            !all(biological0 %in% biological1) |
                !all(scaling0 %in% scaling1) | !all(error0 %in% error1)
        ) {
            stop("H0 is not a special case of H1.")
        }
    }


    value0 <- attr(H0$parameter, "value")
    value1 <- attr(H1$parameter, "value")
    df0 <- attr(H0$parameter, "df")
    df1 <- attr(H1$parameter, "df")

    list(
        llr = value0 - value1,
        statistic = paste0(
            "chisquare with ", df0 - df1, " degrees of freedom."
        ),
        p.value = pchisq(value0 - value1, df = df0 - df1, lower.tail = FALSE)
    )
}
