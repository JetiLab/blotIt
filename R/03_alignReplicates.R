
#' Align time-course data based on an Mixed-Effects alignment model
#'
#' The function deals primarily with time-course data of different
#' targets which have been measured under different experimental
#' conditions and whose measured values might be on a different
#' scale, e.g. because of different amplification. The algorithm
#' determines the different scaling and estimates the time-course on
#' a common scale.
#'
#'
#' @param data data.frame containing the the data to be scaled. Usualy the
#' output of \link{readWide} is used. Obligatory are the columns \code{name},
#' and \code{value}. In the case of time dependet data, a \code{time} column is
#' also necessary. Additionally, \code{data} should have further columns,
#' e.g. characterizing experimental conditions (the fixed effects) and
#' sources of variance between data of the same experimental condition
#' (the scaled variables).
#'
#' @param model character defining the model by which the values in
#' \code{data} can be described, e.g. "yi/sj"
#'
#' @param errorModel character defining a model for the standard
#' deviation of a value, e.g. "sigma0 + value * sigmaR". This model
#' can contain parameters, e.g. "sigma0" and "sigmaR", or numeric
#' variables from \code{data}, e.g. "value" or "time".
#'
#' @param biological two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameters contained in \code{model}, e.g. "yi", and "name1, ..."
#' refers to variables in \code{data}, e.g. "condition".
#' The parameters "par1, ..." are determined specific to the levels of
#' "name1, ...".
#'
#' @param scaling two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameters contained in \code{model}, e.g. "sj", and "name1, ..."
#' refers to variables in \code{data}, e.g. "Experiment".
#' @param error two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameter contained in \code{error}, e.g. "sigma0" and "sigmaR", and
#' "name1, ..." refers to variables in \code{data}. If the same values
#' of "par1, ..." should be assumed for all data, the right hand side of ~ has
#' to reffer to the "name" column of \code{data}.
#' @param fitLogscale logical, defining if the parameters are
#' fitted on a log scale. Computational reasons.
#' @param normalize logical indicating whether the biological effect
#' parameter should be normalized to unit mean.
#'
#' @param averageTechRep logical, indicates if the technical replicates
#' should be averaged
#'
#' @param namesFactored logical, indicates if the \code{name} columns of the
#' output should be parsed as factors instead of characters.
#'
#' @param verbose logical, print out information about each fit
#' @param normalizeInput logical, if TRUE the input will be normalized before
#' scaling, helpful if convergence fails because the data varies for to many
#' orders of magnitude.
#' @param iterlim numerical argument passed to \link{trust}.
#'
#' @param ciProfiles Logical, if \code{TRUE}, the confidence intervals (CI) are
#' calculated based on the profile likelihood (PL). Default is \code{FALSE}, the
#' CI will then be approximated by the Fisher information (FI). Calculation via
#' PL is more exact but roughly 20times slower. The FI approximates the CI
#' conservatively compared to the PL (FI leads to larger CI).
#'
#' @details Alignment of time-course data is achieved by an alignment
#' model which explains the observed data by a function mixing
#' fixed effects, usually parameters reflecting the "underlying"
#' time-course, and latent variables, e.g. scaling parameters taking
#' account for effects like different amplification or loading, etc.
#' Depending on the measurement technique, the data has constant
#' relative error, or constant absolute error or even a combination
#' of those. This error is described by an error function. The error
#' parameters are usually global, i.e. the same parameter values are
#' assumed for all data points.
#'
#' @return Object of class \code{aligned}, i.e. a data frame of the
#' alignment result containing an attribute "outputs":
#'  a list of data frames
#' \describe{
#' \item{aligned}{data.frame with the original column names plus the column
#'     \code{sigma}. Each set of unique biological effects i.e. biological
#'     different condition (e.g. time point, target and treatment) has one set
#'     of \code{value} and \code{sigma}. The values are the estimated true
#'     values i.e. the determined biological parameters. The errors in the
#'     \code{sigma} column are estimated by employing the fisher information
#'     to quantify the uncertainty of the respective fit. Both, the value and
#'     its error are on the common scale.}
#' \item{scaled}{The original measurements scaled to common scale by applying
#'     inverse model. The errors are the result of the evaluation of the error
#'     model and then also scaled to common scale by use of Gaussian error
#'     propagation.}
#' \item{prediction}{Original data with \code{value} replaced by the prediction
#'     (evaluation of the model with the fitted parameters), and \code{sigma}
#'     from the evaluation of the error model. Both are on the original scale.
#'     }
#'
#' \item{original}{The original data as passed as  \code{data}.}
#' \item{originalWithParameters}{The original data but with added columns
#'     containing the estimated parameters}
#'
#' \item{parameter}{Parameter table with the columns: \code{name}: the name of
#'     the current target, \code{level}: the pasted unique set of effects
#'     (biological, scaling or error), \code{parameter}: the parameter
#'     identifier as defined in the (error) model, \code{value} and \code{sigma}
#'     containing the determined values and corresponding errors, \code{nll}:
#'     twice the negative log-likelihood of the fit, \code{noPars} and
#'     \code{noData} containing the number of parameters and data points for
#'     the respected fit. This list entry also has two attributes: \code{value}
#'     containing the final value (residual sum of squares) passed to
#'     \link{trust::trust} by the objective function and \code{df} the degrees
#'     of freedom of the fitting process.}
#'
#' \item{biological}{Names of the columns containing the biological effects}
#' \item{scaling}{Names of the columns containing the scaling effects}
#' \item{error}{Names of the columns containing the error effects}
#' }
#'
#' The estimated parameters are returned by the attribute "parameters".
#' @example inst/examples/exampleAlignReplicates.R
#' @seealso \link{readWide} to read data in a wide column format and
#' get it in the right format for \code{alignReplicates}.
#'
#' @importFrom stats D
#'
#' @export
alignReplicates <- function(data,
                     model = NULL,
                     errorModel = NULL,
                     biological = NULL,
                     scaling = NULL,
                     error = NULL,
                     fitLogscale = TRUE,
                     normalize = TRUE,
                     averageTechRep = FALSE,
                     verbose = FALSE,
                     namesFactored = TRUE,
                     normalizeInput = TRUE,
                     outputScale = "linear",
                     iterlim = 100
                     ) {

    ## check input for mistakes
    inputCheckReport <- inputCheck(
        data = data,
        model = model,
        errorModel = errorModel,
        biological = biological,
        scaling = scaling,
        error = error,
        fitLogscale = fitLogscale,
        outputScale = outputScale
    )

    if (verbose) {
        cat(inputCheckReport, "\n")
    }

    ## Check if data is already blotIt output
    parameterData <- NULL

    if (inherits(data, "aligned")) {
        parameterData <- data$parameters
        data <- data$original
    }

    ## read biological, scaling and error effects from input
    effects <- identifyEffects(
        biological = biological,
        scaling = scaling,
        error = error
    )

    effectsValues <- effects$effectsValues
    effectsPars <- effects$effectsPars

    ## prepare data
    data <- as.data.frame(data)

    toBeScaled <- splitForScaling(
        data,
        effectsValues,
        normalizeInput
    )

    ## Generate unique list of targets
    targets <- make.unique(
        vapply(toBeScaled, function(d) as.character(d$name)[1], ""),
        sep = "_"
    )

    ## Rename names to unique if multiple scalings per target exist
    for (n in seq_along(targets)) {
        toBeScaled[[n]]$name <- targets[n]
    }

    ## Get biological, scaling and error parameter from model and error model
    parameters <- getSymbols(c(model, errorModel), exclude = colnames(data))

    ## Include the normalization term as a constraint if said so in the function
    ## call
    if (normalize) {
        constraint <- paste("1e3 * (mean(", effectsPars[1][1], ") - 1)")
        constStrength <- 1000
    } else {
        constraint <- "0"
        constStrength <- 0
    }

    ## Retrieve the covariates as the "remaining" model parameters, when fixed,
    ## latent and errorparameters are excluded
    covariates <- union(
        getSymbols(
            model,
            exclude = c(effectsPars[1], effectsPars[2], effectsPars[3])
        ),
        getSymbols(
            errorModel,
            exclude = c(effectsPars[1], effectsPars[2], effectsPars[3])
        )
    )
    cat("Covariates:", paste(covariates, sep = ", "), "\n")

    ## Check if parameters from (error-) model and passed expressions coincide
    if (
        length(
            setdiff(
                unlist(c(
                    effectsPars[1],
                    effectsPars[2],
                    effectsPars[3]
                )),
                parameters
            )
        ) > 0
    ) {
        stop("Not all paramters are defined in either arguments
         'scaling', 'biological' or 'error'")
    }

    ## Name the respective parameters fixed, latent and error
    names(parameters)[parameters %in% effectsPars[1]] <- "biological"
    names(parameters)[parameters %in% effectsPars[2]] <- "scaling"
    names(parameters)[parameters %in% effectsPars[3]] <- "error"

    ## parse error model by replacing the "value" by the model
    errorModel <- replaceSymbols(
        "value",
        paste0("(", model, ")"),
        errorModel
    )

    ## Apply transformation according to given 'fitLogscale' parameter
    if (fitLogscale == TRUE) {
        model <- replaceSymbols(parameters, paste0(
            "exp(", parameters,
            ")"
        ), model)
        errorModel <- replaceSymbols(parameters, paste0(
            "exp(",
            parameters, ")"
        ), errorModel)
        constraint <- replaceSymbols(parameters, paste0(
            "exp(",
            parameters, ")"
        ), constraint)
        if (verbose) {
            cat("model, errormodel and constraint are scaled: x <- exp(x)\n")
        }
    } else if (verbose == TRUE) {
        cat("model, errormodel and constraint remain linear scaled.\n")
    }


    ## Calculating the derivative
    modelDeriv <- deparse(
        stats::D(parse(text = model), name = effectsPars[[1]][1])
    )

    ## Calculating the (error) model jacobian
    modelJacobian <- lapply(
        parameters,
        function(p) {
            deparse(stats::D(parse(text = model), name = p))
        }
    )
    errorModelJacobian <- lapply(
        parameters,
        function(p) {
            deparse(stats::D(parse(text = errorModel), name = p))
        }
    )

    ## Construct math function expressions
    constrExpr <- parse(text = constraint)

    ## Replace the parameters used in the (error) model and the derivatives
    ## by placeholders of the respective effect category
    for (n in seq_along(parameters)) {
        model <- replaceSymbols(parameters[n], paste0(
            parameters[n],
            "[", names(parameters)[n], "]"
        ), model)
        errorModel <- replaceSymbols(parameters[n], paste0(
            parameters[n],
            "[", names(parameters)[n], "]"
        ), errorModel)
        modelDeriv <- replaceSymbols(parameters[n], paste0(
            parameters[n],
            "[", names(parameters)[n], "]"
        ), modelDeriv)
        modelJacobian <- lapply(modelJacobian, function(myjac) {
            replaceSymbols(
                parameters[n],
                paste0(
                    parameters[n], "[", names(parameters)[n],
                    "]"
                ), myjac
            )
        })
        errorModelJacobian <- lapply(errorModelJacobian, function(myjac) {
            replaceSymbols(
                parameters[n],
                paste0(
                    parameters[n], "[", names(parameters)[n],
                    "]"
                ), myjac
            )
        })
    }

    cat("Model:        ", model, "\n", sep = "")
    cat("Error Model:  ", errorModel, "\n", sep = "")

    ## Parse functions to an executable expression
    modelExpr <- parse(text = model)
    modelDerivExpr <- parse(text = modelDeriv)
    errorModelExpr <- parse(text = errorModel)

    modelJacobianExpr <- lapply(
        modelJacobian,
        function(myjac) parse(text = myjac)
    )
    errorModelJacobianExpr <- lapply(
        errorModelJacobian,
        function(myjac) parse(text = myjac)
    )


    passParList <- list(
        effectsValues = effectsValues,
        parameterData = parameterData,
        averageTechRep = averageTechRep,
        verbose = verbose,
        covariates = covariates,
        parameters = parameters,
        fitLogscale = fitLogscale,
        effectsPars = effectsPars,
        modelExpr = modelExpr,
        errorModelExpr = errorModelExpr,
        constrExpr = constrExpr,
        modelJacobianExpr = modelJacobianExpr,
        errorModelJacobianExpr = errorModelJacobianExpr,
        constStrength = constStrength,
        normalize = normalize,
        modelDerivExpr = modelDerivExpr,
        outputScale = outputScale,
        iterlim = iterlim
    )

# * Hand over to scaling --------------------------------------------------
    out <- lapply(
        seq_along(toBeScaled),
        function(i,
                 passParList) {
            try({
                cat(
                    "Target ", i, "/", length(toBeScaled), ":",
                    targets[i], "\n", sep = ""
                    )
                out <- scaleTarget(
                    currentData = toBeScaled[[i]],
                    passParList = passParList
                )

                return(out)
            },
            silent = FALSE
            )
        },
        ## additional parameters for scaleTarget
        passParList = passParList
    )

    ## Sanitize results from failed fits
    ## rbind parameter table from not-failed elements of out
    parTable <- do.call(
        rbind,
        lapply(
            out,
            function(o) {
                if (!inherits(o, "try-error")) {
                    o[[6]]
                } else {
                    NULL
                }
            }
        )
    )

    ## add up the "values", i.e. the -2*LL
    attr(parTable, "value") <- do.call(
        sum,
        lapply(
            out,
            function(o) {
                if (!inherits(o, "try-error")) {
                    attr(o[[6]], "value")
                } else {
                    0
                }
            }
        )
    )

    ## add up degrees of freedom
    attr(parTable, "df") <- do.call(
        sum,
        lapply(
            out,
            function(o) {
                if (!inherits(o, "try-error")) {
                    attr(o[[6]], "df")
                } else {
                    0
                }
            }
        )
    )

    ## paste together the other results
    outCombined <- lapply(
        seq_len(5),
        function(i) {
            do.call(
                rbind,
                lapply(
                    out,
                    function(o) {
                        if (!inherits(o, "try-error")) {
                            o[[i]]
                        } else {
                            NULL
                        }
                    }
                )
            )
        }
    )

    names(outCombined) <- c(
        "prediction",
        "scaled",
        "aligned",
        "original",
        "originalWithParameters"
    )

    returnList <- list(
        aligned = outCombined$aligned,
        scaled = outCombined$scaled,
        prediction = outCombined$prediction,
        original = outCombined$original,
        originalWithParameters = outCombined$originalWithParameters,
        parameter = parTable,
        biological = effectsValues[[1]],
        scaling = effectsValues[[2]],
        error = effectsValues[[3]],
        outputScale = outputScale
    )

    if (namesFactored == TRUE) {
        returnList$aligned$name <- factor(
            returnList$aligned$name,
            levels = unique(returnList$aligned$name)
        )
        returnList$scaled$name <- factor(
            returnList$scaled$name,
            levels = unique(returnList$scaled$name)
        )
        returnList$prediction$name <- factor(
            returnList$prediction$name,
            levels = unique(returnList$prediction$name)
        )
        returnList$original$name <- factor(
            returnList$original$name,
            levels = unique(returnList$original$name)
        )
        returnList$originalWithParameters$name <- factor(
            returnList$originalWithParameters$name,
            levels = unique(returnList$originalWithParameters$name)
        )
    }

    return(returnList)
}
