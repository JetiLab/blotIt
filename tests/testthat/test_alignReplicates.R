# alignReplicates() - model parsing ----------------------------------------------

test_that("alignReplicates() - model parsing", {
    simDataWideFile <- system.file(
        "extdata", "simDataWide.csv",
        package = "blotIt"
    )

    input_data <- readWide(simDataWideFile, description = seq_len(3))

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = "one_parameter",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1
        ),
        "Do not pass biological, scaling or error as string."
    )

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel <- "value * sigmaR",
            biological <- yi ~ name + time + condition,
            error = sigmaR ~ name + 1
        ),
        "All of model, errorModel, biological, scaling, error must be set"
    )

    expect_error(
        alignReplicates(
            data = input_data,
            biological = NULL,
            error = left ~ right2
        ),
        "All of model, errorModel, biological, scaling, error must be set."
    )

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = yi ~ condition,
            scaling = wrong ~ ID,
            error = sigmaR ~ name + 1,
            fitLogscale = TRUE
        ),
        "Not all paramters are defined in either arguments
         'scaling', 'biological' or 'error'"
    )

    expect_error(
        alignReplicates(
            data = input_data,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = yi ~ condition,
            scaling = sj ~ ID,
            error = sigmaR ~ name + 1,
            fitLogscale = "string"
        ),
        "'fitLogscale' must be logical"
    )
})
