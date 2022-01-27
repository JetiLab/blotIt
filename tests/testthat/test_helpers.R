
# getSymbols() -----------------------------------------------------------

test_that("getSymbols()", {
    expect_equal(
        getSymbols("left1 ~ right11 + right12"),
        c("left1", "right11", "right12")
    )
    expect_equal(
        getSymbols(
            "left2 ~ (right21 * right22 / (right23 * right24)) - right25"
        ),
        c("left2", "right21", "right22", "right23", "right24", "right25")
    )
})


# splitForScaling() -----------------------------------------------------

test_that("splitForScaling()", {
    simDataWideFile <- system.file(
        "extdata", "simDataWide.csv",
        package = "blotIt"
    )
    sim_data_long <- readWide(simDataWideFile, description = seq_len(3))

    effectsValues <- list(
        biologicaValues = c("name", "time", "condition"),
        scalingValues = c("name", "ID"),
        errorValues = "name"
    )

    effectsValues_1 <- list(
        biologicaValues = c("name", "time", "condition"),
        scalingValues = c("name"),
        errorValues = "name"
    )

    effectsValues_2 <- list(
        biologicaValues = c("name", "time", "condition"),
        scalingValues = c("ID"),
        errorValues = "name"
    )

    expect_equal(
        length(
            splitForScaling(
                data = sim_data_long,
                effectsValues = effectsValues,
                normalizeInput = TRUE
            )
        ),
        14
    )

    expect_equal(
        length(
            splitForScaling(
                data = sim_data_long,
                effectsValues = effectsValues_1,
                normalizeInput = TRUE
            )
        ),
        nrow(unique(sim_data_long["name"]))
    )

    expect_equal(
        length(
            splitForScaling(
                data = sim_data_long,
                effectsValues = effectsValues_2,
                normalizeInput = TRUE
            )
        ),
        1
    )
})



# replaceSymbols() -------------------------------------------------------

test_that("replaceSymbols()", {
    expect_equal(
        replaceSymbols(
            what = "before",
            by = "after",
            x = "this ~ before"
        ),
        "this~after"
    )

    expect_equal(
        replaceSymbols(
            what = c("before1", "before2", "before3"),
            by = c("after1", "after2", "after3"),
            x = "this ~ before1*before2/before3"
        ),
        "this~after1*after2/after3"
    )
})



# analyzeBlocks() --------------------------------------------------------

test_that("analyzeBlocks()", {
    test_block_matrix <- matrix(
        c(
            1, 0, 0, 0,
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 1, 1, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
            0, 0, 0, 1
        ),
        nrow = 7,
        byrow = TRUE
    )
    expect_equal(
        analyzeBlocks(
            test_block_matrix
        ),
        list(
            c(1, 2),
            c(3, 4, 5),
            c(6, 7)
        )
    )
    expect_equal(
        length(
            analyzeBlocks(
                test_block_matrix
            )
        ),
        3
    )
})


# inputCheck() -----------------------------------------------------------

test_that("inputCheck()", {
    simDataWideFile <- system.file(
        "extdata", "simDataWide.csv",
        package = "blotIt"
    )
    sim_data_long <- readWide(simDataWideFile, description = seq_len(3))

    expect_equal(
        inputCheck(
            data = sim_data_long,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = yi ~ name + time + condition,
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            fitLogscale = TRUE
        ),
        "All input checks passed."
    )

    expect_error(
        inputCheck(
            data = sim_data_long,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = yi ~ name + time + condition + wrong,
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            fitLogscale = TRUE
        ),
        paste0(
            "Not all column names set in 'biological', 'scaling' and 'error' ",
            "are present in 'data'."
        )
    )

    expect_error(
        inputCheck(
            data = sim_data_long,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            # biological = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            fitLogscale = TRUE
        ),
        "All of model, errorModel, biological, scaling, error must be set."
    )

    expect_error(
        inputCheck(
            data = sim_data_long,
            # model = "yi / sj",
            errorModel = "value * sigmaR",
            # biological = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            # error = sigmaR ~ name + 1,
            fitLogscale = TRUE
        ),
        "All of model, errorModel, biological, scaling, error must be set."
    )

    expect_error(
        inputCheck(
            data = sim_data_long,
            model = "yi / sj",
            errorModel = "value * sigmaR",
            biological = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            fitLogscale = "wrong"
        ),
        "'fitLogscale' must be logical"
    )
})


# genInitPars() -------------------------------------------------
test_that("genInitPars()", {
    parameters <- c(
        "biological" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    levelsList <- list(
        biological = c(
            "biological_1",
            "biological_2",
            "biological_3",
            "biological_4",
            "biological_5",
            "biological_6",
            "biological_7"
        ),
        scaling = c(
            "scaling_1",
            "scaling_2"
        ),
        error = "error1"
    )


    expect_equal(
        genInitPars(
            parameters = parameters,
            fitLogscale = FALSE,
            levelsList
        ),
        c(
            "yi_1" = 1,
            "yi_2" = 1,
            "yi_3" = 1,
            "yi_4" = 1,
            "yi_5" = 1,
            "yi_6" = 1,
            "yi_7" = 1,
            "sj_8" = 1,
            "sj_9" = 1,
            "sigmaR_10" = 1
        )
    )

    expect_equal(
        genInitPars(
            parameters = parameters,
            fitLogscale = TRUE,
            levelsList
        ),
        c(
            "yi_1" = 0,
            "yi_2" = 0,
            "yi_3" = 0,
            "yi_4" = 0,
            "yi_5" = 0,
            "yi_6" = 0,
            "yi_7" = 0,
            "sj_8" = 0,
            "sj_9" = 0,
            "sigmaR_10" = 0
        )
    )

    expect_equal(
        genInitPars(
            parameters = parameters,
            fitLogscale = TRUE,
            levelsList = list(
                biological = c(
                    "pAKT_0_0Uml Epo",
                    "pAKT_5_0Uml Epo",
                    "pAKT_10_0Uml Epo"
                ),
                scaling = c(
                    "scaling_1",
                    "scaling_2"
                ),
                error = "error1"
            )
        ),
        c(
            "yi_1" = 0,
            "yi_2" = 0,
            "yi_3" = 0,
            "sj_4" = 0,
            "sj_5" = 0,
            "sigmaR_6" = 0
        )
    )

    expect_equal(
        genInitPars(
            parameters = parameters,
            fitLogscale = FALSE,
            levelsList = list(
                biological = c(
                    "biological_a",
                    "biological_b",
                    "biological_c"
                ),
                scaling = c(
                    "scaling_a",
                    "scaling_b",
                    "scaling_c"
                ),
                error = c(
                    "error_a",
                    "error_b",
                    "error_c"
                )
            )
        ),
        c(
            "yi_1" = 1,
            "yi_2" = 1,
            "yi_3" = 1,
            "sj_4" = 1,
            "sj_5" = 1,
            "sj_6" = 1,
            "sigmaR_7" = 1,
            "sigmaR_8" = 1,
            "sigmaR_9" = 1
        )
    )
})


# objFunction() ----------------------------------------------------
test_that("objFunction()", {
    test_initial_parameters <- c(
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        sj = 1,
        sj = 1,
        sigmaR = 1
    )

    test_data_fit <- data.frame(
        name = c(
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT"
        ),
        time = c(
            0, 5, 10, 20, 30, 60, 240, 0
        ),
        value = c(
            1.0210929,
            1.2130989,
            0.8659901,
            0.9323670,
            1.0052727,
            0.9780511,
            1.1585450,
            0.8255822
        ),
        sigma = c(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN
        ),
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo",
            "pAKT_0_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_2"
        ),
        error = c(
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1"
        )
    )
    row.names(test_data_fit) <- c(
        "pAKT1",
        "pAKT2",
        "pAKT3",
        "pAKT4",
        "pAKT5",
        "pAKT6",
        "pAKT7",
        "pAKT8"
    )




    test_parameters <- c(
        "biological" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    test_levels_list <- list(
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test", "pAKT_2"
        ),
        error = c(
            "pAKT_1"
        )
    )

    test_effectsPars <- list(
        biologicalPars = "yi",
        scalingPars = "sj",
        errorPars = "sigmaR"
    )

    test_c_strenth <- 1000

    test_all_levels <- c(
        "pAKT_0_0Uml Epo",
        "pAKT_5_0Uml Epo",
        "pAKT_10_0Uml Epo",
        "pAKT_20_0Uml Epo",
        "pAKT_30_0Uml Epo",
        "pAKT_60_0Uml Epo",
        "pAKT_240_0Uml Epo",
        "pAKT_test",
        "pAKT_2",
        "pAKT_1"
    )

    test_mask <- genMask(
        test_initial_parameters,
        test_parameters,
        test_all_levels,
        test_data_fit
    )

    test_modelExpr <- parse(text = "yi[biological]/sj[scaling]")
    test_errorModel <- parse(
        text = "(yi[biological]/sj[scaling])* sigmaR[error]"
    )

    test_modelJacobianExpr <- list(
        biological = parse(text = "1/sj[scaling]"),
        scaling = parse(text = "-(yi[biological]/sj[scaling]^2)"),
        error = parse(text = "0")
    )

    test_errorModelJacobianExpr <- list(
        biological = parse(text = "1/sj[scaling]*sigmaR[error]"),
        scaling = parse(
            text = "-(yi[biological]/sj[scaling]^2*sigmaR[error])"
        ),
        error = parse(text = "(yi[biological]/sj[scaling])")
    )

    test_constrExpr <- parse(text = "1e3 * (mean( yi ) - 1)")





    test_passParList <- list(
        parameters = test_parameters,
        fitLogscale = TRUE,
        constStrength = test_c_strenth,
        modelExpr = test_modelExpr,
        errorModelExpr = test_errorModel,
        modelJacobianExpr = test_modelJacobianExpr,
        errorModelJacobianExpr = test_errorModelJacobianExpr,
        constrExpr = test_constrExpr
    )

    test_passParList2 <- list(
        dataFit = test_data_fit,
        levelsList = test_levels_list,
        mask = test_mask,
        effectsPars = test_effectsPars
    )


    # * actual tests ----------------------------------------------------------


    expect_equal(
        objFunction(
            currentPars = test_initial_parameters,
            passParList = test_passParList,
            passParList2 = test_passParList2,
            calcDeriv = TRUE
        )$value,
        0.12445657
    )


    expect_equal(
        round(
            objFunction(
                currentPars = test_initial_parameters,
                passParList = test_passParList,
                passParList2 = test_passParList2,
                calcDeriv = TRUE
            )$gradient,
            digits = 8
        ),
        c(
            4.244916840,
            1.482979920,
            2.232102490,
            2.126117550,
            1.989399000,
            2.042934290,
            1.632636970,
            -13.463094600,
            -2.287992460,
            15.751086860
        )
    )
})


# splitForScaling -------------------------------------------------------

test_that("splitForScaling()", {
    simDataWideFile <- system.file(
        "extdata", "simDataWide.csv",
        package = "blotIt"
    )
    sim_data_long <- readWide(simDataWideFile, description = seq_len(3))

    test_effect_values <- list(
        biologicaValues = c("name", "time", "condition"),
        scalingValues = c("name", "ID"),
        errorValues = c("name")
    )

    expected_values <- c(
      116.83827140,
      138.80850037,
      99.09067869,
      106.68583790,
      115.02805006,
      111.91323007,
      132.56618485,
      94.46701749,
      122.44084967,
      51.58390606
    )

    expect_equal(
      splitForScaling(
        data = sim_data_long,
        effectsValues = test_effect_values,
        normalizeInput = FALSE
      )[[1]]$value,
      expected = expected_values
    )


    expect_equal(
      splitForScaling(
        data = sim_data_long,
        effectsValues = test_effect_values,
        normalizeInput = TRUE
      )[[1]]$value,
      expected_values / mean(expected_values)
    )

    expect_equal(
        splitForScaling(
            data = sim_data_long,
            effectsValues = test_effect_values,
            normalizeInput = FALSE
        )[[1]]$value,
        expected_values
    )
})


# resolveFunction() -----------------------------------------------------
test_that("resolveFunction()", {
    test_initial_parameters <- c(
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        sj = 1,
        sj = 1,
        sigmaR = 1
    )

    test_data_fit <- data.frame(
        name = c(
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT"
        ),
        time = c(
            0, 5, 10, 20, 30, 60, 240, 0
        ),
        value = c(
            1.0210929,
            1.2130989,
            0.8659901,
            0.9323670,
            1.0052727,
            0.9780511,
            1.1585450,
            0.8255822
        ),
        sigma = c(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN
        ),
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo",
            "pAKT_0_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_2"
        ),
        error = c(
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1"
        )
    )
    row.names(test_data_fit) <- c(
        "pAKT1",
        "pAKT2",
        "pAKT3",
        "pAKT4",
        "pAKT5",
        "pAKT6",
        "pAKT7",
        "pAKT8"
    )




    test_parameters <- c(
        "biological" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    test_levels_list <- list(
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test", "pAKT_2"
        ),
        error = c(
            "pAKT_1"
        )
    )

    test_effectsPars <- list(
        biologicalPars = "yi",
        scalingPars = "sj",
        errorPars = "sigmaR"
    )

    test_c_strenth <- 1000

    test_all_levels <- c(
        "pAKT_0_0Uml Epo",
        "pAKT_5_0Uml Epo",
        "pAKT_10_0Uml Epo",
        "pAKT_20_0Uml Epo",
        "pAKT_30_0Uml Epo",
        "pAKT_60_0Uml Epo",
        "pAKT_240_0Uml Epo",
        "pAKT_test",
        "pAKT_2",
        "pAKT_1"
    )

    test_mask <- genMask(
        test_initial_parameters,
        test_parameters,
        test_all_levels,
        test_data_fit
    )

    test_modelExpr <- parse(text = "yi[biological]/sj[scaling]")
    test_errorModel <- parse(
        text = "(yi[biological]/sj[scaling])* sigmaR[error]"
    )

    test_modelJacobianExpr <- list(
        biological = parse(text = "1/sj[scaling]"),
        scaling = parse(text = "-(yi[biological]/sj[scaling]^2)"),
        error = parse(text = "0")
    )

    test_errorModelJacobianExpr <- list(
        biological = parse(text = "1/sj[scaling]*sigmaR[error]"),
        scaling = parse(
            text = "-(yi[biological]/sj[scaling]^2*sigmaR[error])"
        ),
        error = parse(text = "(yi[biological]/sj[scaling])")
    )

    test_constrExpr <- parse(text = "1e3 * (mean( yi ) - 1)")





    test_passParList <- list(
        parameters = test_parameters,
        fitLogscale = TRUE,
        constStrength = test_c_strenth,
        modelExpr = test_modelExpr,
        errorModelExpr = test_errorModel,
        modelJacobianExpr = test_modelJacobianExpr,
        errorModelJacobianExpr = test_errorModelJacobianExpr,
        constrExpr = test_constrExpr
    )

    test_passParList2 <- list(
        dataFit = test_data_fit,
        levelsList = test_levels_list,
        mask = test_mask,
        effectsPars = test_effectsPars
    )



    # * actual tests ----------------------------------------------------------


    expect_equal(
        resolveFunction(
            currentPars = test_initial_parameters,
            passParList = test_passParList,
            passParList2 = test_passParList2,
            calcDeriv = FALSE
        )$residuals,
        c(
            -0.0210929,
            -0.2130989,
            0.1340099,
            0.0676330,
            -0.0052727,
            0.0219489,
            -0.1585450,
            0.1744178,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            0.0000000
        )
    )
})
