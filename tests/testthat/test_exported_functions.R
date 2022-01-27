
# readWide() -------------------------------------------------------------


test_that("readWide works properly", {
    simDataWideFile <- system.file(
        "extdata", "simDataWide.csv",
        package = "blotIt"
    )
    simDataWide <- read.csv(simDataWideFile)

    expect_equal(
        length(
            unique(readWide(simDataWideFile, description = seq_len(3))$name)
        ),
        length(simDataWide) - 3
    )

    expect_error(
        readWide(simDataWideFile),
        "Specify columns containing descriptions."
    )

    expect_warning(
        readWide(
            simDataWideFile,
            description = c("time", "does_not_exist")
        ),
        paste0(
            "Not all columns proposed by argument 'description' are available",
            " in file.\nTaking the available ones."
        )
    )
})
