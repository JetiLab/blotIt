## Get example dataset
exampleWB <- system.file(
  "extdata", "INFdata_WB_wide.csv",
  package = "blotIt"
)

## read in example data by use of 'readWide()'
exampleData <- readWide(exampleWB, description = seq_len(3))

## execute alignReplicates
out <- alignReplicates(
  data = exampleData,
  model = "yi / sj",
  errorModel = "value * sigmaR",
  biological = yi ~ name + time + condition,
  scaling = sj ~ name + experiment,
  error = sigmaR ~ name,
  normalize = TRUE,
  averageTechRep = FALSE,
  verbose = FALSE,
  normalizeInput = TRUE
)
