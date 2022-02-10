
# BlotIt - Optimal alignment of western blot and qPCR experiments

The present package bases on [blotIt2](https://github.com/dkaschek/blotIt2) by Daniel Kaschek. The aim of this toolbox is to align biological replicate data to a common scale, making the quantitative data of different experiments comparable.

## System preperation

blotIt requires the `R` packages `utils, MASS, data.table, ggplot2, rootSolve, data.table` and `trust`. Additionally, the package `devtools` is needed to install blotIt from github. If not already done, the required packages can be installed by executing

```r
install.packages(c("utils", "MASS", "data.table", "ggplot2", "rootSolve", "trust", "devtools", "data.table"))
```
blotIt is then installed via `devtools`:
```r
devtools::install_github("JetiLab/blotIt")
```
## Usage

### Data import
First, the package is imported
```r
library(blotIt)
```
It is assumed that measurements are present in human-readable wide formatted .csv files. An example data set can be found at
```r
exampleDataPath <- system.file(
                "extdata", "simDataWide.csv",
                package = "blotIt"
            )
```
The file has the structure

|time|	condition|	experiment|	pAKT|	pERK|...|
|--- | --- | --- | --- | ---|--- | 
|0	|0Uml Epo	| 1 |80.36|	138.14| 
|5	|0Uml Epo	| 1	|92.40|	130.79| 
|...|...|...|...|...|...
10	|0Uml Epo	| 2	|22.32|		|44.22
30	|0Uml Epo	| 2	|	|	|33.13
|...|...|...|...|...|...


The first three columns contain description data: time points, measurement conditions and experiments. All following columns contain the measurements of different targets, with the first row containing the names and the following rows comprising the measurement values corresponding to the time, condition and experiment stated in the first columns.

The information, which columns contain descriptions, has to be passed to `readWide()`:
```r
importedData <- readWide(
    file = exampleDataPath, # path to the example file
    description = seq(1, 3), # Indices of columns containing the information
    sep = ",",
    dec = "."
)
```
The result is then a long table of the form

|    |time| condition| experiment|  name |    value|
|--- | --- | --- | --- | ---|--- |
pAKT1|       0|  0Uml Epo|  1|  pAKT| 80.360534
pAKT2|       5|  0Uml Epo|  1|  pAKT| 92.402752
pAKT3|      10|  0Uml Epo|  1|  pAKT| 62.321235
pAKT4|      20|  0Uml Epo|  1|  pAKT| 85.094192
pAKT5|      30|  0Uml Epo|  1|  pAKT| 65.610456
|...|...| ...| ...|...|...|

The wide data table is melted into a `long` format, where the column `name` contains the different targets (i.e. the names of the columns not defined as `description` in `readWide`) and the column `value` comprises the respective values.

Data in this long format can then be passed to `alignReplicates()` for scaling.

### Scale data
Scaling biological replicates to one common scale has the advantage that the values, although on arbitrary scale, are comparable. The scaling approach implemented in blotIt is based on optimization:

```r
out <- alignReplicates(
  data = importedData,
  model = "yi / sj",
  errorModel = "value * sigmaR",
  biological = yi ~ name + time + condition,
  scaling = sj ~ name + experiment,
  error = sigmaR ~ name + 1,
  normalize = TRUE,
  averageTechRep = FALSE,
  verbose = FALSE,
  normalizeInput = TRUE
)
```

The resulting `out` file consists of a list containing the aligned and scaled data as well as some additional datasets.

### Plotting
The data can be plotted via
```r
P <- plotIt(
  inputList = out,
  plotPoints = "aligned",
  plotLine = "aligned",
  scales = "free",
  alignZeros = TRUE,
  plotCaption = FALSE
)
```
In `plotPoints` and `plotLine` the respective dataset chosen for plotting is defined. Many more options are available, see `?plotIt` for more details.
