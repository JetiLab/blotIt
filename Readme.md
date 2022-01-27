# blotIt3 - a framework for alignment of biological replicate data

The present package is a rewritten version of [blotIt2](https://github.com/dkaschek/blotIt2) by Daniel Kaschek. The aim of this toolbox is to scale biological replicate data to a common scale, making the quantitative data of different gels comparable.

Please note that blotIt3 and blotIt2 *can* be used in parallel. All functions have different names, so they can not only be installed but also loaded and used simultaneously (great for double checking).

## System preperation

blotIt3 requires the `R` packages `utils, MASS, data.table, ggplot2, rootSolve` and `trust`. Additionally, the package `devtools` is needed to install blotIt3 from github. If not already done, the required packages can be installed by executing

```r
install.packages(c("utils", "MASS", "data.table", "ggplot2", "rootSolve", "trust", "devtools"))
```
blotIt3 then is installed via `devtools`:
```r
devtools::install_github("SeverinBang/blotIt3")
```
## Usage

### Data import
First, the package is imported
```r
library(blotIt3)
```
A .csv file is imported and is formatted by the function `read_wide`. An example data file is supplied. It can be accessed by 
```r
example_data_path <- system.file(
                "extdata", "sim_data_wide.csv",
                package = "blotIt3"
            )
```
This reads out the provided example file, transfers it to a temporary location and stores the path to this temporary location in `example_data_path`.
The example file is structured as follows
|time|	condition|	ID|	pAKT|	pEPOR|	pJAK2|...|
|--- | --- | --- | --- | ---|--- | ---|
|0	|0Uml Epo	|1.1	|116.838271399017|	295.836863524109| |...
|5	|0Uml Epo	|1.1	|138.808500374087|	245.229971713582| |...
|...|...|...|...|...|...|...
0	|0Uml Epo	|2	|94.4670174938645|		|293.604761934545|	...
5	|0Uml Epo	|2	|	|	|398.958892340432|	...
|...|...|...|...|...|...|...

The first three columns contain description data: time points, measurement conditions and IDs (e.g. the IDs of the different gels). All following columns contain the measurements of different targets, with the first row containing the names and the following the measurement values corresponding to the time, condition and ID stated in the first columns.

The information which columns contain descriptions has to be passed to `read_wide`:
```r
imported_data <- read_wide(
    file = example_data_path, # path to the example file
    description = seq(1,3), # Indices of columns containing the information
    sep = ",", # sign seperating the colums
    dec = "." # decimal sign
)
```
The result is then a long table of the form

|    |time| condition| ID|  name |    value|
|--- | --- | --- | --- | ---|--- |
pAKT1|       0|  0Uml Epo|  1|  pAKT| 116.83827
pAKT2|       5|  0Uml Epo|  1|  pAKT| 138.80850
pAKT3|      10|  0Uml Epo|  1|  pAKT|  99.09068
pAKT4|      20|  0Uml Epo|  1|  pAKT| 106.68584
pAKT5|      30|  0Uml Epo|  1|  pAKT| 115.02805
pAKT6|      60|  0Uml Epo|  1|  pAKT| 111.91323
pAKT7|     240|  0Uml Epo|  1|  pAKT| 132.56618
|...|...| ...| ...|...|...|

While the first (nameless) columns just contains (unique) row names. New are the columns `name` and `value`. While the column names of the original file are pasted in the former, the latter contains the respective values.
The data.frame `imported_data` can now be passed to the main function.
### Scale data
The full function call is
```r
scaled_data <- align_me(
  data = imported_data,
  model = "yi / sj",
  error_model = "value * sigmaR",
  biological = yi ~ name + time + condition,
  scaling = sj ~ name + ID,
  error = sigmaR ~ name + 1,
  parameter_fit_scale = "log",
  normalize = TRUE,
  average_techn_rep = FALSE,
  verbose = FALSE,
  normalize_input = TRUE
)
```
We will go now through the parameters individually:
- `data` A long table, usually the output of `read_wide`
- `model` A formula like describing the model used for aligning. The present one `yi / sj` means that the measured values `Y_i` are the real values `yi` scaled by scaling factors `sj`. The model therefore is the real value divided by the corresponding scaling factor.
- `error_model` A description of which errors affect the data. Here, only a relative error is present, where the parameter `sigmaR` is scaled by the respective `value`
- `biological` Description of which parameter (left hand side of the tilde) represented by which columns (right hand side of the tilde) contain the "biological effects". In the present example, the model states that the real value is represented by `yi` -- which is the left hand side of the present `biological` entry. The present right hand side is "name", "time" and "condition".
In short: we state that the entries "name", "time" and "condition" contain _real_, biological differences.
- `scaling` Same as above, but here is defined which columns contain identificators of different scaling. Here it is "name" and "ID", meaning that measurements with differ in this effects, (but have the same `biological` effects) are scaled upon another.
- `error` Describes how the error affects the values individually. The present formulation means, that the error parameter is _not_ individually adjusted.
- `parameter_fit_scale` Describes the scale on which the parameter are fitted. `align_me()` accepts "linear", "log", "log2" and "log10". The default is "Linear".
- `average_techn_rep` A logical parameter that indicates, if technical replicates should be averaged before the scaling.
- `verbose` If set to `TRUE` additional information will be printed in the console.
- `normalize_input` If set to `TRUE`, the data will be scaled before the actual scaling. This means that the raw input will be scaled to a common order of magnitude before the scaling parameters will be calculated. This is only a computational aid, to eliminate a rare fail of convergence when the different values differ by many orders of magnitude. Setting this to `TRUE` makes only sense (and is only supported) for `parameter_fit_scale = "linear"`.

The result of `align_me()` is a list with the entries
- `aligned` A `data.frame` with the columns containing the biological effects as well as the columns `value` containing the "estimated true values" and `sigma` containing the uncertainty of the fits. Both are on common 
- `scaled` The original data but with the values scaled to common scale and errors from the evaluation of the error model, also scaled to common scale (obeying Gaussian error propagation).
- `prediction` The scales and sigma are from the evaluation of the respective models (on original scale).
- `original` Just the original parameters
- `original_with_parameters` As above but with additional columns for the estimated parameters. 
- `biological` Names of the columns defined to contain the `biological` effects.
- `scaling` Names of the columns defined to contain the `scaling` effects.

### Plot Data
`blotIt3` provides _one_ plotting function `plot_align_me()` which data set will be plotted can be specified per parameter
```r
plot_align_me(
    out_list = scaled_data,
    plot_points = "aligned",
    plot_line = "aligned",
    spline = FALSE,
    scales = "free",
    align_zeros = TRUE,
    plot_caption = TRUE,
    ncol = NULL,
    my_colors = NULL,
    duplicate_zero_points = FALSE,
    my_order = NULL
)
```
The parameters again are:
- `out_list` the result of `align_me()`
- `plot_points` It can separately specified which data sets should be plotted as dots and as line. Here the data set for the dots is defined. It can be either of `original`, `scaled`, `prediction` or `aligned`.
- `plot_line` Same above but for the line.
- `spline` Logical parameter, if set to `TRUE`, the line plotted will be not straight lines connecting points but a smooth spline.
- `scales` String passed as `scales` argument to `facet_wrap`.
- `align_zeros` Logical parameter, if set to `TRUE` the zero ticks will be aligned throughout all the sub plots, although the axis can have different scales.
- `plot_caption` Logical parameter, indicating if a caption describing which data is plotted should be added to the plot.
- `ncol` Numerical passed as `ncol` argument to `facet_wrap`.
- `my_colors` list of custom color values as taken by the `values` argument in the `scale_color_manual` method for `ggplot` objects, if not set the default `ggplot` color scheme is used.
- `duplicate_zero_points` Logical, if set `TRUE` all zero time points are assumed to belong to the first condition. E.g. when the different conditions consist of treatments added at time zero. Default is `FALSE`.
- `my_order` Optional list of target names in the custom order that will be used for faceting
- `...` Logical expression used for subsetting the data frames, e.g. `name == "pAKT" & time < 60`


### Licence:
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)
