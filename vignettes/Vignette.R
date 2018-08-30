### R code from vignette source '/lustre/data/FF/Carmen/BitBucket/Bioc_submission/HiCLegos/HiCLegos/vignettes/Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Vignette.Rnw:125-126
###################################################
options(width=80)


###################################################
### code chunk number 2: Vignette.Rnw:157-159 (eval = FALSE)
###################################################
## 
## Lego_list_mcool_normalisations(names.only = TRUE)


###################################################
### code chunk number 3: Vignette.Rnw:166-168 (eval = FALSE)
###################################################
## 
## Lego_list_mcool_normalisations(names.only = FALSE)


###################################################
### code chunk number 4: Vignette.Rnw:174-177 (eval = FALSE)
###################################################
## 
## mcoolName="H1-hESC-HiC-4DNFI7JNCNFB.mcool"
## Lego_list_mcool_resolutions(mcool = filename)


###################################################
### code chunk number 5: Vignette.Rnw:206-225 (eval = FALSE)
###################################################
## 
## Output.lego <- "H1-hESC-HiC-4DNFI7JNCNFB-10000-ICE-normalised-chr1.lego"
## mcool <- mcoolName
## 
## CreateLego_from_mcool(Lego = Output.lego, 
##     mcool = mcool, 
##     binsize = 10000, 
##     chrs = "chr1")
## 
## Lego_load_data_from_mcool(Lego = Output.lego, 
##     mcool = mcool, 
##     chr1 = "chr1", 
##     chr2 = "chr1", 
##     binsize = 10000,
##     cooler.batch.size = 1000000, 
##     matrix.chunk = 2000, 
##     dont.look.for.chr2 = TRUE, 
##     remove.prior = TRUE,
##     norm.factor = "Iterative-Correction")


###################################################
### code chunk number 6: Vignette.Rnw:253-290
###################################################

library("HiCLegos")

Bintable.path <- system.file("extdata", 
    "Bintable_40kb.txt", 
    package = "HiCLegos")

Chromosomes <- "chr19"
CreateLego(ChromNames = Chromosomes, 
    BinTable = Bintable.path, 
    bin.delim = " ", 
    Output.Filename = "test.hdf", 
    exec = "cat", 
    remove.existing = TRUE)

Test.mat <- matrix(NA,nrow = 800, ncol = 800)
Row <- row(Test.mat)
Col <- col(Test.mat)
Dist <- Col - Row
Matrix.file <- "Test_matrix.txt"

write.table(x = Dist, 
    file = Matrix.file, 
    sep = " ", 
    quote = FALSE,
    row.names = FALSE, 
    col.names = FALSE)

Lego.file <- "test.hdf"

Lego_load_matrix(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19",
    matrix.file = Matrix.file, 
    delim = " ", 
    exec = "cat",
    remove.prior = TRUE)


###################################################
### code chunk number 7: Vignette.Rnw:296-303
###################################################

Lego_load_cis_matrix_till_distance(Lego = Lego.file, 
    chr = "chr19", 
    matrix.file = Matrix.file, 
    delim = " ", 
    distance = 100, 
    remove.prior = TRUE)


###################################################
### code chunk number 8: Vignette.Rnw:312-318
###################################################

Lego.file <- system.file("extdata", 
    "test.hdf", 
    package = "HiCLegos")

Lego_list_rangekeys(Lego.file)


###################################################
### code chunk number 9: Vignette.Rnw:325-327
###################################################

Lego_get_bintable(Lego.file)


###################################################
### code chunk number 10: Vignette.Rnw:333-336
###################################################

Lego_get_ranges(Lego = Lego.file, 
    rangekey = "Bintable")


###################################################
### code chunk number 11: Vignette.Rnw:342-346
###################################################

Lego_get_ranges(Lego = Lego.file, 
    rangekey = "Bintable", 
    chr = "chr19")


###################################################
### code chunk number 12: Vignette.Rnw:354-358
###################################################

testRun=Lego_return_region_position(Lego = Lego.file, 
    region = "chr19:5000000:10000000")
head(testRun)


###################################################
### code chunk number 13: Vignette.Rnw:368-373
###################################################

testRun=Lego_fetch_range_index(Lego = Lego.file, 
    chr = "chr19", 
    start = 5000000, 
    end = 10000000)


###################################################
### code chunk number 14: Vignette.Rnw:395-399
###################################################

Values <- Lego_get_values_by_distance(Lego = Lego.file, 
    chr = "chr19", 
    distance = 4)


###################################################
### code chunk number 15: Vignette.Rnw:404-414
###################################################

Failsafe_median_log10 <- function(x){
    x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
    return(median(log10(x+1)))
}

Lego_get_values_by_distance(Lego = Lego.file, 
    chr = "chr19", 
    distance = 4, 
    FUN = Failsafe_median_log10)


###################################################
### code chunk number 16: Vignette.Rnw:421-432
###################################################

Failsafe_median_log10 <- function(x){
    x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
    return(median(log10(x+1)))
}

Lego_get_values_by_distance(Lego = Lego.file, 
    chr = "chr19", 
    distance = 4, 
    constrain.region = "chr19:1:5000000", 
    FUN = Failsafe_median_log10)


###################################################
### code chunk number 17: Vignette.Rnw:440-445
###################################################

Sub.matrix <- Lego_get_matrix_within_coords(Lego = Lego.file, 
    x.coords="chr19:5000001:10000000", 
    force = TRUE,
    y.coords = "chr19:5000001:10000000")


###################################################
### code chunk number 18: Vignette.Rnw:450-459
###################################################

x.axis <- 5000000/40000
y.axis <- 10000000/40000

Sub.matrix <- Lego_get_matrix(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19", 
    x.vector = c(x.axis:y.axis), 
    y.vector = c(x.axis:y.axis))


###################################################
### code chunk number 19: Vignette.Rnw:472-483
###################################################

Coordinate <- c("chr19:1:40000","chr19:40001:80000")
Lego.file <- system.file("extdata", 
    "test.hdf", 
    package = "HiCLegos")
Test_Run <- Lego_fetch_row_vector(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19", 
    by = "ranges", 
    vector = Coordinate,
    regions = c("chr19:1:1000000", "chr19:40001:2000000"))


###################################################
### code chunk number 20: Vignette.Rnw:488-499
###################################################

Coordinate <- c(1,2)
Lego.file <- system.file("extdata", 
    "test.hdf", 
    package = "HiCLegos")
Test_Run <- Lego_fetch_row_vector(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19", 
    by = "position", 
    vector = Coordinate,
    regions = c("chr19:1:1000000", "chr19:40001:2000000"))


###################################################
### code chunk number 21: Vignette.Rnw:522-524
###################################################

Lego_list_matrix_mcols()


###################################################
### code chunk number 22: Vignette.Rnw:529-538
###################################################

Lego.file <- system.file("extdata", 
    "test.hdf", 
    package = "HiCLegos")
testRun <- Lego_get_matrix_mcols(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19", 
    what = "row.sums")



###################################################
### code chunk number 23: Vignette.Rnw:548-552
###################################################

Lego_matrix_isdone(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 24: Vignette.Rnw:557-561
###################################################

Lego_matrix_issparse(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 25: Vignette.Rnw:566-570
###################################################

Lego_matrix_maxdist(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 26: Vignette.Rnw:575-579
###################################################

Lego_matrix_minmax(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 27: Vignette.Rnw:584-588
###################################################

Lego_matrix_minmax(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 28: Vignette.Rnw:593-597
###################################################

Lego_matrix_dimensions(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 29: Vignette.Rnw:602-606
###################################################

Lego_matrix_filename(Lego = Lego.file, 
    chr1 = "chr19", 
    chr2 = "chr19")


###################################################
### code chunk number 30: Vignette.Rnw:620-636
###################################################

Lego.file <- system.file("extdata", 
    "test.hdf", 
    package = "HiCLegos")


Chromosome <- "chr19"
di_window <- 10
lookup_window <- 30
TAD_ranges <- Lego_local_score_differentiator(Lego = Lego.file, 
    chrs = Chromosome, 
    di.window = di_window, 
    lookup.window = lookup_window, 
    strict = TRUE, 
    fill.gaps=TRUE, 
    chunk.size = 500)


###################################################
### code chunk number 31: Vignette.Rnw:654-662
###################################################

Name <- paste("LSD",
    di_window,
    lookup_window,
    Chromosome,sep = "_")
Lego_add_ranges(Lego = "test.hdf", 
    ranges = TAD_ranges, 
    rangekey = Name)


###################################################
### code chunk number 32: Vignette.Rnw:672-676
###################################################

Lego_list_rangekeys(Lego = "test.hdf")
TAD_ranges <- Lego_get_ranges(Lego = "test.hdf", 
    rangekey = Name)


###################################################
### code chunk number 33: chr19-5MB-10MB-normal (eval = FALSE)
###################################################
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal.pdf", 
##     Legos = Lego.file, 
##     x.coords = "chr19:5000000:10000000", 
##     y.coords = "chr19:5000000:10000000", 
##     palette = "Reds", 
##     width = 10, 
##     height = 11, 
##     return.object=TRUE)


###################################################
### code chunk number 34: chr19-5MB-10MB-normal
###################################################
Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal.pdf", 
    Legos = Lego.file, 
    x.coords = "chr19:5000000:10000000", 
    y.coords = "chr19:5000000:10000000", 
    palette = "Reds", 
    width = 10, 
    height = 11, 
    return.object=TRUE)


###################################################
### code chunk number 35: chr19-5MB-10MB-normal2 (eval = FALSE)
###################################################
## 
## Failsafe_log10 <- function(x){
##     x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
##     return(log10(x+1))
## }
## 
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal2.pdf", 
##     Legos = Lego.file, 
##     x.coords = "chr19:5000000:10000000", 
##     y.coords = "chr19:5000000:10000000", 
##     FUN = Failsafe_log10, 
##     legend.title = "Log10 Hi-C signal", 
##     palette = "Reds", 
##     width = 10, 
##     height = 11,
##     return.object=TRUE)


###################################################
### code chunk number 36: chr19-5MB-10MB-normal2
###################################################

Failsafe_log10 <- function(x){
    x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
    return(log10(x+1))
}

Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal2.pdf", 
    Legos = Lego.file, 
    x.coords = "chr19:5000000:10000000", 
    y.coords = "chr19:5000000:10000000", 
    FUN = Failsafe_log10, 
    legend.title = "Log10 Hi-C signal", 
    palette = "Reds", 
    width = 10, 
    height = 11,
    return.object=TRUE)


###################################################
### code chunk number 37: chr19-5MB-10MB-normal3 (eval = FALSE)
###################################################
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal3.pdf", 
##     Legos = Lego.file, 
##     x.coords = "chr19:5000000:10000000", 
##     y.coords = "chr19:5000000:10000000", 
##     FUN = Failsafe_log10, 
##     value.cap = 0.99, 
##     legend.title = "Log10 Hi-C signal", 
##     palette = "Reds", 
##     width = 10, 
##     height = 11, 
##     return.object=TRUE)


###################################################
### code chunk number 38: chr19-5MB-10MB-normal3
###################################################
Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal3.pdf", 
    Legos = Lego.file, 
    x.coords = "chr19:5000000:10000000", 
    y.coords = "chr19:5000000:10000000", 
    FUN = Failsafe_log10, 
    value.cap = 0.99, 
    legend.title = "Log10 Hi-C signal", 
    palette = "Reds", 
    width = 10, 
    height = 11, 
    return.object=TRUE)


###################################################
### code chunk number 39: chr19-5MB-10MB-normal4 (eval = FALSE)
###################################################
## 
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal4.pdf", 
##     Legos = Lego.file, 
##     x.coords = "chr19:5000000:10000000", 
##     y.coords = "chr19:5000000:10000000", 
##     FUN = Failsafe_log10, 
##     value.cap = 0.99, 
##     legend.title = "Log10 Hi-C signal", 
##     palette = "Reds", 
##     width = 10, 
##     height = 11, 
##     rotate = TRUE,
##     return.object=TRUE)
## 


###################################################
### code chunk number 40: chr19-5MB-10MB-normal4
###################################################

Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal4.pdf", 
    Legos = Lego.file, 
    x.coords = "chr19:5000000:10000000", 
    y.coords = "chr19:5000000:10000000", 
    FUN = Failsafe_log10, 
    value.cap = 0.99, 
    legend.title = "Log10 Hi-C signal", 
    palette = "Reds", 
    width = 10, 
    height = 11, 
    rotate = TRUE,
    return.object=TRUE)



###################################################
### code chunk number 41: chr19-5MB-10MB-normal5 (eval = FALSE)
###################################################
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal5.pdf", 
##     Legos = Lego.file, 
##     x.coords = "chr19:5000000:10000000", 
##     y.coords = "chr19:5000000:10000000", 
##     FUN = Failsafe_log10, 
##     value.cap = 0.99, 
##     legend.title = "Log10 Hi-C signal", 
##     palette = "Reds", 
##     width = 15, 
##     height = 5, 
##     rotate = TRUE,
##     return.object=TRUE)
## 


###################################################
### code chunk number 42: chr19-5MB-10MB-normal5
###################################################
Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal5.pdf", 
    Legos = Lego.file, 
    x.coords = "chr19:5000000:10000000", 
    y.coords = "chr19:5000000:10000000", 
    FUN = Failsafe_log10, 
    value.cap = 0.99, 
    legend.title = "Log10 Hi-C signal", 
    palette = "Reds", 
    width = 15, 
    height = 5, 
    rotate = TRUE,
    return.object=TRUE)



###################################################
### code chunk number 43: chr19-5MB-10MB-normal6 (eval = FALSE)
###################################################
## 
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal6.pdf", 
##    Legos = Lego.file, 
##    tad.ranges = TAD_ranges, 
##    x.coords = "chr19:5000000:10000000", 
##    y.coords = "chr19:5000000:10000000", 
##    colours = "#E0CA3C", 
##    FUN = Failsafe_log10, 
##    value.cap = 0.99, 
##    legend.title = "Log10 Hi-C signal", 
##    palette = "Reds", 
##    width = 10, 
##    height = 11,
##    return.object=TRUE)


###################################################
### code chunk number 44: chr19-5MB-10MB-normal6
###################################################

Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal6.pdf", 
   Legos = Lego.file, 
   tad.ranges = TAD_ranges, 
   x.coords = "chr19:5000000:10000000", 
   y.coords = "chr19:5000000:10000000", 
   colours = "#E0CA3C", 
   FUN = Failsafe_log10, 
   value.cap = 0.99, 
   legend.title = "Log10 Hi-C signal", 
   palette = "Reds", 
   width = 10, 
   height = 11,
   return.object=TRUE)


###################################################
### code chunk number 45: chr19-5MB-10MB-normal7 (eval = FALSE)
###################################################
## Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal.pdf", 
##     Legos = Lego.file, 
##     tad.ranges = TAD_ranges,
##     x.coords = "chr19:5000000:10000000", 
##     y.coords = "chr19:5000000:10000000", 
##     colours = "#E0CA3C",
##     FUN = Failsafe_log10, 
##     value.cap = 0.99, 
##     legend.title = "Log10 Hi-C signal", 
##     palette = "Reds", 
##     width = 15, 
##     height = 9,
##     cut.corners = TRUE,
##     rotate = TRUE,
##     return.object=TRUE) 


###################################################
### code chunk number 46: chr19-5MB-10MB-normal7
###################################################
Lego_vizart_plot_heatmap(File = "chr19-5MB-10MB-normal5.pdf", 
    Legos = Lego.file, 
    x.coords = "chr19:5000000:10000000", 
    y.coords = "chr19:5000000:10000000", 
    FUN = Failsafe_log10, 
    value.cap = 0.99, 
    legend.title = "Log10 Hi-C signal", 
    palette = "Reds", 
    width = 15, 
    height = 5, 
    rotate = TRUE,
    return.object=TRUE)



