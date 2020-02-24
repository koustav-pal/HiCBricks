# HDF structure

#' Create the entire HDF5 structure and load the bintable
#'
#' `Create_many_Bricks` creates the HDF file and returns a BrickContainer
#'
#' This function creates the complete HDF data structure, loads the binning
#' table associated to the Hi-C experiment, creates a 2D matrix
#' layout for all specified chromosome pairs and creates a json file for the
#' project. At the end, this function will return a S4 object of class 
#' BrickContainer.  **Please note**, the binning table must be a 
#' discontinuous one (first range end != secode range start), 
#' as ranges overlaps using the "any" form will routinely identify adjacent 
#' ranges with the same end and start to be in the overlap. Therefore, this 
#' criteria is enforced as default behaviour.
#'
#' @param BinTable \strong{Required}
#' A string containing the path to the file to load as the binning table for
#' the Hi-C experiment. The number of entries per chromosome defines the
#' dimension of the associated Hi-C data matrices. For example, if chr1
#' contains 250 entries in the binning table, the _cis_ Hi-C data matrix for
#' chr1 will be expected to contain 250 rows and 250 cols. Similary, if the
#' same binning table contained 150 entries for chr2, the _trans_ Hi-C
#' matrices for chr1,chr2 will be a matrix with dimension 250 rows and
#' 150 cols.
#'
#' There are no constraints on the bintable format. As long as the table is
#' in a delimited format, the corresponding table columns can be outlined with
#' the associated parameters. The columns of importance are chr, start and end.
#'
#' It is recommended to always use binning tables where the end and start of
#' consecutive ranges are not the same. If they are the same, this may lead to
#' \strong{unexpected behaviour} when using the GenomicRanges "any" overlap
#' function.
#'
#' @param output_directory \strong{Required}
#' A string specifying the location where the HDF files will be created.
#'
#' @param file_prefix \strong{Required}
#' A string specifying the prefix that is concatenated to the hdf files stored
#' in the output_directory.
#' 
#' @param bin_delim \strong{Optional}. Defaults to tabs.
#' A character vector of length 1 specifying the delimiter used in the file
#' containing the binning table.
#'
#' @param col_index \strong{Optional}. Default "c(1,2,3)".
#' A character vector of length 3 containing the indexes of the required
#' columns in the binning table. the first index, corresponds to the chr
#' column, the second to the start column and the third to the end column.
#'
#' @param impose_discontinuity \strong{Optional}. Default TRUE.
#' If TRUE, this parameter ensures a check to make sure that required the end
#' and start coordinates of consecutive entries are not the same per
#' chromosome.
#'
#' @param hdf_chunksize \strong{Optional}.
#' A numeric vector of length 1. If provided, the HDF dataset will use this
#' value as the chunk size, for all matrices. By default, the ChunkSize is
#' set to matrix dimensions/100.
#'
#' @param remove_existing \strong{Optional}. Default FALSE.
#' If TRUE, will remove the HDF file with the same name and create a new one.
#' By default, it will not replace existing files.
#' 
#' @param link_existing \strong{Optional}. Default FALSE.
#' If TRUE, will re-add the HDF file with the same name.
#' By default, this parameter is set to FALSE.
#'
#' @param experiment_name \strong{Optional}.
#' If provided, this will be the experiment name for the BrickContainer.
#' 
#' @param resolution \strong{required}.
#' A value of length 1 of class character or numeric specifying the resolution 
#' of the Hi-C data loaded.
#' 
#' @param type \strong{optional}. Default any
#' A value from one of any, cis, trans specifying the type of matrices to load.
#' Any will load both cis (intra-choromosomal, e.g. chr1 vs chr1) and trans (
#' inter-chromosomal, e.g. chr1 vs chr2) Hi-C matrices. Whereas cis and trans 
#' will load either cis or trans Hi-C matrices.
#' 
#' @details The structure of the HDF file is as follows:
#' The structure contains three major groups which are then hierarchically
#' nested with other groups to finally lead to the corresponding datasets. 
#' \itemize{
#'    \item Base.matrices - \strong{group} For storing Hi-C matrices
#'    \itemize{
#'        \item chromosome - \strong{group}
#'        \item chromosome - \strong{group}
#'        \itemize{
#'            \item attributes - \strong{attribute}
#'            \itemize{
#'                \item Filename - Name of the file
#'                \item Min - min value of Hi-C matrix
#'                \item Max - max value of Hi-C matrix
#'                \item sparsity - specifies if this is a sparse matrix
#'                \item distance - max distance of data from main diagonal
#'                \item Done - specifies if a matrix has been loaded
#'            }
#'            \item matrix - \strong{dataset} - contains the matrix
#'            \item chr1_bin_coverage - \strong{dataset} - proportion of row 
#' cells with values greater than 0
#'            \item chr1_row_sums - \strong{dataset} - total sum of all values
#' in a row
#'            \item chr2_col_sums - \strong{dataset} - total sum of all values
#' in a col
#'            \item chr2_bin_coverage - \strong{dataset} - proportion of col 
#' cells with values greater than 0
#'            \item sparsity - \strong{dataset} - proportion of non-zero cells
#' near the diagonal
#'        }
#'    }
#'    \item Base.ranges - \strong{group}, Ranges tables for quick and easy
#' access. Additional ranges tables are added here under separate group names.
#'    \itemize{
#'        \item Bintable - \strong{group} - The main binning table associated
#' to a Brick.
#'        \itemize{
#'            \item ranges - \strong{dataset} - Contains the three main columns
#' chr, start and end.
#'            \item offsets - \strong{dataset} - first occurence of any given
#' chromosome in the ranges dataset.
#'            \item lengths - \strong{dataset} - Number of occurences of that
#' chromosome
#'            \item chr.names - \strong{dataset} - What chromosomes are present
#' in the given ranges table.
#'        }
#'    }
#'    \item Base.metadata - \strong{group}, A place to store metadata info
#'    \itemize{
#'        \item chromosomes - \strong{dataset} - Metadata information
#' specifying the chromosomes present in this particular Brick file.
#'        \item other metadata tables.
#'    }
#'}
#' @return This function will generate the target Brick file. Upon completion,
#' the function will return an object of class BrickContainer.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' out_dir <- file.path(tempdir(), "Creator_test")
#' dir.create(out_dir)
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' @details
#' Keep in mind that if the end coordinates and start coordinates of adjacent
#' ranges are not separated by at least a value of 1, then
#' impose.discontinuity = TRUE will likely cause an error to occur.
#' This may seem obnoxious, but GenomicRanges by default will consider an
#' overlap of 1 bp as an overlap. Therefore, to be certain that ranges which
#' should not be, are not being targeted during retrieval operations, a check
#' is initiated to make sure that adjacent ends and starts are not
#' overlapping.
#' To load continuous ranges, use impose.discontinuity = FALSE.
#'
#' Also note, that col.index determines which columns to use for chr, start
#' and end. Therefore, the original binning table may have 10 or 20 columns,
#' but it only requires the first three in order of chr, start and end.
#'
Create_many_Bricks <- function(BinTable, bin_delim="\t", col_index=c(1,2,3), 
    impose_discontinuity=TRUE, hdf_chunksize=NULL, output_directory = NA, 
    file_prefix = NA, remove_existing = FALSE, 
    link_existing = FALSE, experiment_name = NA, resolution = NA, 
    type = c("both", "cis", "trans")){

    H5close()
    type = match.arg(type)
    Reference.object <- GenomicMatrix$new()
    if(is.na(output_directory)){
        stop("Parameter output_directory must be provided")
    } else if(output_directory == "."){
        Err.msg <- paste("Unable to create in directory", output_directory, 
            "without an explicit path", "definition!", "If you want to create",
            "files inside the", "current working directory,",
            "please provide", "the complete path. Or, you can",
            "provide as an", "argument,", 
            paste("file.path(getwd())", sep = ""), 
        sep = " ")
        stop(._format_by_charlen(string = Err.msg))
    } else if(!dir.exists(output_directory)){
        stop("output_directory does not exist")
    }
    output_directory <- normalizePath(output_directory)
    if(is.na(resolution)){
        stop("Parameter resolution must be provided")
    }
    if(is.na(file_prefix)){
        stop("Parameter file_prefix must be provided")
    }
    if(is.null(BinTable)){
        stop("Variable Bintable cannot be empty. Binning information must be ",
            "provided at startup\n")
    }
    if(is.na(experiment_name)){
        Dated <- format(Sys.Date(),"%Y-%m-%d")
        experiment_name = paste("A HiCBricks object created on",Dated)
    }

    resolution <- .format_resolution(resolution)
    Config_filepath <- .make_configuration_path(output_directory)
    Configuration_matrix_list <- list()

    Bintable.list <- Read_bintable(Filename = BinTable, 
        read.delim = bin_delim, col.index = col_index, 
        impose.discontinuity = impose_discontinuity)
    bintable_df <- Bintable.list[['main.tab']]
    # Create the 0 level directories in the HDF file
    ChromosomeList <- unique(bintable_df[,'chr'])
    Chrom_info_df <- return_chrominfo_df(bintable_df = bintable_df, 
        chromosomes = ChromosomeList)
    chromosome_lengths <- Chrom_info_df$size
    Resolutions <- resolution

    if(file.exists(Config_filepath)){
        message("Sourcing Container parameters from existing json ",
            "at location ", output_directory)
        message("Ignoring new parameters")
        Container <- load_BrickContainer(Config_filepath)
        Configuration_header <- return_configuration_header(Container)

        file_prefix <- Configuration_header$file_prefix
        output_directory <- Configuration_header$project_directory
        experiment_name <- BrickContainer_list_experiment_name(Container)
        Resolutions <- BrickContainer_list_resolutions(Container)
        chromosome_df <- BrickContainer_list_chromosomes(Container, 
            lengths = TRUE)
        ChromosomeList <- chromosome_df$chrom
        chromosome_lengths <- chromosome_df$lengths
        # Files_list <- BrickContainer_list_files(Container)
        Configuration_matrix_list <- return_configuration_matrix_info(
            Container)
        if(!remove_existing){
            if(any(grepl(pattern = resolution, 
                    x = Resolutions, 
                    ignore.case = TRUE))){
                stop("Resolution ",resolution,
                    " is already present.")
            }
        }
        Resolutions <- unique(c(Resolutions, resolution))
    }
    
    Configuration_header <- .create_configuration_header(
        file_prefix = file_prefix, 
        output_directory = output_directory, 
        experiment_name = experiment_name, 
        resolution = Resolutions, 
        chromosomes = ChromosomeList, 
        chromosome_lengths = chromosome_lengths)

    Chromosome.pairs.list <- return_chromosome_pairs(
        chromosomes = ChromosomeList,
        type = type)
    Root.folders <- Reference.object$GetRootFolders()

    for (chrom1 in names(Chromosome.pairs.list)) {
        for (chrom2 in Chromosome.pairs.list[[chrom1]]) {
            hdf_filename <- paste(paste(file_prefix, 
                resolution,
                chrom1, "vs", chrom2, sep = "_"), 
            Reference.object$brick.extension, sep = ".")
            Configuration_matrix_list[[paste(chrom1, chrom2, 
                resolution, sep = "_")]] <- 
                .create_brick(output_directory = output_directory, 
                    filename = hdf_filename, 
                    chrom1 = chrom1, 
                    chrom2 = chrom2, 
                    resolution = resolution, 
                    bintable_df = bintable_df, 
                    hdf_chunksize = hdf_chunksize, 
                    remove_existing = remove_existing,
                    link_existing = link_existing)
        }
    }
    Container <- .prepare_BrickContainer(Configuration_header, 
        Configuration_matrix_list, 
        Config_filepath)
    .write_configuration_file(Container, Config_filepath)
    return(Container)
}

#' Create the entire HDF5 structure and load the bintable from a mcool file
#'
#' `Create_many_Bricks_from_mcool` is a wrapper on Create_many_Bricks which 
#' creates the Brick data structure from an mcool file.
#'
#' mcool are a standard 4D nucleome data structure for Hi-C data. Read more
#' about the 4D nucleome project \href{https://data.4dnucleome.org/}{here}.
#'
#' @inheritParams Brick_load_data_from_mcool
#'
#' @inheritParams Create_many_Bricks
#'
#' @return This function will generate the target Brick file. Upon completion,
#' the function will provide the path to the created/tracked HDF file.
#'
#' @examples
#'
#' \dontrun{
#' require(curl)
#' out_dir <- file.path(tempdir(),"mcool_test_dir")
#' dir.create(path = out_dir)
#' curl_download(url = paste("https://data.4dnucleome.org/",
#' "files-processed/4DNFI7JNCNFB/",
#' "@@download/4DNFI7JNCNFB.mcool", sep = ""),
#' destfile = file.path(out_dir,"H1-hESC-HiC-4DNFI7JNCNFB.mcool"))
#'
#' mcool <- file.path(out_dir,"H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' 
#' Create_many_Bricks_from_mcool(output_directory = out_dir,
#' file_prefix = "Test",
#' mcool = mcool,
#' resolution = 50000,
#' experiment_name = "A random 4DN dataset")
#'
#' }
#'
#' @seealso \code{\link{Brick_load_data_from_mcool}} to load data from
#' the mcool to a Brick store.
#'
Create_many_Bricks_from_mcool <- function(output_directory = NA, 
    file_prefix = NA, mcool = NULL, resolution = NULL, 
    experiment_name = NA, remove_existing = FALSE){

    Reference.object <- GenomicMatrix$new()
    if(is.null(mcool)){
        stop("mcool must be provided as mcool= /path/to/something")
    }
    if(!file.exists(mcool)){
        stop("mcool not found!")   
    }
    resolutions <- Brick_list_mcool_resolutions(mcool = mcool)
    if(!is.null(resolutions)){
        mcool.version <- GetAttributes(Path = Create_Path(
            c(Reference.object$mcool.resolutions.name, resolutions[1])), 
            File=mcool, Attributes="format-version", on = "group",
            ignore.fun.cast = TRUE)[,"format-version"]
        if(is.null(resolution)){
            stop("resolution cannot be NULL when resolutions are present..\n")
        }
        if(length(resolution) > 1){
            stop("resolution cannot have more than one value\n")
        }
        resolution <- .format_resolution(resolution)
        if(!(resolution %in% resolutions)){
            stop("all resolutions were not found in this mcool file. See all",
                " resolutions available with Brick_list_mcool_resolutions\n")
        }
    }else{
        mcool.version <- GetAttributes(Path = NULL, File=mcool, 
            Attributes="format-version", on = "file",
            ignore.fun.cast = TRUE)[,"format-version"]
    }
    mcool.version <- as.numeric(as.character(mcool.version))
    message("Provided mcool is a version ", mcool.version," file.")
    cooler.remap.chrom <- ._mcool_remap_chromosomes(File = mcool,
        mcool.version = mcool.version, resolution = !is.null(resolutions),
        binsize = resolution)
    ChromNames <- cooler.remap.chrom[,"chr.name"]
    mcool_bintable_ranges <- ._mcool_bintable_ranges(mcool.file = mcool,
        resolution = !is.null(resolutions),
        mcool.remap.chrom = cooler.remap.chrom, 
        binsize = resolution,
        mcool.version = mcool.version)
    RetVar <- Create_many_Bricks(BinTable = mcool_bintable_ranges, 
        output_directory = output_directory, file_prefix = file_prefix, 
        resolution = resolution, experiment_name = experiment_name, 
        remove_existing = remove_existing)
    return(RetVar)
}

#' Get all available normalisations in an mcool file.
#'
#' `Brick_list_mcool_resolutions` lists all available resolutions in the mcool
#' file.
#'
#' @param mcool \strong{Required}.
#' A parameter specifying the name of an mcool file
#'
#' @return A named vector listing all possible resolutions in the file.
#' 
#' @examples
#' 
#' \dontrun{
#' require(curl)
#' out_dir <- file.path(tempdir(),"mcool_test_dir")
#' dir.create(path = out_dir)
#' 
#' curl_download(url = paste("https://data.4dnucleome.org/",
#' "files-processed/4DNFI7JNCNFB/",
#' "@@download/4DNFI7JNCNFB.mcool", sep = ""),
#' destfile = file.path(out_dir,"H1-hESC-HiC-4DNFI7JNCNFB.mcool"))
#'
#' mcool <- file.path(out_dir,"H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' 
#' Brick_list_mcool_resolutions(mcool)
#'
#' }
Brick_list_mcool_resolutions <- function(mcool){
    return(mcool_list_resolutions(mcool = mcool))
}

#' Get all available normalisations in an mcool file.
#'
#' `Brick_list_mcool_normalisations` lists the names available for
#' accessing the various normalisation factors in an mcool file. Please note,
#' this only lists the mapping of the columns to their respective names.
#' It does not check for the availability of that particular column in
#' the mcool file
#'
#' @param names.only \strong{Optional}. Default FALSE
#' A parameter specifying whether to list only the human readable names without
#' their respective column names in the mcool file.
#'
#' @return A named vector listing all possible normalisation factors.
#'
#' @examples
#' Brick_list_mcool_normalisations()
#'
Brick_list_mcool_normalisations <- function(names.only = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(names.only){
        return(names(Reference.object$mcool.available.normalisations()))
    }
    return(Reference.object$mcool.available.normalisations())
}

#' Check if a normalisation exists in an mcool file.
#'
#' `Brick_mcool_normalisation_exists` checks if a particular normalisation
#' exists in an mcool file.
#'
#' @inheritParams Brick_load_data_from_mcool
#'
#' @param norm_factor \strong{Required}.
#' The normalization factor to use for normalization from an mcool file.
#' norm_factor currently accepts one of "Iterative-Correction", "Knight-Ruitz",
#' "Vanilla-coverage", "Vanilla-coverage-square-root".
#' 
#' @return A boolean vector of length 1
#'
#' @examples
#'
#' \dontrun{
#'
#' require(curl)
#' out_dir <- file.path(tempdir(), "mcool_test_dir")
#' dir.create(path = out_dir)
#' curl_download(url = paste("https://data.4dnucleome.org/",
#' "files-processed/4DNFI7JNCNFB/",
#' "@@download/4DNFI7JNCNFB.mcool", sep = ""),
#' destfile = file.path(out_dir, "H1-hESC-HiC-4DNFI7JNCNFB.mcool"))
#'
#' mcool <- file.path(out_dir, "H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' Brick_mcool_normalisation_exists(mcool = mcool,
#' norm_factor = "Iterative-Correction",
#' resolution = 50000)
#'
#' }
#'
Brick_mcool_normalisation_exists <- function(mcool, norm_factor = NULL,
    resolution = NULL){
    Reference.object <- GenomicMatrix$new()
    Norm.factors <- Brick_list_mcool_normalisations()
    Norm.factor <- Norm.factors[norm_factor]
    if(length(Norm.factor)!= 1){
        stop("Please check the available norm factors with ",
            "Brick_list_mcool_normalisations.\n")
    }
    names(Norm.factor) <- NULL
    mcool.version <- GetAttributes(Path = NULL, File=mcool,
        Attributes="format-version", on = "file",
        ignore.fun.cast = TRUE)[,"format-version"]
    mcool.version <- as.numeric(as.character(mcool.version))
    Bintable.keys <- Reference.object$mcool.bintable.keys(
        version = mcool.version)
    Bintable.group <- Bintable.keys[1]
    resolutions <- Brick_list_mcool_resolutions(mcool = mcool)
    if(!is.null(resolutions) & is.null(resolution)){
        stop("resolution must be provided when",
            " different resolutions are present",
            " in an mcool file.\n")
    }
    resolution <- .format_resolution(resolution)
    if(!is.null(resolutions) & !is.null(resolution)){
        if(!(resolution %in% resolutions)){
            stop("resolution not found in mcool file.",
                " Please check available resolutions",
                " with Brick_list_mcool_resolutions.\n")
        }
    }
    if(!is.null(resolutions)){
        Bintable.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,
                resolution,Bintable.group)) 
    }else{
        Bintable.group.path <- Create_Path(Bintable.group)
    }
    Handler <- ._Brick_Get_Something_(Group.path = Bintable.group.path,
        Brick = mcool, return.what = "group_handle")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    CloseH5Con(Handle = Handler, type = "group")
    return(Norm.factor %in% GroupList)
}

#' Get the chrominfo for the Hi-C experiment.
#'
#' `Brick_get_chrominfo` fetches the associated chrominfo table for the
#' Brick it is associated to.
#'
#' @param Brick \strong{Required}.
#' A string specifying the path to the Brick store created with 
#' Create_many_Brick.
#'
#' @inheritParams Brick_add_ranges
#' 
#' @return A three column data.frame containing chromosomes, nrows and length.
#'
#' chromosomes corresponds to all chromosomes in the provided bintable.
#'
#' nrows corresponds to the number of entries in the bintable or dimension
#' for that chromosome in a Hi-C matrix.
#'
#' Length is the total bp length of the same
#' chromosome (max value for that chromosome in the bintable).
#'
#' @examples
#' 
#' Bintable_path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "HiCBricks_chrominfo_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable_path, 
#' bin_delim=" ", remove_existing=TRUE, output_directory = out_dir, 
#' file_prefix = "HiCBricks_vignette_test", resolution = 100000,
#' experiment_name = "HiCBricks vignette test")
#' 
#' Brick_get_chrominfo(Brick = My_BrickContainer, resolution = 100000)
#' 
Brick_get_chrominfo <- function(Brick, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    configuration_null_check(resolution, "resolution")
    configuration_na_check(resolution, "resolution")
    configuration_length_check(resolution, "resolution", 1)
    resolution <- .format_resolution(resolution)
    Matrix_info <- return_configuration_matrix_info(Brick)
    current_resolution <- vapply(Matrix_info, function(a_list){
        a_list$resolution == resolution
    },TRUE)
    if(!any(current_resolution)){
        stop(resolution," not found in provided BrickContainer")
    }
    chrom1_binned_length <- vapply(Matrix_info[current_resolution], 
        function(a_list){
        a_list$dimensions[1]
    }, 100)
    chrom1s <- vapply(Matrix_info[current_resolution], 
        function(a_list){
        a_list$chrom1
    }, "chr1")
    chrom1_max_sizes <- vapply(Matrix_info[current_resolution], 
        function(a_list){
        a_list$lengths[1]
    }, 100) 
    chrom1_not_duplicated <- !duplicated(chrom1s)
    Chrom_info_df <- data.frame(chr = chrom1s[chrom1_not_duplicated],
    nrow = chrom1_binned_length[chrom1_not_duplicated],
    size = chrom1_max_sizes[chrom1_not_duplicated],
    stringsAsFactors = FALSE)
    rownames(Chrom_info_df) <- NULL
    return(Chrom_info_df)
}

#' Creates a ranges object from provided vectors.
#'
#' `Brick_make_ranges` creates a GRanges object from the provided arguments
#' 
#' @param chrom \strong{Required}.
#' A 1 dimensional character vector of size N specifying the chromosomes in the
#' ranges.
#'
#' @param start \strong{Required}.
#' A 1 dimensional numeric vector of size N specifying the start positions in
#' the ranges.
#'
#' @param end \strong{Required}.
#' A 1 dimensional numeric vector of size N specifying the end positions in
#' the ranges. Must be less than Start.
#'
#' @param strand \strong{Optional}.
#' A 1 dimensional character vector of size N specifying the strand of the
#' ranges. If not provided, this will be set to the default *.
#'
#' @param names \strong{Optional}.
#' A 1 dimensional character vector of size N specifying the names of the
#' ranges. If not provided, this will be set to the default chr:start:end.
#'
#' @return A GenomicRanges object with the previous sort order being preserved
#'
#' @examples
#'
#' Chrom <- c("chrS","chrS","chrS","chrS","chrS")
#' Start <- c(10000,20000,40000,50000,60000)
#' End <- c(10001,20001,40001,50001,60001)
#' Test_ranges <- Brick_make_ranges(chrom = Chrom, start = Start, end = End)
#'
Brick_make_ranges = function(chrom, start, end, strand = NA,
    names = NA){
    Reference.object <- GenomicMatrix$new()
    if(length(names) == 0 | any(is.na(names)) | any(is.null(names))){
        names <- paste(chrom, as.integer(start), as.integer(end),
            sep = Reference.object$Ranges.separator)
    }
    if(all(is.na(strand))){
        strand <- rep("*", length(chrom))
    }
    Object<-GenomicRanges::GRanges(
        seqnames=Rle(chrom),
        ranges=IRanges(start,end= end,names= names),
        strand=Rle(strand(strand)))
    return(Object)
}

#' Store a ranges object in the Brick store.
#'
#' `Brick_add_ranges` loads a GRanges object into the Brick store.
#'
#' With this function it is possible to associate other ranges objects with the
#' Brick store. If metadata columns are present, the are also loaded into the
#' Brick store. Although not explicitly asked for, the metadata columns should
#' not be of type list as this may create complications down the line. We
#' ask for ranges objects, so if the same ranges object is later retrieved
#' two additional columns will be present. These are the strand and width
#' columns that are obtained when a ranges is converted into a data.frame.
#' Users can ignore these columns.
#' 
#' @param ranges \strong{Required}.
#' An object of class ranges specifying the ranges to store in the Brick.
#'
#' @param rangekey \strong{Required}.
#' The name to use for the ranges within the Brick store.
#'
#' @param resolution \strong{Optional}. Default NA
#' When an object of class BrickContainer is provided, resolution defines the 
#' resolution on which the function is executed
#' 
#' @param all_resolutions \strong{Optional}. Default FALSE
#' If resolution is not defined and all_resolutions is TRUE, the resolution 
#' parameter will be ignored and the function is executed on all files listed 
#' in the Brick container
#' 
#' @param num_cpus \strong{Optional}. Default 1
#' When an object of class BrickContainer is provided, num_cpus defines the
#' maximum number of parallel jobs that will be run.
#' 
#' @inheritParams Brick_get_chrominfo
#'
#' @return Returns TRUE if completed successfully.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "add_ranges_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#'
#' Chrom <- c("chrS","chrS","chrS","chrS","chrS")
#' Start <- c(10000,20000,40000,50000,60000)
#' End <- c(10001,20001,40001,50001,60001)
#' Test_ranges <- Brick_make_ranges(chrom = Chrom, start = Start, end = End)
#' Brick_add_ranges(Brick = My_BrickContainer, ranges = Test_ranges,
#' rangekey = "test_ranges", all_resolutions = TRUE)
#'
Brick_add_ranges = function(Brick, ranges, rangekey, resolution = NA,
    all_resolutions = FALSE, num_cpus = 1){
    Reference.object <- GenomicMatrix$new()
    if(!(class(ranges) %in% "GRanges") | ("list" %in% class(ranges))){
        stop("Object of class Ranges expected")
    }
    BrickContainer_resolution_check(resolution, all_resolutions)
    configuration_length_check(rangekey, "rangekey", 1)
    .check_if_rangekey_exists_and_do_something(Brick, rangekey, 
        issue_stop = TRUE, 
        resolution = resolution,
        all_resolutions = all_resolutions,
        message = paste("rangekey already exists!",
            "Unable to overwrite existing rangekeys!"))
    Brick_paths <- BrickContainer_get_path_to_file(Brick = Brick,
        resolution = resolution)
    ranges_list <- .prepare_ranges_df_for_add(ranges)
    Ranges.df.coords <- ranges_list[["ranges"]]
    Metadata.list <- ranges_list[["metadata"]]
    Param <- .get_instance_biocparallel(workers = num_cpus)
    bplapply(Brick_paths, function(Brick_path){
        ._Brick_Add_Ranges_(Brick = Brick_path,
        Group.path = Create_Path(c(Reference.object$hdf.ranges.root,
            rangekey)),
        ranges.df = Ranges.df.coords,
        name = rangekey,
        mcol.list = Metadata.list)
    }, BPPARAM = Param)
    return(TRUE)
}

#' List the matrix pairs present in the Brick store.
#'
#' `Brick_list_matrices` will list all chromosomal pair matrices from the 
#' Brick store, with their associated filename, value range, done status 
#' and sparse
#'
#' @inheritParams Brick_get_chrominfo
#' @inheritParams Brick_load_matrix
#' @inheritParams Brick_add_ranges
#'
#' @return Returns a data.frame object with columns chr1, chr2 corresponding
#' to chromosome pairs, and the associated attributes. filename corresponds to
#' the name of the file that was loaded for the pair. min and max specify the
#' minimum and maximum values in the matrix, done is a logical value
#' specifying if a matrix has been loaded and sparsity specifies if a matrix
#' is defined as a sparse matrix.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_matrices_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Brick_list_matrices(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#'
Brick_list_matrices = function(Brick, chr1 = NA, chr2 = NA, resolution = NA, 
    all_resolutions = FALSE){
    Reference.object <- GenomicMatrix$new()
    configuration_na_check(chr1, "chr1")
    configuration_na_check(chr2, "chr2")
    BrickContainer_resolution_check(resolution, all_resolutions)
    Brick_filepaths <- BrickContainer_list_files(Brick = Brick, 
    chr1 = chr1, chr2 = chr2, resolution = resolution)
    Colnames <- Reference.object$matrices.chrom.attributes
    chr1.list <- lapply(seq_len(nrow(Brick_filepaths)), function(i){
        Row <- Brick_filepaths[i,]
        Values <- GetAttributes(
            Path = Create_Path(
                c(Reference.object$hdf.matrices.root,
                    Row$chrom1,
                    Row$chrom2)),
            File = Row$filepaths,
            Attributes = Colnames,
            on = "group")
        temp.df <- cbind(data.frame("chr1" = Row$chrom1, 
            "chr2" = Row$chrom2, resolution = Row$resolution), Values)
    })
    Matrix.list.df <- do.call(rbind,chr1.list)
    rownames(Matrix.list.df) <- NULL
    return(Matrix.list.df)
}


#' List the ranges tables stored within the Brick.
#'
#' `Brick_list_rangekeys` lists the names of all ranges associated to a Brick.
#' 
#' @return A one dimensional character vector of length x specifying the names
#' of all ranges currently present in the file.
#'
#' @inheritParams Brick_get_chrominfo
#' @inheritParams Brick_add_ranges
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_rangekeys_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Brick_list_rangekeys(Brick = My_BrickContainer, resolution = 100000)
#'
Brick_list_rangekeys = function(Brick, resolution = NA, 
    all_resolutions = FALSE){
    Reference.object <- GenomicMatrix$new()
    BrickContainer_resolution_check(resolution, all_resolutions)
    Brick_filepaths <- BrickContainer_get_path_to_file(Brick,
    resolution = resolution)
    Handler <- ._Brick_Get_Something_(
        Group.path = Create_Path(Reference.object$hdf.ranges.root),
        Brick = Brick_filepaths[1], return.what = "group_handle")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    return(GroupList)
}

#' Check to see if the Brick contains a ranges with a certain name.
#'
#' `Brick_rangekey_exists` checks for the presence of a particular ranges with
#' a certain name.
#'
#' @param rangekey \strong{Required}.
#' A string specifying the name of the ranges to check for.
#'
#' @inheritParams Brick_get_chrominfo
#' @inheritParams Brick_list_rangekeys
#'
#' @return A logical vector of length 1 with either TRUE or FALSE values.
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_rangekeys_exists_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' Brick_rangekey_exists(Brick = My_BrickContainer, rangekey = "Bintable",
#' resolution = 100000)
#'
Brick_rangekey_exists = function(Brick, rangekey, resolution = NA, 
    all_resolutions = FALSE){
    Keys <- Brick_list_rangekeys(Brick = Brick, 
        resolution = resolution, 
        all_resolutions = all_resolutions)
    return(rangekey %in% Keys)
}

#' Find out what metadata columns are associated to a ranges with a certain
#' name
#'
#' `Brick_list_ranges_mcols` will list the metadata columns of the specified
#' ranges if it is present in the Brick store.
#'
#' @param rangekey \strong{Optional}.
#' A string specifying the name of the ranges. If not present, the metadata
#' columns of all ranges will be listed.
#' 
#' @inheritParams Brick_get_chrominfo
#'
#' @return if no metadata columns are present, NA. If metadata columns are
#' present, a data.frame object containing the name of the ranges and the
#' associated metadata column name.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_ranges_mcols_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Brick_list_ranges_mcols(Brick = My_BrickContainer, rangekey = "Bintable",
#' resolution = 100000)
#'
Brick_list_ranges_mcols = function(Brick, rangekey = NULL, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    configuration_na_check(resolution, "resolution")
    configuration_length_check(resolution, "resolution", 1)
    Brick_filepath <- BrickContainer_get_path_to_file(Brick,
        resolution = resolution)[1]
    RangeKeys <- Brick_list_rangekeys(Brick = Brick, resolution = resolution)
    if(!is.null(rangekey)){
        .check_if_rangekey_not_exists_and_do_something(Brick = Brick, 
            rangekey = rangekey,
            issue_stop = TRUE,
            resolution = resolution,
            message = paste(rangekey,"does not exist! Unable to retrieve",
                "metadata information for non-existent rangekey"))
        RangeKeys <- RangeKeys[RangeKeys %in% rangekey]
    }
    mcol.df <- .prepare_ranges_metadata_mcols(Brick = Brick_filepath, 
        rangekeys = RangeKeys)
    return(mcol.df)
}

#' Fetch the ranges associated to a rangekey or chromosome.
#'
#' `Brick_get_ranges` will get a ranges object if present in the Brick store
#' and return a GRanges object.
#'
#' If a rangekey is present, the ranges will be retrieve and a GRanges
#' constructed. Metadata columns will also be added. If these are rangekeys
#' other than "Bintable", and had been added using Brick_add_ranges the width
#' and Strand columns may appear as metadata columns. These will most likely
#' be artifacts from converting the original ranges object to a data.frame.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @param chr \strong{Optional}.
#' A chr string specifying the chromosome to select from the ranges.
#'
#' @param rangekey \strong{Required}.
#' A string specifying the name of the ranges.
#'
#' @return Returns a GRanges object with the associated metadata columns that
#' may have been present in the Ranges object.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_get_ranges_test")
#' 
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Brick_get_ranges(Brick = My_BrickContainer, chr = "chr2L", 
#' rangekey = "Bintable", resolution = 100000)
#'
Brick_get_ranges = function(Brick = NA, chr = NA, rangekey = NA, 
    resolution = NA){
    Reference.object <- GenomicMatrix$new()
    BrickContainer_class_check(Brick)
    configuration_na_check(resolution, "resolution")
    configuration_length_check(resolution, "resolution", 1)
    configuration_length_check(chr, "chr", 1)
    configuration_na_check(rangekey, "rangekey")
    configuration_length_check(rangekey, "rangekey", 1)
    if(!Brick_rangekey_exists(Brick = Brick, rangekey = rangekey,
        resolution = resolution)){
        stop("rangekey not found!")
    }
    Brick_filepath <- BrickContainer_get_path_to_file(Brick, 
        resolution = resolution)[1]
    Start <- NULL
    Stride <- NULL
    Count <- NULL
    Index <- NULL
    if(!is.na(chr)){
        chromosomes <- ._Brick_Get_Something_(Group.path = Create_Path(
            c(Reference.object$hdf.ranges.root, rangekey)),
        Brick = Brick_filepath,
        Name = Reference.object$hdf.ranges.chr.name,
        return.what = "data")
        if(any(!(chr %in% chromosomes))){
            stop("chr not found!")
        }
        Starts <- ._Brick_Get_Something_(
            Group.path = Create_Path(
                c(Reference.object$hdf.ranges.root, rangekey)),
            Brick = Brick_filepath,
            Name = Reference.object$hdf.ranges.offset.name,
            return.what = "data")
        Lengths <- ._Brick_Get_Something_(Group.path = Create_Path(
            c(Reference.object$hdf.ranges.root, rangekey)),
        Brick = Brick_filepath,
        Name = Reference.object$hdf.ranges.lengths.name,
        return.what = "data")
        Which.one <- chromosomes == chr
        Start <- Starts[Which.one]
        Stride <- 1
        Count <- Lengths[Which.one]
    }
    Dataset <- ._Brick_Get_Something_(Group.path = Create_Path(
        c(Reference.object$hdf.ranges.root, rangekey)),
        Brick = Brick_filepath, 
        Name = Reference.object$hdf.ranges.dataset.name,
        Start = Start, Stride = Stride,
        Count = Count, return.what = "data")
    Dataset <- Brick_make_ranges(chrom = Dataset[,'chr'],
        start = Dataset[,'start'], end = Dataset[,'end'])
    MCols <- Brick_list_ranges_mcols(Brick = Brick, rangekey = rangekey,
        resolution = resolution)
    if(is.data.frame(MCols)){
        Fltr <- MCols$m.col %in%
        Reference.object$genomic.ranges.protected.names
        GRangesCols <- MCols$m.col[Fltr]
        MCols.col <- MCols$m.col[!Fltr]
        m.start <- Start[1]
        m.stride <- Stride[1]
        m.count <- Count[1]
        if(length(GRangesCols) > 0){
            genomic.ranges.FUN.names <- c("strand", "seqlevels", 
                "seqlengths", "width")
            FUN.names <- genomic.ranges.FUN.names[genomic.ranges.FUN.names
            %in% GRangesCols]
            for(FUN_name in FUN.names){
                a_value <- ._Brick_Get_Something_(
                    Group.path = Create_Path(
                    c(Reference.object$hdf.ranges.root, rangekey)),
                Brick = Brick_filepath, Name = FUN_name, Start = m.start, 
                Stride = m.stride, Count = m.count, return.what = "data")
                if(FUN_name == "strand"){
                    strand(Dataset) <- a_value
                }else if(FUN_name == "seqlevels"){
                    seqlevels(Dataset) <- a_value
                }else if(FUN_name == "seqlengths"){
                    seqlengths(Dataset) <- a_value
                }else if(FUN_name == "width"){
                    width(Dataset) <- a_value
                }
            }
        }
        if(length(MCols.col) > 0){
            MCols.DF.list <- lapply(MCols.col,function(x){
                Dataset <- ._Brick_Get_Something_(
                    Group.path = Create_Path(c(
                        Reference.object$hdf.ranges.root,
                        rangekey)), Brick = Brick_filepath,
                    Name = x, Start = m.start, Stride = m.stride,
                    Count = m.count, return.what = "data")
                DF <- DataFrame(Temp = Dataset)
                colnames(DF) <- x
                DF
            })
            MCols.DF <- do.call(cbind,MCols.DF.list)
            mcols(Dataset) <- MCols.DF           
        }
    }
    return(Dataset)
}

#' Returns the binning table associated to the Hi-C experiment.
#'
#' `Brick_get_bintable` makes a call to \code{\link{Brick_get_ranges}} to
#' retrieve the binning table of the associated Brick store. This is equivalent
#' to passing the argument rangekey = "bintable" in
#' \code{\link{Brick_get_ranges}}
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @param chr \strong{Optional}.
#' A chr string specifying the chromosome to select from the ranges.
#'
#' @return Returns a GRanges object containing the binning
#' table associated to the Brick store.
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_get_bintable_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Brick_get_bintable(Brick = My_BrickContainer, resolution = 100000)
#'
#' @seealso Brick_get_ranges
Brick_get_bintable = function(Brick, chr = NA, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    Table <- Brick_get_ranges(Brick = Brick, chr = chr,
        rangekey = Reference.object$hdf.bintable.ranges.group,
        resolution = resolution)
    return(Table)
}

#' Returns the position of the supplied ranges in the binning table
#' associated to the Hi-C experiment.
#'
#' `Brick_fetch_range_index` constructs a ranges object using
#' \code{\link{Brick_make_ranges}}, creates an overlap operation using
#' \code{GenomicRanges::findOverlaps}, where the constructed ranges is
#' the \emph{subject} and the Hi-C experiment associated binning table is the
#' \emph{query}. The return of this object is a list of ranges with their
#' corresponding indices in the binning table.
#' 
#' @inheritParams Brick_get_chrominfo
#'
#' @param chr \strong{Required}.
#' A character vector of length N specifying the chromosomes to select from
#' the ranges.
#'
#' @param start \strong{Required}.
#' A numeric vector of length N specifying the start positions in the
#' chromosome
#'
#' @param end \strong{Required}.
#' A numeric vector of length N specifying the end positions in the chromosome
#'
#' @param names \strong{Optional}.
#' A character vector of length N specifying the names of the chromosomes. If
#' absent, names will take the form chr:start:end.
#'
#' @param type \strong{Optional}. Default any
#' Type of overlap operation to do. It should be one of two, any or within.
#' any considers any overlap (atleast 1 bp) between the provided ranges and
#' the binning table.
#'
#' @return Returns a GenomicRanges object of same length as the chr, start, end
#' vectors provided. The object is returned with an additional column, Indexes.
#' Indexes is a column of class \code{IRanges::IntegerList}, which is
#' part of the larger \code{IRanges::AtomicList} superset. This
#' "Indexes" column can be accessed like a normal GRanges column with the
#' additional list accessor [[]] in place of the normal vector accessor [].
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "fetch_range_index_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#'
#' Chrom <- c("chr2L","chr2L")
#' Start <- c(1,40000)
#' End <- c(1000000,2000000)
#'
#' Test_Run <- Brick_fetch_range_index(Brick = My_BrickContainer, 
#' chr = Chrom, start = Start, end = End, resolution = 100000)
#' Test_Run$Indexes[[1]]
#'
Brick_fetch_range_index = function(Brick = NA, chr = NA, start = NA, 
    end = NA, names = NA, resolution = NA, type = "any"){
    AllTypes<-c("any","within")
    configuration_na_check(chr, "chr")
    configuration_na_check(start, "start")
    configuration_na_check(end, "end")
    configuration_na_check(chr, "Brick")
    if( any(!(type %in% AllTypes)) ){
        stop("type takes one of two arguments: c(\"any\",\"within\")")
    }
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    if(any(!(chr %in% ChromInfo[,'chr']))){
        stop("Provided chr does not exist in chromosome list.\n")
    }
    if(any(!is.character(chr) | !._Check_numeric(start) |
        !._Check_numeric(end))){
        stop("Provided chr, start, end do not match expected class ",
            "definitions of character, numeric, numeric\n")
    }
    Unique.chromosomes <- unique(chr)
    OverlapByChromosome.list <- lapply(Unique.chromosomes, function(cur.chr){
        Filter <- chr==cur.chr
        Cur.Chrom <- chr[Filter]
        Cur.Start <- start[Filter]
        Cur.end <- end[Filter]
        Cur.Names <- names[Filter]
        SubjectRanges <- Brick_get_bintable(Brick = Brick, chr = cur.chr,
            resolution = resolution)
        if( any(!(Cur.end <= max(end(SubjectRanges)) &
            Cur.Start >= min(start(SubjectRanges)))) ){
            stop("Start or end is out of ranges for Bintable\n")
        }
        QueryRanges <- Brick_make_ranges(chrom=Cur.Chrom, start=Cur.Start,
            end=Cur.end, names=Cur.Names)
        elementMetadata(QueryRanges)[["Indexes"]] <- IntegerList(NA)
        HitsObject <- findOverlaps(SubjectRanges,QueryRanges,type=type)
        UniqueQueries <- seq_along(QueryRanges)
        elementMetadata(QueryRanges)[["Indexes"]][UniqueQueries] <- lapply(
            UniqueQueries, function(x){
            A.Query <- UniqueQueries[x]
            MatchingQueries <- queryHits(HitsObject)[
            subjectHits(HitsObject)==A.Query]
            MatchingQueries
            })
        QueryRanges
    })
    OverlapByChromosome <- do.call(c,unlist(OverlapByChromosome.list,
        use.names = FALSE))
    return(OverlapByChromosome)
}

#' Provides the overlapping position (within) from the bintable.
#'
#' `Brick_return_region_position` takes as input a human-readable coordinate
#' format of the form chr:start:end and outputs the overlapping bintable
#' positions. This module does a "within" operation. So only bins which overlap
#' completely with the region will be returned. This is not an iterable module,
#' so the user has to make iterative calls to the module itself.
#'
#' @section Design choice:
#' This may seem to be a poor design choice at first glance, but I do not
#' think this to be the case. By not being iterable, this function circumvents
#' the problem of how to structure the data for the user. If one more element
#' was accepted, the return object would have become a list, which increases
#' the data structure complexity significantly for users who are just starting
#' out with R. Therefore this problem is left for the users themselves to
#' deal with.
#' 
#' @inheritParams Brick_get_chrominfo
#'
#' @param region \strong{Required}.
#' A character vector of length 1 specifying the region to overlap. It must
#' take the form chr:start:end.
#'
#' @return Returns a 1 dimensional vector containing the position of the
#' overlapping regions in the bintable associated the Brick store.
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"),
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "region_position_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Coordinate <- "chr2L:1:1000000"
#' 
#' Test_Run <- Brick_return_region_position(Brick = My_BrickContainer,
#' region = Coordinate, resolution = 100000)
#'
Brick_return_region_position = function(Brick, region, resolution = NA){
    if(!is.character(region) | length(region) > 1){
        stop("region must be a character vector of length 1")
    }
    Coord.Split<- Split_genomic_coordinates(Coordinate=region)
    region.chr<-Coord.Split[[1]][1]
    region.start<-as.numeric(Coord.Split[[1]][2])
    region.stop<-as.numeric(Coord.Split[[1]][3])
    region.ranges<-Brick_fetch_range_index(Brick = Brick, chr=region.chr,
        start=region.start, end=region.stop, resolution = resolution,
        type="within")
    Vector.coordinates <- region.ranges$Indexes[[1]]
    return(Vector.coordinates)
}

#' Load a NxM dimensional matrix into the Brick store.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @param chr1 \strong{Required}.
#' A character vector of length 1 specifying the chromosome corresponding to
#' the rows of the matrix
#'
#' @param chr2 \strong{Required}.
#' A character vector of length 1 specifying the chromosome corresponding to
#' the columns of the matrix
#'
#' @param matrix_file \strong{Required}.
#' A character vector of length 1 specifying the name of the file to load as a
#' matrix into the Brick store.
#'
#' @param delim \strong{Optional}. Default " "
#' The delimiter of the matrix file.
#'
#' @param remove_prior \strong{Optional}. Default FALSE
#' If a matrix was loaded before, it will not be replaced. Use remove_prior to
#' override and replace the existing matrix.
#'
#' @param num_rows \strong{Optional}. Default 2000
#' Number of rows to read, in each chunk.
#'
#' @param is_sparse \strong{Optional}. Default FALSE
#' If true, designates the matrix as being a sparse matrix, and computes the
#' sparsity.index. The sparsity index measures the proportion of non-zero rows
#' or columns at a certain distance from the diagonal (100) in cis interaction
#' matrices.
#'
#' @param sparsity_bins \strong{Optional}. Default 100
#' With regards to computing the sparsity.index, this parameter decides the
#' number of bins to scan from the diagonal.
#'
#' @return Returns TRUE if all went well.
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_load_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#'
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#'
Brick_load_matrix = function(Brick = NA, chr1 = NA, chr2 = NA, resolution = NA,
    matrix_file = NA, delim = " ", remove_prior = FALSE, num_rows = 2000, 
    is_sparse = FALSE, sparsity_bins = 100){
    Reference.object <- GenomicMatrix$new()
    BrickContainer_class_check(Brick)
    ListVars <- list(chr1 = chr1, chr2 = chr2, 
        resolution = resolution, matrix_file = matrix_file, 
        delim = delim, remove_prior = remove_prior, 
        num_rows = num_rows, is_sparse = is_sparse, 
        sparsity_bins = sparsity_bins)

    lapply(seq_along(ListVars), function(x){
        configuration_length_check(ListVars[[x]], names(ListVars[x]), 1)
    })
    lapply(seq_along(ListVars), function(x){
        configuration_na_check(ListVars[[x]], names(ListVars[x]))
    })
    Brick_filepath <- BrickContainer_get_path_to_file(Brick, 
        chr1 = chr1, chr2 = chr2, resolution = resolution)
    if(length(Brick_filepath) == 0){
        stop("Did not find a Brick file corresponding to", chr1, chr2, 
            " pair. Note that, chr2 position is expected to be greater",
            " than or equal to than chr1 position in the chromosome list.",
            " Please check the chromosome list using",
            " BrickContainer_list_chromosomes")
    }
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    if(Brick_matrix_isdone(Brick = Brick, chr1 = chr1,
        chr2 = chr2, resolution = resolution) && !remove_prior){
        stop("A matrix was preloaded before. ",
            "Use remove_prior = TRUE to force value replacement\n")
    }
    Chrom.info.df <- Brick_get_chrominfo(Brick = Brick, 
        resolution = resolution)
    chr1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr1,"nrow"]
    chr2.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr2,"nrow"]
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,chr1,chr2))
    compute_sparsity <- FALSE
    if(is_sparse && chr1 == chr2){
        compute_sparsity <- TRUE
    }
    RetVar <- ._ProcessMatrix_(Brick = Brick_filepath, 
        Matrix.file = matrix_file, delim = delim, Group.path = Group.path, 
        chr1.len = chr1.len, chr2.len = chr2.len, num.rows = num_rows, 
        is.sparse = is_sparse, compute.sparsity = compute_sparsity, 
        sparsity.bins = sparsity_bins)
    return(RetVar)
}


#' Load a NxN dimensional sub-distance \emph{cis} matrix into
#' the Brick store.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @param chr \strong{Required}.
#' A character vector of length 1 specifying the chromosome corresponding to
#' the rows and cols of the matrix
#'
#' @param distance \strong{Required}. Default NULL.
#' For very high-resolution matrices, read times can become extremely slow and
#' it does not make sense to load the entire matrix into the data structure, as
#' after a certain distance, the matrix will become extremely sparse. This
#' ensures that only interactions upto a certain distance from the main
#' diagonal will be loaded into the data structure.
#'
#' @param num_rows \strong{Optional}. Default 2000
#' Number of rows to insert per write operation in the HDF file.
#'
#' @return Returns TRUE if all went well.
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_load_dist_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#'
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_cis_matrix_till_distance(Brick = My_BrickContainer, 
#' chr = "chr2L", resolution = 100000, matrix_file = Matrix_file, 
#' delim = " ", distance = 30, remove_prior = TRUE)
#'
Brick_load_cis_matrix_till_distance = function(Brick = NA, chr = NA, 
    resolution = NA, matrix_file, delim = " ", distance, remove_prior = FALSE,
    num_rows = 2000, is_sparse = FALSE, sparsity_bins = 100){

    Reference.object <- GenomicMatrix$new()
    ListVars <- list(chr = chr, resolution = resolution, 
        matrix_file = matrix_file, delim = delim, distance = distance, 
        remove_prior = remove_prior, num_rows = num_rows, 
        is_sparse = is_sparse, sparsity_bins = sparsity_bins)
    lapply(seq_along(ListVars), function(x){
        configuration_length_check(ListVars[[x]], names(ListVars[x]), 1)
    })
    lapply(seq_along(ListVars), function(x){
        configuration_na_check(ListVars[[x]], names(ListVars[x]))
    })
    Brick_filepath <- BrickContainer_get_path_to_file(Brick, 
        chr1 = chr, chr2 = chr, resolution = resolution)
    if(length(Brick_filepath) == 0){
        stop("Did not find a Brick file corresponding to", chr, chr, 
            " pair. Note that, chr2 position is expected to be greater",
            " than or equal to than chr1 position in the chromosome list.",
            " Please check the chromosome list using",
            " BrickContainer_list_chromosomes")
    }
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr,
        chr2 = chr, resolution = resolution)){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    if(Brick_matrix_isdone(Brick = Brick, chr1 = chr,
        chr2 = chr, resolution = resolution) && !remove_prior){
        stop("A matrix was preloaded before. Use remove_prior = TRUE to ",
            "force value replacement\n")
    }
    Chrom.info.df <- Brick_get_chrominfo(Brick = Brick, 
        resolution = resolution)
    chr1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr,"nrow"]
    chr2.len <- chr1.len
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,chr,chr))
    compute_sparsity <- FALSE
    if(is_sparse){
        compute_sparsity <- TRUE
    }
    RetVar <- ._Process_matrix_by_distance(Brick = Brick_filepath,
        Matrix.file = matrix_file, delim = delim, Group.path = Group.path,
        chr1.len = chr1.len, chr2.len = chr2.len, num.rows = num_rows,
        distance = distance, is.sparse = is_sparse,
        compute.sparsity = compute_sparsity, sparsity.bins = sparsity_bins)
    return(RetVar)
}


#' Load a NxN dimensional matrix into the Brick store from an mcool file.
#'
#' Read an mcool contact matrix coming out of 4D nucleome projects into a
#' Brick store.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @inheritParams Create_many_Bricks_from_mcool
#' 
#' @param mcool \strong{Required}. Path to an mcool file.
#'
#' @param norm_factor \strong{Optional}. Default "Iterative-Correction".
#' The normalization factor to use for normalization from an mcool file.
#' norm_factor currently accepts one of "Iterative-Correction", "Knight-Ruitz",
#' "Vanilla-coverage", "Vanilla-coverage-square-root" and NULL. If NULL,
#' the function will load only counts from the mcool file.
#'
#' @param matrix_chunk \strong{Optional}. Default 2000.
#' The nxn matrix square to fill per iteration in a mcool file.
#'
#' @param cooler_read_limit \strong{Optional}. Default 10000000.
#' cooler_read_limit sets the upper limit for the number of records per matrix
#' chunk. If the number of records per chunk is higher than this value, the 
#' matrix_chunk value will be re-evaluated dynamically.
#' 
#' @return Returns TRUE if all went well.
#'
#' @examples
#'
#' \dontrun{
#'
#' require(curl)
#' out_dir <- file.path(tempdir(),"mcool_load_test")
#' dir.create(path = out_dir)
#' curl_download(url = paste("https://data.4dnucleome.org/",
#' "files-processed/4DNFI7JNCNFB/",
#' "@@download/4DNFI7JNCNFB.mcool", sep = ""),
#' destfile = file.path(out_dir,"H1-hESC-HiC-4DNFI7JNCNFB.mcool"))
#'
#' mcool <- file.path(out_dir,"H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' 
#' My_BrickContainer <- Create_many_Bricks_from_mcool(
#' output_directory = out_dir,
#' file_prefix = "Test",
#' mcool = mcool,
#' resolution = 50000,
#' experiment_name = "A random 4DN dataset")
#'
#' Brick_load_data_from_mcool(Brick = My_BrickContainer, mcool = mcool,
#' resolution = 50000, matrix_chunk = 2000, remove_prior = TRUE,
#' norm_factor = "Iterative-Correction")
#'
#' }
#'
#'
#' @seealso \code{\link{Create_many_Bricks_from_mcool}} to create matrix from 
#' an mcool file, \code{\link{Brick_list_mcool_resolutions}} to list available
#' resolutions in an mcool file, \code{\link{Brick_list_mcool_normalisations}}
#' to list available normalisation factors in the mcool file.
#'
Brick_load_data_from_mcool <- function(Brick, mcool, resolution = NULL, 
    matrix_chunk = 2000, cooler_read_limit = 10000000, remove_prior = FALSE,
    norm_factor = "Iterative-Correction"){
    Reference.object <- GenomicMatrix$new()
    
    resolutions <- Brick_list_mcool_resolutions(mcool = mcool)
    if(!is.null(resolutions) & is.null(resolution)){
        stop("resolution must be provided when",
            " different resolutions are present in an mcool file.\n")
    }
    resolution <- .format_resolution(resolution)
    if(!is.null(resolutions) & !is.null(resolution)){
        if(!(resolution %in% resolutions)){
            stop("resolution not found in mcool file.",
                " Please check available resolutions",
                " with Brick_list_mcool_resolutions.\n")
        }
    }
    Norm_factor <- NULL
    if(!is.null(norm_factor)){
        All.factors <- Brick_list_mcool_normalisations(names.only=TRUE)
        if(!Brick_mcool_normalisation_exists(mcool = mcool,
            norm_factor = norm_factor, resolution = resolution)){
            Factors.present <- All.factors[vapply(All.factors, 
                function(x){
                    Brick_mcool_normalisation_exists(mcool = mcool,
                    norm_factor = x, resolution = resolution)
                }, TRUE)]
            stop(norm_factor," was not found in this mcool file.\n",
                "Normalisation factors present in the mcool file are: ",
                paste(Factors.present, collapse = ", "))
        }
        Norm_factors <- Brick_list_mcool_normalisations()
        Norm_factor <- Norm_factors[norm_factor]
        names(Norm_factor) <- NULL
    }
    RetVar <- .process_mcool(Brick = Brick, mcool_path = mcool, 
        resolution = resolution, matrix_chunk = matrix_chunk, 
        norm_factor = Norm_factor, cooler_read_limit = cooler_read_limit,
        has_resolution = !is.null(resolutions))
    return(RetVar)
}

#' Check if a matrix has been loaded for a chromosome pair.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns a logical vector of length 1, specifying if a matrix has
#' been loaded or not.
#'
#' @examples 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_isdone_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_isdone(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#'
Brick_matrix_isdone = function(Brick, chr1, chr2, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Matrix.list <- Brick_list_matrices(Brick = Brick, chr1 = chr1, 
        chr2 = chr2, resolution = resolution)
    return(Matrix.list[Matrix.list$chr1 == chr1 &
        Matrix.list$chr2 == chr2, "done"])
}

#' Check if a matrix for a chromosome pair is sparse.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns a logical vector of length 1, specifying if a matrix was
#' loaded as a sparse matrix.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_issparse_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_issparse(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#'
Brick_matrix_issparse = function(Brick, chr1, chr2, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Matrix.list <- Brick_list_matrices(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)
    return(Matrix.list[Matrix.list$chr1 == chr1 &
        Matrix.list$chr2 == chr2, "sparsity"])
}


#' Get the maximum loaded distance from the diagonal of any matrix.
#'
#' If values beyond a certain distance were not loaded in the matrix, this
#' distance parameter is useful. This package by default will check this param
#' to make sure that it is not returning non-existent data.
#' 
#' `Brick_matrix_maxdist` will return this parameter.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns an integer vector of length 1, specifying the maximum
#' distance loaded for that matrix
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_maxdist_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_maxdist(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#'
Brick_matrix_maxdist = function(Brick, chr1, chr2, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop("chr1 chr2 pairs were not found\n")
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)){
        stop("chr1 chr2 pairs were not loaded\n")
    }
    Matrix.list <- Brick_list_matrices(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)
    return((Matrix.list[Matrix.list$chr1 == chr1 &
        Matrix.list$chr2 == chr2, "distance"]))
}

#' Check if a chromosome pair exists.
#'
#' Matrices are created when the bintable is loaded and the chromosome names
#' are provided. If a user is in doubt regarding whether a matrix is present or
#' not it is useful to check this function. If the Bintable did not contain a
#' particular chromosome, any matrices for that chromosome would not be present
#' in the file
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns a logical vector of length 1, specifying if the matrix
#' exists or not.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_exists_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_exists(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#' 
Brick_matrix_exists = function(Brick, chr1, chr2, resolution = NA){
    configuration_na_check(resolution, "resolution")
    Brick_list <- BrickContainer_list_files(Brick = Brick, 
        chr1 = chr1, chr2 = chr2, resolution = resolution)
    return(nrow(Brick_list) > 0)
}

#' Return the value range of the matrix
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns a numeric vector of length 2, specifying the minimum and
#' maximum finite real values in the matrix.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_minmax_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_minmax(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#' 
Brick_matrix_minmax = function(Brick, chr1, chr2, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Matrix.list <- Brick_list_matrices(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- c(Matrix.list[Filter, "min"],Matrix.list[Filter, "max"])
    return(Extent)
}

#' Return the dimensions of a matrix
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns the dimensions of a Hi-C matrix for any given
#' chromosome pair.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_dimension_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_dimensions(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#' 
Brick_matrix_dimensions = function(Brick, chr1, chr2, resolution = NA){
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Reference.object <- GenomicMatrix$new()
    Brick_filepath <- BrickContainer_get_path_to_file(Brick = Brick, 
        chr1 = chr1, chr2 = chr2, resolution = resolution)
    Extents <- ._GetDimensions(group.path = Create_Path(
        c(Reference.object$hdf.matrices.root,chr1,chr2)),
        dataset.path = Reference.object$hdf.matrix.name,
        File = Brick_filepath, return.what = "size")
    return(Extents)
}


#' Return the filename of the loaded matrix
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns a character vector of length 1 specifying the filename of
#' the currently loaded matrix.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "matrix_filename_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_matrix_filename(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000)
#' 
Brick_matrix_filename = function(Brick, chr1, chr2, resolution = NA){
    Reference.object <- GenomicMatrix$new()
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Matrix.list <- Brick_list_matrices(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- Matrix.list[Filter, "filename"]
    return(Extent)
}

#' Return values separated by a certain distance.
#'
#' `Brick_get_values_by_distance` can fetch values with or without
#' transformation or subsetted by a certain distance. Please note,
#' this module is not an iterable module.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @param chr \strong{Required}.
#' A string specifying the chromosome for the cis Hi-C matrix from which values
#' will be retrieved at a certain distance.
#'
#' @param distance \strong{Required}. 0 based.
#' Fetch values separated by distance.
#' 
#' @param constrain_region \strong{Optional}.
#' A character vector of length 1 with the form chr:start:end specifying the
#' region for which the distance values must be retrieved.
#'
#' @param batch_size \strong{Optional}. Default 500
#' A numeric vector of length 1 specifying the size of the chunk to retrieve
#' for diagonal selection.
#' 
#' @param FUN \strong{Optional}.
#' If provided a data transformation with FUN will be applied before values
#' are returned.
#' 
#' @return Returns a numeric vector of length N depending on the presence of
#' constrain_region, FUN and distance from the main diagonal.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "val_by_dist_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_get_values_by_distance(Brick = My_BrickContainer, chr = "chr2L",
#' distance = 0, resolution = 100000)
#'
#' Failsafe_median <- function(x){
#'      x[is.nan(x) | is.infinite(x) | is.na(x)] <- 0
#'      return(median(x))
#' }
#'
#' Brick_get_values_by_distance(Brick = My_BrickContainer, chr = "chr2L", 
#' resolution = 100000, distance = 4, FUN = Failsafe_median)
#'
#' @seealso \code{\link{Brick_get_matrix_within_coords}} to get matrix by
#' using matrix coordinates, \code{\link{Brick_fetch_row_vector}} to get values
#' in a certain row/col and subset them, \code{\link{Brick_get_vector_values}}
#' to get values using matrix coordinates, \code{\link{Brick_get_matrix}} to
#' get matrix by using matrix coordinates.
#'
Brick_get_values_by_distance = function(Brick, chr, distance, 
    resolution, constrain_region = NULL, batch_size=500, FUN=NULL){
    if(any(vapply(list(Brick,chr,distance),is.null,TRUE))) {
        stop("Brick, chr, distance cannot be NULL.\n")
    }
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr, 
        chr2 = chr, resolution = resolution) &
        !Brick_matrix_isdone(Brick = Brick, chr1 = chr, 
            chr2 = chr, resolution = resolution)){
        stop("Chromosome is not listed or has not been ",
            "loaded in this HDF file.\n")
    }
    if(any(vapply(list(Brick,chr,distance),length,1) > 1)){
        stop("Brick, chr and distance can only have values of length 1.\n")
    }
    Nrow <- (ChromInfo[ChromInfo$chr == chr,"nrow"]-1)
    if(any(distance > Nrow, distance < 0)){
        stop("distance must range between 0 and",Nrow,"\n")
    }
    Max.dist <- Brick_matrix_maxdist(Brick = Brick, chr1 = chr, chr2 = chr,
        resolution = resolution)
    if(distance > Max.dist){
        stop(paste("The farthest pixel loaded for ",
            "this matrix was at a distance of "
            ,Max.dist,"bins from the diagonal. ",
            " The current selection subsets out-of-bounds data.\n"))
    }
    Root.folders <- Reference.object$GetRootFolders()
    Path <- Create_Path(c(Root.folders['matrices'],chr,chr))
    Vector.start <- 1
    Vector.stop <- ChromInfo[ChromInfo$chr==chr,"nrow"]
    if(!is.null(constrain_region)){
        Vector.coordinates <- Brick_return_region_position(Brick = Brick,
            region=constrain_region, resolution = resolution)
        if(is.null(Vector.coordinates)){
            stop("Overlap operation was unsuccessful! ",
                "Please check coordinates "
                ,constrain_region)
        }
        Vector.start <- min(Vector.coordinates)
        Vector.stop <- max(Vector.coordinates)
    }
    Starting.col <- Vector.start + distance
    Count <- ((Vector.stop - Vector.start) + 1) - distance
    # tot.mat.extracted <- Count * Count
    # cat("Boo ",Count,"\n")
    Start <- c(Vector.start,Starting.col)
    Stride <- c(1,1)
    Counts <- Count
    CumSums <- 0
    Groups <- c(Reference.object$hdf.matrices.root,chr)
    if(Count > batch_size){
        repeated <- floor(Count/batch_size)
        Counts <- rep(batch_size,repeated)
        if(repeated != ceiling(Count/batch_size)){
            cumulative <- sum(Counts)
            Counts <- c(Counts,(Count-cumulative))
        }
        CumSums <- cumsum(c(0,Counts[seq_len(length(Counts)-1)]))
    }
    Brick_filepath <- BrickContainer_get_path_to_file(Brick = Brick, 
    chr1 = chr, chr2 = chr, resolution = resolution)
    DistancesVector.list <- lapply(seq_len(length(Counts)),function(x){
        Count <- Counts[x]
        Offset <- CumSums[x]
        cur.start <- Start+Offset
        diag(._Brick_Get_Something_(Group.path = Path, Brick = Brick_filepath,
            Name = Reference.object$hdf.matrix.name,
            Start= cur.start, Stride = Stride, Count = c(Count,Count),
            return.what = "data"))
    })

    DistancesVector <- do.call(c,DistancesVector.list)
    if(is.null(FUN)){
        return(DistancesVector)
    }else{
        return(FUN(DistancesVector))
    }
    # Preparing Start stride and
}

#' Return a matrix subset between two regions.
#'
#' `Brick_get_matrix_within_coords` will fetch a matrix subset after
#' creating an overlap operation between both regions and the bintable
#' associated to the Brick store.
#' This function calls \code{\link{Brick_get_matrix}}.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @param x_coords \strong{Required}.
#' A string specifying the region to subset on the rows. It takes the form
#' chr:start:end. An overlap operation with the associated bintable will be
#' done to identify the bins to subset on the row
#'
#' @param y_coords \strong{Required}.
#' A string specifying the region to subset on the rows. It takes the form
#' chr:start:end. An overlap operation with the associated bintable will be
#' done to identify the bins to subset on the column
#'
#' @param force \strong{Optional}. Default FALSE
#' If true, will force the retrieval operation when matrix contains loaded
#' data until a certain distance.
#'
#' @param FUN \strong{Optional}.
#' If provided a data transformation with FUN will be applied before
#' the matrix is returned.
#'
#' @return Returns a matrix of dimension x_coords binned length by y_coords
#' binned length. This may differ based on FUN.
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_matrix_coords_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_get_matrix_within_coords(Brick = My_BrickContainer,
#' x_coords = "chr2L:1:1000000",
#' y_coords = "chr2L:1:1000000",
#' resolution = 100000)
#'
#' Brick_get_matrix_within_coords(Brick = My_BrickContainer,
#' x_coords = "chr2L:1:1000000",
#' y_coords = "chr2L:1:1000000",
#' resolution = 100000,
#' FUN = mean)
#'
#' 
#' @seealso \code{\link{Brick_get_matrix}} to get matrix by using matrix
#' coordinates, \code{\link{Brick_get_values_by_distance}} to get values
#' separated at a certain distance, \code{\link{Brick_fetch_row_vector}}
#' to get values in a certain row/col and subset them,
#' \code{\link{Brick_get_vector_values}} to get values using matrix
#' coordinates.
Brick_get_matrix_within_coords = function(Brick, x_coords,
    y_coords, resolution, force = FALSE, FUN=NULL){
    type <- "within"
    if( (is.null(x_coords)) | (is.null(y_coords)) ){
        stop("x_coords, y_coords and Brick cannot be NULL")
    }
    if(!(length(x_coords)==1) | !(length(y_coords)==1)){
        stop("This function processes single process calls at a time.",
            " Setup an Iterator for more functionality")
    }
    if( !is.character(x_coords) | !is.character(y_coords) ){
        stop("Two string variables were expected for x_coords & y_coords,",
            " found x_coords class ", class(x_coords), " and y_coords class ",
            class(y_coords))
    }
    xcoords.split <- Split_genomic_coordinates(Coordinate=x_coords)
    chr1 <- xcoords.split[[1]][1]
    chr1.start <- as.numeric(xcoords.split[[1]][2])
    chr1.stop <- as.numeric(xcoords.split[[1]][3])
    chr2.split <- Split_genomic_coordinates(Coordinate=y_coords)
    chr2 <- chr2.split[[1]][1]
    chr2.start <- as.numeric(chr2.split[[1]][2])
    chr2.stop <- as.numeric(chr2.split[[1]][3])
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    ChromosomeList <- ChromInfo[,"chr"]
    if( any(!(c(chr1,chr2) %in% ChromosomeList)) ){
        stop("Provided chromosomes were not found in chromosome list.")
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop(chr1," ",chr2," matrix is yet to be loaded into the class.")
    }
    chr1.ranges <- Brick_fetch_range_index(Brick = Brick, chr = chr1,
        start = chr1.start, end = chr1.stop, names=NA, 
        resolution = resolution, type=type)
    chr2.ranges <- Brick_fetch_range_index(Brick = Brick, chr = chr2,
        start = chr2.start, end = chr2.stop, names=NA, 
        resolution = resolution, type=type)
    if(is.null(chr1.ranges$Indexes[[1]]) | is.null(chr2.ranges$Indexes[[1]])){
        stop("Overlap operation was unsuccessful! Please check coordinates ",
            x_coords," & ",y_coords)
    }
    x_vector <- chr1.ranges$Indexes[[1]]
    y_vector <- chr2.ranges$Indexes[[1]]
    Matrix <- Brick_get_matrix(Brick = Brick, chr1 = chr1, chr2 = chr2,
        x_coords = x_vector, y_coords = y_vector, resolution = resolution, 
        force = force, FUN = FUN)
    return(Matrix)
}

#' Return a matrix subset.
#'
#' `Brick_get_matrix` will fetch a matrix subset between row values
#' ranging from min(x_coords) to max(x_coords) and column values ranging from
#' min(x_coords) to max(x_coords)
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @param x_coords \strong{Required}.
#' A one-dimensional numeric vector specifying the rows to subset.
#'
#' @param y_coords \strong{Required}.
#' A one-dimensional numeric vector specifying the columns to subset.
#'
#' @param FUN \strong{Optional}.
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#'
#' @inheritParams Brick_get_matrix_within_coords
#'
#' @return Returns a matrix of dimension x_coords length by y_coords length.
#' This may differ based on the operations with FUN.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_matrix_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_get_matrix(Brick = My_BrickContainer, chr1 = "chr2L", chr2 = "chr2L",
#' x_coords = c(1:10), y_coords = c(1:10), resolution = 100000)
#'
#' @seealso \code{\link{Brick_get_matrix_within_coords}} to get matrix by using
#' matrix genomic coordinates, \code{\link{Brick_get_values_by_distance}} to
#' get values separated at a certain distance,
#' \code{\link{Brick_fetch_row_vector}} to getvalues in a certain row/col and
#' subset them, \code{\link{Brick_get_vector_values}} to get values using
#' matrix coordinates.
Brick_get_matrix = function(Brick, chr1, chr2, x_coords,
    y_coords, resolution, force = FALSE, FUN=NULL){
    # cat(" Rows: ",x_coords," Cols: ",y_coords,"\n")
    if(any(!._Check_numeric(x_coords) | !._Check_numeric(y_coords))){
        stop("x_coords and y_coords must be numeric.\n")
    }
    if(is.null(chr1) | is.null(chr2) | is.null(x_coords) | is.null(y_coords)){
        stop("Either of chr1, chr2, x_coords or y_coords ",
            "were provided as NULL values.\n")
    }
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    ChromosomeList <- ChromInfo[,"chr"]
    if( any(!(c(chr1,chr2) %in% ChromosomeList)) ){
        stop("Provided chromosomes were not found in chromosome list.\n")
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)){
        stop(chr1,chr2," matrix is yet to be loaded into the class.\n")
    }
    chr1.len <- ChromInfo$nrow[ChromInfo$chr == chr1]
    chr2.len <- ChromInfo$nrow[ChromInfo$chr == chr2]
    if(any(x_coords > chr1.len) |
        any(y_coords > chr2.len) |
        min(x_coords,y_coords) < 1 ) {
        stop("x_coords or y_coords falls outside ",
            "the bounds of loaded Bintables")
    }
    Matrix <- Brick_get_vector_values(Brick = Brick, chr1=chr1, chr2=chr2,
        xaxis=x_coords, yaxis=y_coords, resolution = resolution, force = force)
    if(is.null(FUN)){
        return(Matrix)             
    }else{
        return(FUN(Matrix))
    }
}

#' Return row or col vectors.
#'
#' `Brick_fetch_row_vector` will fetch any given rows from a matrix. If
#' required, the rows can be subsetted on the columns and transformations
#' applied. Vice versa is also true, wherein columns can be retrieved and
#' rows subsetted.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @param by \strong{Required}.
#' One of two possible values, "position" or "ranges". A one-dimensional
#' numeric vector of length 1 specifying one of either position or ranges.
#'
#' @param vector \strong{Required}.
#' If by is position, a 1 dimensional numeric vector containing the rows to be
#' extracted is expected. If by is ranges, a 1 dimensional character vector
#' containing the names of the bintable is expected.
#' This function does not do overlaps. Rather it returns any given row or
#' column based on their position or names in the bintable.
#'
#' @param regions \strong{Optional}. Default NULL
#' A character vector of length vector is expected. Each element must be of the
#' form chr:start:end. These regions will be converted back to their original
#' positions and the corresponding rows will be subsetted by the corresponding
#' region element. If the length of regions does not match, the subset
#' operation will not be done and all elements from the rows will be returned.
#'
#' @inheritParams Brick_get_matrix_within_coords
#' 
#' @param flip \strong{Optional}. Default FALSE
#' If present, will flip everything. This is equivalent to selecting columns,
#' and subsetting on the rows.
#'
#' @param FUN \strong{Optional}. Default NULL
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#'
#' @return Returns a list of length vector. Each list element will be of length
#' chr2 binned length or if regions is present the corresponding region length.
#' This may differ based on the operations with FUN.
#'
#' @examples
#'
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_row_vector_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Coordinate <- c("chr2L:1:100000","chr2L:100001:200000")
#' Test_Run <- Brick_fetch_row_vector(Brick = My_BrickContainer,
#' chr1 = "chr2L", chr2 = "chr2L", resolution = 100000, 
#' by = "ranges", vector = Coordinate,
#' regions = c("chr2L:1:1000000", "chr2L:40001:2000000"))
#'
#' @seealso \code{\link{Brick_get_matrix_within_coords}} to get matrix by
#' using matrix genomic coordinates, \code{\link{Brick_get_values_by_distance}}
#' to get values separated at a certain distance,
#' \code{\link{Brick_fetch_row_vector}} to get values in a certain row/col and
#' subset them, \code{\link{Brick_get_matrix}} to get matrix by using
#' matrix coordinates.
#'
Brick_fetch_row_vector = function(Brick, chr1, chr2, resolution, 
    by=c("position", "ranges"), vector, regions=NULL, force = FALSE, 
    flip = FALSE, FUN=NULL){
    Chrom.all <- c(chr1,chr2)
    if(is.null(chr1) | is.null(chr2) | is.null(by) | is.null(vector)) {
        stop("Either of chr1, chr2, by or vector were provided as NULL values")
    }
    if(!(by %in% c("position","ranges")) | length(by) != 1){
        stop("by expects a vector of type character, length 1 and ",
            "takes either one of position or ranges as values")
    }
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    ChromosomeList <- ChromInfo[,"chr"]
    max.dist <- Brick_matrix_maxdist(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)
    if(!is.character(chr1) | !is.character(chr2)){
        stop("Provided Chromosomes does not appear to be of class character")
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop(chr1,chr2," matrix is yet to be loaded.")
    }
    if(by=="position"){
        Class.type <- ._Check_numeric
        Class.exp <- c("numeric","integer")
    }
    if(by=="ranges"){
        Class.type <- is.character
        Class.exp <- "character"
    }
    if(!Class.type(vector)){
        stop("vector must be of class ",
            ifelse(length(Class.exp)>1,paste(Class.exp,collapse=" or "),
                paste(Class.exp))," when by has value ",by)
    }
    if(length(Chrom.all)>2){
        stop("This module is not iterable")
    }
    chr1.ranges <- Brick_get_bintable(Brick = Brick, chr=chr1, 
        resolution = resolution)
    chr2.ranges <- Brick_get_bintable(Brick = Brick, chr=chr2,
        resolution = resolution)
    if(flip){
        chr1.ranges <- Brick_get_bintable(Brick = Brick, chr=chr2,
            resolution = resolution)
        chr2.ranges <- Brick_get_bintable(Brick = Brick, chr=chr1,
            resolution = resolution)
    }
    if(by=="ranges"){
        if(any(!(vector %in% names(chr1.ranges)))){
            stop("All ranges not found in ",chr1," Bintable")
        }
        Positions<-which(names(chr1.ranges) %in% vector)
        names(Positions) <- vector
    }else{
        if(any(max(vector)>max(length(chr1.ranges)) | min(vector)<1)) {
            stop("Position vector falls outside the bounds of ",
                chr1," Bintable")
        }              
        Positions <- vector
        names(Positions) <- names(chr1.ranges[vector])
    }
    PikaPika<-FALSE
    if(!is.null(regions) & length(regions)==length(vector)){
        PikaPika<-TRUE
    }
    Seq.indices <- seq_along(Positions)
    names(Seq.indices) <- names(Positions)
    Vector.values<-lapply(Seq.indices,function(ind){
        x <- Positions[ind]
        if(PikaPika){
            region <- regions[ind]
            y <- Brick_return_region_position(Brick = Brick, region = region,
                resolution = resolution)
        }else{
            y <- seq_len(length(chr2.ranges))
        }
        chr2.names <- names(chr2.ranges[y])
        if(flip){
            chr2.f <- chr1
            chr1 <- chr2
            chr2 <- chr2.f
            x1 <- y
            y <- x
            x <- x1
        }
        Values <- Brick_get_vector_values(Brick = Brick, chr1=chr1, chr2=chr2, 
            resolution = resolution, xaxis=x, yaxis=y, FUN=FUN, force = force)
        return(Values)
    })
    return(Vector.values)
}

#' Return a N dimensional vector selection.
#'
#' `Brick_get_vector_values` is the base function being used by all
#' other matrix retrieval functions.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @param xaxis \strong{Required}.
#' A 1 dimensional vector containing the rows to retrieve. Gaps in this vector
#' may result in unexpected behaviour as the values which are considered are
#' min(xaxis) and max(xaxis) for retrieval.
#'
#' @param yaxis \strong{Required}.
#' A 1 dimensional vector containing the columns to retrieve. Gaps in this
#' vector may result in unexpected behaviour as the values which are considered
#' are min(yaxis) and max(yaxis) for retrieval.
#'
#' @param FUN \strong{Optional}. Default NULL
#' If provided a data transformation with FUN will be applied before the vector
#' is returned.
#'
#' @inheritParams Brick_get_matrix_within_coords
#'
#' @return Returns a vector of length yaxis if length of xaxis is 1. Else
#' returns a matrix of dimension xaxis length by yaxis length.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_vector_val_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_get_vector_values(Brick = My_BrickContainer, chr1 = "chr2L",
#' chr2 = "chr2L", resolution = 100000, xaxis = c(1:10), yaxis = c(1:10))
#'
#' @section Note: Whatever the length of xaxis or yaxis may be, the coordinates
#' under consideration will range from min(xaxis) to max(xaxis) on the rows or
#' min(yaxis) to max(yaxis) on the columns.
#'
Brick_get_vector_values = function(Brick, chr1, chr2, resolution,
    xaxis, yaxis, FUN=NULL, force = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(is.null(chr1) | is.null(chr2)){
        stop("chr1 and chr2 keys cannot be empty!")
    }
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    ChromosomeList <- ChromInfo[,"chr"]
    if(!is.character(chr1) | !is.character(chr2)){
        stop("Provided Chromosomes does not appear to be of class character")
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop(chr1,chr2," matrix is yet to be loaded.")
    }
    if(is.null(xaxis) & is.null(yaxis)){
        stop("Both xaxis and yaxis cannot be null")
    }
    Start <- c(min(xaxis),min(yaxis))
    Stride <- c(1,1)
    Count <- c(length(xaxis),length(yaxis))
    Max.dist <- Brick_matrix_maxdist(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)
    if(max(min(xaxis)-max(yaxis),(max(xaxis) - min(yaxis))) > Max.dist &
    !force){
        stop(paste("The farthest pixel loaded for ",
            "this matrix was at a distance of ",
            Max.dist,"bins from the diagonal. ",
            "The current selection subsets out-of-bounds data.\n"))
    }
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,
        chr1, chr2))
    Brick_filepath <- BrickContainer_get_path_to_file(Brick = Brick, 
        chr1 = chr1, chr2 = chr2, resolution = resolution)
    Vector <- ._Brick_Get_Something_(Group.path = Group.path, 
        Brick = Brick_filepath, Name = Reference.object$hdf.matrix.name, 
        Start = Start, Stride = Stride, Count = Count, return.what = "data")
    if(is.null(FUN)){
        return(Vector)
    }else{
        return(FUN(Vector))
    }
}

#' Return an entire matrix for provided chromosome pair for a resolution.
#'
#' `Brick_get_entire_matrix` will return the entire matrix for the entire 
#' chromosome pair provided an object of class BrickContainer, and values for 
#' chr1, chr2 and resolution values.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @return Returns an object of class matrix with dimensions corresponding to
#' chr1 binned length by chr2 binned length.
#'
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_vector_val_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Entire_matrix <- Brick_get_entire_matrix(Brick = My_BrickContainer, 
#' chr1 = "chr2L", chr2 = "chr2L", resolution = 100000)
#'
Brick_get_entire_matrix = function(Brick, chr1, chr2, resolution){
    Reference_object <- GenomicMatrix$new()
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, resolution = resolution)
    ChromosomeList <- ChromInfo[,"chr"]
    if(!is.character(chr1) | !is.character(chr2)){
        stop("Provided Chromosomes does not appear to be of class character")
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop(chr1,chr2," matrix is yet to be loaded.")
    }
    Group_path <- Create_Path(c(Reference_object$hdf.matrices.root, 
        chr1, chr2))
    Brick_filepath <- BrickContainer_get_path_to_file(Brick = Brick, 
        chr1 = chr1, chr2 = chr2, resolution = resolution)
    dataset_handle <- ._Brick_Get_Something_(Group.path = Group_path, 
        Brick = Brick_filepath, Name = Reference_object$hdf.matrix.name, 
        return.what = "dataset_handle")
    entire_matrix <- dataset_handle[]
    CloseH5Con(Handle = dataset_handle, type = "dataset")
    return(entire_matrix)
}

#' Get the matrix metadata columns in the Brick store.
#'
#' `Brick_get_matrix_mcols` will get the specified matrix metadata column for
#' a chr1 vs chr2 Hi-C data matrix. Here, chr1 represents the rows and chr2
#' represents the columns of the matrix. For cis Hi-C matrices, where 
#' chr1 == chr2, chr2_bin_coverage and chr2_col_sums equals chr1_bin_coverage 
#' and chr1_row_sums respectively.
#' 
#' These metadata columns are: 
#' - chr1_bin_coverage: Percentage of rows containing non-zero values
#' - chr2_bin_coverage: Percentage of columns containing non-zero values
#' - chr1_row_sums: Total signal (if normalised) or number of reads 
#' (if counts) in each row.
#' - chr2_col_sums: Total signal (if normalised) or number of reads 
#' (if counts) in each column.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @inheritParams Brick_load_matrix
#'
#' @param what \strong{Required}
#' A character vector of length 1 specifying the matrix metric to retrieve
#' 
#' @return Returns a 1xN dimensional vector containing the specified matrix
#' metric
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_matrix_mcols_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_get_matrix_mcols(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", resolution = 100000, what = "chr1_bin_coverage")
#' 
Brick_get_matrix_mcols = function(Brick, chr1, chr2, resolution, 
    what = c("chr1_bin_coverage", "chr2_bin_coverage", 
        "chr1_row_sums", "chr2_col_sums")){
    what = match.arg(what)
    Reference.object <- GenomicMatrix$new()
    Meta.cols <- Reference.object$hdf.matrix.meta.cols()
    BrickContainer_class_check(Brick)
    if(any(is.null(c(chr1,chr2,what)))){
        stop("Brick, chr1, chr2, what cannot be NULL.\n")
    }
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1, chr2 = chr2, 
        resolution = resolution)){
        stop("Matrix for this chromsome pair does not exist.\n")  
    }
    if(!Brick_matrix_isdone(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution)){
        stop("Matrix for this chromsome pair is yet to be loaded.\n")  
    }
    if(length(what) >1){
        stop("What must be a character vector of length 1\n")     
    }
    if(!Brick_matrix_issparse(Brick = Brick, chr1 = chr1, chr2 = chr2,
        resolution = resolution) & what == Meta.cols["sparse"]){
        stop("This matrix is not a sparse matrix.",
            " So sparsity.index was not calculated\n")
    }
    if(what == Meta.cols["sparse"] & chr1 != chr2){
        stop("sparsity.index only applies to cis matrices (chr1 == chr2).\n")
    }
    Brick_filepath <- BrickContainer_get_path_to_file(Brick = Brick, 
        chr1 = chr1, chr2 = chr2, resolution = resolution)
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root, chr1,
        chr2))
    Vector <- ._Brick_Get_Something_(Group.path = Group.path, 
        Brick = Brick_filepath, Name = what, return.what = "data")
    return(Vector)
}

#' List the matrix metadata columns in the Brick store.
#'
#' `Brick_get_matrix_mcols` will list the names of all matrix metadata 
#' columns.
#' 
#' @return Returns a vector containing the names of all matrix metadata columns
#'
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "list_matrix_mcols_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Brick_list_matrix_mcols()
Brick_list_matrix_mcols = function(){
    Reference.object <- GenomicMatrix$new()
    Meta.cols <- Reference.object$hdf.matrix.meta.cols()
    return(Meta.cols)
}


#' Export an entire resolution from a given BrickContainer as a 
#' upper triangle sparse matrix
#'
#' `Brick_export_to_sparse` will accept as input an object of class 
#' BrickContainer, a string of length 1 as resolution and a path specifying
#' the output file to write. It writes the content of the all loaded Brick
#' objects as a upper triangle sparse matrix (col > row) containing 
#' non-zero values.
#'
#' @inheritParams Brick_get_chrominfo
#'
#' @param out_file Path to the output file to write.
#' 
#' @param sep column delimiter in output file. Default single space.
#' 
#' @param remove_file Default FALSE. If a file by the same name is present
#' that file will be removed.
#'
#' @return Returns a data.frame corresponding to the head of the output file
#' 
#' @examples
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "write_file")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_export_to_sparse(Brick = My_BrickContainer, 
#' out_file = file.path(out_dir, "example_out.txt"), 
#' resolution = 100000)
#' 
Brick_export_to_sparse <- function(Brick, out_file, remove_file = FALSE, 
    resolution, sep = " "){
    Reference.object <- GenomicMatrix$new()
    Chrominfo_df <- Brick_get_chrominfo(Brick, resolution = resolution)

    end_positions <- cumsum(Chrominfo_df[,"nrow"])
    start_positions <- c(1, end_positions[-length(end_positions)] + 1)
    names(start_positions) <- Chrominfo_df$chr
    chromosome_coord_pairs <- .make_coords_list(Brick, resolution)
    if(file.exists(out_file) & !remove_file){
        stop(out_file, " already exists at path")
    }
    open_connection <- file(out_file, open = "w")
    temp_df <- data.frame(chr1 = "chr1", 
            chr2 = "chr2", 
            row_coord = "row_coord",
            col_coord = "col_coord",
            value = "value")
    write.table(x = temp_df, file = open_connection, sep = sep, 
        row.names = FALSE, col.names = FALSE, quote = FALSE)
    lapply(seq_len(nrow(chromosome_coord_pairs)), function(x){
        a_row <- chromosome_coord_pairs[x,]
        a_vector <- .fetch_upper_tri_value_by_row(
            Brick_filepath = a_row$filepath,
            chr1 = a_row$chr1, 
            chr2 = a_row$chr2, 
            row = a_row$chr1_start,
            col = a_row$chr2_start, 
            chr2_length = a_row$chr2_length)
        a_vector <- a_vector[1,]
        non_zero_filter <- a_vector != 0
        row_seq <- rep(a_row$chr1_start, times = length(a_vector))
        row_seq <- row_seq[non_zero_filter] + (start_positions[a_row$chr1] - 1)
        col_seq <- seq(from = a_row$chr2_start, by = 1, 
            length.out = a_row$chr2_length)
        col_seq <- col_seq[non_zero_filter] + (start_positions[a_row$chr2] - 1)
        temp_df <- data.frame(chr1 = a_row$chr1, 
            chr2 = a_row$chr2, 
            row_coord = row_seq,
            col_coord = col_seq,
            value = a_vector[non_zero_filter],
            stringsAsFactors = FALSE)
        write.table(x = temp_df, file = open_connection, sep = sep, 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
    })
    close(open_connection)
    return(read.table(file = out_file, header = TRUE, sep = sep, 
            stringsAsFactors = FALSE, nrows = 100))
}


#' Identify compartments in the Hi-C data
#'
#' `Brick_call_compartments` identifies compartments in Hi-C data. Reference
#' Lieberman-Aiden et al. 2009.
#' 
#' @inheritParams Brick_get_values_by_distance
#' 
#' @return A dataframe containing the chromosome genomic coordinates and the 
#' first three principal components.
#' 
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_vector_val_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Compartments_df <- Brick_call_compartments(Brick = My_BrickContainer, 
#' chr = "chr2L", resolution = 100000)
#' 
Brick_call_compartments <- function(Brick, chr, resolution){
    configuration_length_check(chr, "chr", 1)
    configuration_length_check(resolution, "resolution", 1)
    BrickContainer_class_check(Brick)
    a_matrix <- .remove_nas(Brick_get_entire_matrix(Brick = Brick, 
        chr1 = chr, chr2 = chr, resolution = resolution))
    normalised_matrix <- .normalize_by_distance_values(a_matrix)
    correlation_matrix <- cor(normalised_matrix)
    correlation_matrix <- .remove_nas(correlation_matrix)
    bintable_df <- as.data.frame(Brick_get_bintable(Brick = Brick, 
        chr = chr, resolution = resolution))
    bintable_df <- bintable_df[,c("seqnames", "start", "end")]
    pca_list <- prcomp(correlation_matrix)
    bintable_df$pc1 <- pca_list[["x"]][,"PC1"]
    bintable_df$pc2 <- pca_list[["x"]][,"PC2"]
    bintable_df$pc3 <- pca_list[["x"]][,"PC3"]
    return(bintable_df)
}

#' Identify compartments in the Hi-C data
#'
#' `Brick_load_data_from_table` loads data from a table like file, such as
#' sparse matrices.
#' 
#' @inheritParams Brick_load_matrix
#' 
#' @param table_file Path to the file that will be loaded
#' 
#' @param col_index \strong{Optional}. Default "c(1,2,3)".
#' A character vector of length 3 containing the indexes of the required
#' columns in the table file. the first index, corresponds to bin1, the 
#' second to bin2 and the third to the signal value.
#' 
#' @return A dataframe containing the chromosome genomic coordinates and the 
#' first three principal components.
#' 
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "get_vector_val_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
#' chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_export_to_sparse(Brick = My_BrickContainer, 
#' out_file = file.path(out_dir, "example_out.txt"), 
#' remove_file = TRUE, sep = " ", 
#' resolution = 100000)
#' 
#' Brick_load_data_from_sparse(Brick = My_BrickContainer, 
#' table_file = file.path(out_dir, "example_out.txt"), 
#' delim = " ", resolution = 100000)
#' 
Brick_load_data_from_sparse <- function(Brick, table_file, delim = " ", 
    resolution = NULL, batch_size = 1000000, matrix_chunk = 2000, 
    remove_prior = FALSE) {
    Reference.object <- GenomicMatrix$new()
    col_index <- c(1, 2, 3)
    BrickContainer_class_check(Brick)
    Resolutions <- BrickContainer_list_resolutions(Brick)
    resolution <- .format_resolution(resolution)

    if(!(resolution %in% Resolutions)) {
        stop("resolution does not exist in the BrickContainer")
    }

    RetVar <- .process_tsv(Brick = Brick, table_file = table_file, 
        delim = delim, resolution = resolution, matrix_chunk = matrix_chunk, 
        batch_size = batch_size, remove_prior = remove_prior, 
        col_index = col_index, is_sparse = FALSE)
    return(RetVar)
}