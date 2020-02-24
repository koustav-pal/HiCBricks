#' Create a BrickContainer object from a JSON file
#'
#' `load_BrickContainer` creates a BrickContainer object from a JSON file
#'
#' @param config_file Default NULL
#' A character string of length 1 specifying the path to the path to the 
#' configuration json created using `Create_many_bricks`
#' 
#' @param project_dir Default NULL
#' A character string of length 1 specifying the path to the path to the 
#' configuration json created using `Create_many_bricks`
#' 
#' @return An object of class BrickContainer
#' 
#' @examples
#' 
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' out_dir <- file.path(tempdir(), "BrickContainer_load_test")
#' dir.create(out_dir)
#' Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' My_BrickContainer <- load_BrickContainer(project_dir = out_dir)
#' 
load_BrickContainer <- function(config_file = NULL, project_dir = NULL){
    Reference.object <- GenomicMatrix$new()
    if(is.null(config_file) & is.null(project_dir)){
        stop(" One of config_file or project_dir must be provided.")
    }
    if(!is.null(project_dir)){
        config_file <- file.path(project_dir, 
            Reference.object$brick.config.name)
    }
    if(!file.exists(config_file)){
        stop(config_file," does not exist")
    }
    Config_header <- return_configuration_header(config_file = config_file)
    Configuration_matrix_list <- return_configuration_matrix_info(
        config_file = config_file)
    Container <- .prepare_BrickContainer(Config_header, 
        Configuration_matrix_list, 
        config_file)
    return(Container)
}

#' Get the list of HDF files present in the Brick container.
#'
#' `BrickContainer_list_files` fetches the list of HDF files associated to a 
#' particular BrickContainer
#'
#' @param Brick \strong{Required}.
#' A string specifying the path to the Brick store created 
#' with Create_many_Bricks.
#'
#' @inheritParams Brick_load_matrix
#' @inheritParams Create_many_Bricks
#' 
#' @param type 
#' A value from one of cis, trans specifying the type of files to list
#' cis will list intra-choromosomal file paths and trans will list 
#' inter-chromosomal file paths.
#' 
#' @return A 5 column tibble containing chromosome pairs, Hi-C resolution, 
#' the type of Hi-C matrix and the path to a particular Hi-C matrix file.
#'
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' out_dir <- file.path(tempdir(), "BrickContainer_list_file_test")
#' dir.create(out_dir)
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'    bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'    experiment_name = "Vignette Test", resolution = 100000, 
#'    remove_existing = TRUE)
#' 
#' BrickContainer_list_files(Brick = My_BrickContainer, chr1 = "chr2L",
#' chr2 = NA)
#' 
BrickContainer_list_files <- function(Brick = NULL,
    chr1 = NA, chr2 = NA, type = NA, resolution = NA){
    BrickContainer_class_check(Brick)
    Files_tib <- .return_file_list(Brick)
    ChromosomeList <- return_chromosomes(Brick)
    ResolutionList <- return_resolutions(Brick)
    out_dir <- return_output_directory(Brick)
    if(!is.na(resolution)){
        resolution <- .format_resolution(resolution)
        if(!any(resolution %in% ResolutionList)){
            stop("all resolutions were not found in resolution list")
        }
        Files_tib <- Files_tib[Files_tib$resolution %in% resolution,]
    }
    if(all(!is.na(chr1))){
        if(!any(chr1 %in% ChromosomeList)){
            stop("chr1 was not found in chromosome list")
        }
        Files_tib <- Files_tib[Files_tib$chrom1 %in% chr1,]
    }
    if(all(!is.na(chr2))){
        if(!any(chr2 %in% ChromosomeList)){
            stop("chr2 was not found in chromosome list")
        }
        Files_tib <- Files_tib[Files_tib$chrom2 %in% chr2,]
    }
    if(!is.na(type)){
        if(!(type %in% c("cis", "trans"))){
            stop("type must be one of cis or trans")
        }
        Files_tib <- Files_tib[Files_tib$mat_type %in% type,]
    }
    Files_tib$filepaths <- file.path(out_dir, Files_tib$filename)
    # message("Found ", nrow(Files_tib), " files...")
    return(Files_tib)
}

#' Get the path to HDF files present in the Brick container.
#'
#' `BrickContainer_get_path_to_file` fetches the list of HDF file paths 
#' associated to a particular BrickContainer
#'
#' @inheritParams BrickContainer_list_files
#' 
#' @return A vector containing filepaths
#'
#' @examples
#' 
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_list_filepath_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_get_path_to_file(Brick = My_BrickContainer, chr1 = "chr2L",
#' chr2 = "chr2L", resolution = 100000)
#' 
BrickContainer_get_path_to_file <- function(Brick = NULL,
    chr1 = NA, chr2 = NA, type = NA, resolution = NA){
    File_list <- BrickContainer_list_files(Brick = Brick,
        chr1 = chr1, 
        chr2 = chr2, 
        type = type, 
        resolution = resolution)
    return(File_list$filepaths)
}

#' Change the location of HDF files in the BrickContainer object
#'
#' `BrickContainer_change_experiment_name` changes the location of 
#' name of the experiment
#'
#' @param Brick \strong{Required}.
#' A string specifying the path to the BrickContainer created using 
#' `Create_many_Bricks` or `Load_BrickContainer`
#'
#' @param experiment_name \strong{Required}. Default NULL
#' A string specifying the new experiment name
#'
#' @return An object of class BrickContainer where the experiment_name 
#' has been changed
#' 
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_expt_name_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_change_experiment_name(Brick = My_BrickContainer, 
#' experiment_name = "I change my mind")
#' 
BrickContainer_change_experiment_name <- function(Brick = NULL, 
    experiment_name = NULL) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    configuration_null_check(value = experiment_name, 
        name = "experiment_name")
    configuration_length_check(value = experiment_name,
        name = "experiment_name", length_threshold = 1)
    Brick@headers$experiment_name <- experiment_name
    Config_filepath <- .make_configuration_path(
        return_output_directory(Brick))
    .write_configuration_file(Brick, Config_filepath)
    Brick <- load_BrickContainer(Config_filepath)
    return(Brick)
}

#' Change the location of HDF files in the BrickContainer object
#'
#' `BrickContainer_change_output_directory` changes the location of 
#' associated HDF files
#'
#' @param Brick \strong{Required}.
#' A string specifying the path to the Brick store created with CreateBrick.
#'
#' @param output_directory \strong{Required}. Default NULL
#' A string specifying new location of the output_directory. Please note, 
#' that the location of the HDF files will not be changed.
#'
#' @return An object of class BrickContainer where the output_directory 
#' has been changed
#' 
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_out_dir_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_change_output_directory(Brick = My_BrickContainer, 
#' output_directory = tempdir())
BrickContainer_change_output_directory <- function(Brick = NULL, 
    output_directory = NULL) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    configuration_null_check(value = output_directory, 
        name = "output_directory")
    configuration_length_check(value = output_directory, 
        name = "output_directory", length_threshold = 1)
    Brick@headers$project_directory <- normalizePath(output_directory)
    Config_filepath <- .make_configuration_path(output_directory)
    .write_configuration_file(Brick, Config_filepath)
    Brick <- load_BrickContainer(Config_filepath)
    return(Brick)
}

#' Remove links to all Hi-C matrices for a given resolution.
#'
#' `BrickContainer_unlink_resolution` removes links to all files associated 
#' to a given resolution
#'
#' @param Brick \strong{Required}.
#' A string specifying the path to the Brick store created with CreateBrick.
#'
#' @param resolution \strong{Required}
#' A string specifying the resolution to remove. This string must match the 
#' resolution values listed by `BrickContainer_list_resolutions`
#'
#' @return An object of class BrickContainer where the resolution 
#' and links to its associated files have been removed
#' 
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_unlink_res_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 40000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_unlink_resolution(Brick = My_BrickContainer, 
#' resolution = 40000)
#' 
BrickContainer_unlink_resolution <- function(Brick = NULL, resolution = NULL) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    Resolutions <- BrickContainer_list_resolutions(Brick)
    configuration_length_check(resolution, "resolution", 1)
    configuration_null_check(value = resolution, name = "resolution")
    resolution <- .format_resolution(resolution)
    if(!(resolution %in% Resolutions)){
        stop("Provided resolution is not present in the BrickContainer")
    }
    Header <- return_configuration_header(Brick)
    Files_tib <- BrickContainer_list_files(Brick = Brick)
    original_length <- nrow(Files_tib)
    Keep_which <- which(!(Files_tib$resolution %in% resolution))
    Matrix_info <- return_configuration_matrix_info(Brick)
    Matrix_info <- Matrix_info[Keep_which]
    Files_tib <- .create_file_list(Matrix_info)
    Header$resolutions <- Resolutions[!(Resolutions %in% resolution)]
    Brick@file_list <- Files_tib
    Brick@matrix_info <- Matrix_info
    Brick@headers <- Header
    Config_filepath <- .make_configuration_path(
        return_output_directory(Brick))
    .write_configuration_file(Brick, Config_filepath)
    Return_Brick <- load_BrickContainer(config_file = Config_filepath)
    message("Removed links to ", original_length - nrow(Files_tib),
        " files under resolution ", resolution)
    message("Original files are still present at location! ", 
        "Please, remove them manually.")
    return(Return_Brick)
}

#' Return the descriptive name of the BrickContainer
#'
#' `BrickContainer_list_experiment_name` returns the descriptive name of a
#' BrickContainer
#'
#' @inheritParams BrickContainer_list_files
#'
#' @return A character string specifying the descriptive name of the 
#' BrickContainer
#' 
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_list_expt_name_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_list_experiment_name(My_BrickContainer)
BrickContainer_list_experiment_name <- function(Brick = NULL) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    return(return_experiment_name(Brick))
}

#' Return the output directory of the BrickContainer
#'
#' `BrickContainer_list_output_directory` returns the location of the 
#' associated HDF files
#'
#' @inheritParams BrickContainer_list_files
#'
#' @return A character string specifying the descriptive name of the 
#' BrickContainer
#' 
#' @examples
#' 
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_list_out_dir_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_list_output_directory(My_BrickContainer)
BrickContainer_list_output_directory <- function(Brick = NULL) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    return(return_output_directory(Brick))
}


#' Return the descriptive name of the BrickContainer
#'
#' `BrickContainer_list_resolutions` returns the resolutions available in 
#' the BrickContainer
#'
#' @inheritParams BrickContainer_list_files
#'
#' @return A character string specifying the descriptive name of the 
#' BrickContainer
#' 
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_list_resolution_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_list_resolutions(My_BrickContainer)
BrickContainer_list_resolutions <- function(Brick = NULL) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    return(return_resolutions(Brick))
}

#' Return the descriptive name of the BrickContainer
#'
#' `BrickContainer_list_chromosomes` returns the chromosomes available in 
#' the BrickContainer
#'
#' @inheritParams BrickContainer_list_files
#'
#' @param lengths Default FALSE
#' If TRUE, will also return the chromosomal lengths
#' 
#' @return If lengths is FALSE, only the chromosome names are returned. If 
#' lengths is TRUE, a data.frame containing the chromosome names and their 
#' lengths is provided.
#' 
#' @examples
#' Bintable.path <- system.file("extdata",
#' "Bintable_100kb.bins", package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "BrickContainer_list_chromosome_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
#'   experiment_name = "Vignette Test", resolution = 100000, 
#'   remove_existing = TRUE)
#' 
#' BrickContainer_list_chromosomes(My_BrickContainer)
BrickContainer_list_chromosomes <- function(Brick = NULL, lengths = FALSE) {
    configuration_null_check(Brick, "Brick")
    BrickContainer_class_check(Brick)
    ChromosomeList <- return_chromosomes(Brick)
    if(!lengths){
        return(ChromosomeList)
    }
    ChromosomeLengths <- return_chromosome_lengths(Brick)
    return(data.frame(chrom = ChromosomeList, lengths = ChromosomeLengths,
        stringsAsFactors = FALSE))
}

# #' Return cis interactions separated by a certain distance
# #'
# #' `BrickContainer_get_values_by_distance` returns cis interactions separated
# #' by a certain distance. This is a parallel version of 
# #' Brick_get_values_by_distance. It also lacks the constrain region feature
# #' that is found in the non-parallel version. 
# #'
# #' @inheritparam Brick_get_values_by_distance
# #' @inheritparam BrickContainer_list_files
# #'
# #' @param chrs 
# #' List of chromosomes for which the distance vectors will be fetched
# #' 
# #' @param all_resolutions Default FALSE
# #' If TRUE, distance values will be fetched for all resolutions listed in the
# #' BrickContainer
# #' 
# #' @param num_cpus Default 1
# #' Maximum number of parallel jobs that can be run
# #' 
# #' @return A tibble containing colnames chr, resolution and fun_value. The
# #' three columns specify the chromosome name, resolution of the matrix and
# #' the distance value being fetched.
# #' 
# #' @examples
# #' Brick.file = system.file("extdata", "test.json", package = "HiCBricks")
# #' My_brick_container <- load_BrickContainer(Brick.file)
# #' Brick_parallel_get_values_by_distance(My_brick_container)
# Brick_parallel_get_values_by_distance <- function(Brick, chrs, distance,
#     batch.size=500, FUN=NULL, resolution = NA, all_resolutions = FALSE, 
#     num_cpus = 1){
#     Reference.object <- GenomicMatrix$new()
#     if(!BrickContainer_class_check(Brick)) {
#         stop("Brick must be of class BrickContainer")
#     }
#     BrickContainer_resolution_check(resolution, all_resolutions)
#     Files_df <- BrickContainer_list_files(Brick = Brick, type = "cis",
#         chr1 = chrs, resolution = resolution)
#     Files_max_distance <- do.call(c, bplapply(seq_len(nrow(Files_df)), 
#             function(file_path_x){
#             Brick_matrix_maxdist(Brick = Files_df[file_path_x,"filepaths"],
#                 chr1 = Files_df$chr1[file_path_x],
#                 chr2 = Files_df$chr1[file_path_x]) > distance
#         })
#     )
#     if(any(Files_max_distance)){
#         message("Skipping out of bounds chromosomes, ", 
#             paste(Files_df$chr1[Files_max_distance], collapse = ", "))
#     }
#     Files_df <- Files_df[!Files_max_distance,]
#     if(nrow(Files_df) == 0 & any(Files_max_distance)){
#         stop("No files available for processing. ",
#             "Chromosomes were out of bounds")
#     }
#     Result_tib_list <- bplapply(seq_along(Files_df_split), function(x){
#         Brick_path <- cur_files_df$filepaths[x]
#         chr <- cur_files_df$chr1[x]
#         cur_resolution <- cur_files_df$resolution[x]
#         Value <- Brick_get_values_by_distance(Brick = Brick_path, 
#             chr = chr,
#             distance = distance,
#             batch.size = batch.size,
#             FUN = FUN)
#         tibble(chr = chr,
#             resolution = cur_resolution,
#             fun_value = Value)
#     }, BPPARAM = .get_instance_biocparallel(workers = num_cpus))
#     Result_tib <- do.call(rbind, Result_tib_list)
#     return(Result_tib)
# }


# #' Return a list of sub-matrices for a given set of genomic loci. 
# #'
# #' Provided a list of genomic loci, `BrickContainer_get_matrix_within_coords` 
# #' returns a list of sub-matrices by calling Brick_get_matrix_within_coords
# #' in parallel.
# #'
# #' @inheritparam Brick_get_values_by_distance
# #' @inheritparam BrickContainer_list_files
# #'  
# #' @return A named ordered list of submatrices or submatrices with FUN 
# #' applied with names x.coords:y.coords. The order provided is that of the 
# #' x.coords and y.coords vector provided.
# #' 
# #' @examples
# #' Brick.file = system.file("extdata", "test.json", package = "HiCBricks")
# #' My_brick_container <- load_BrickContainer(Brick.file)
# #' Brick_parallel_get_matrix_within_coords(My_brick_container)
# Brick_parallel_get_matrix_within_coords <- function(Brick, x.coords,
#     y.coords, force = FALSE, FUN=NULL, resolution = NA, num_cpus = 1){
#     Reference.object <- GenomicMatrix$new()
#     if(length(x.coords)!=length(y.coords)){
#         stop("x.coords and y.coords length do not match")
#     }
#     ChromosomeList <- BrickContainer_list_chromosomes(Brick)
#     configuration_na_check(resolution, "resolution")
#     configuration_length_check(resolution, "resolution", 1)
#     xcoords.split <- Split_genomic_coordinates(Coordinate=x.coords)
#     ycoords.split <- Split_genomic_coordinates(Coordinate=y.coords)
#     chr1 <- vapply(xcoords.split, function(x){ x[1] },"chr1")
#     chr2 <- vapply(ycoords.split, function(x){ x[1] },"chr1")
#     chr1_positions <- match(chr1, ChromosomeList)
#     chr2_positions <- match(chr2, ChromosomeList)
#     if(any(is.na(chr1_positions) | is.na(chr2_positions))){
#         Message_1 <- NULL
#         Message_2 <- NULL
#         if(any(is.na(chr1_positions))){
#             Message_1 <- paste("Problem in x.coords position",
#                 paste(which(is.na(chr1_positions)), collapse = ", "))
#         }
#         if(any(is.na(chr2_positions))){
#             Message_2 <- paste("Problem in y.coords position",
#                 paste(which(is.na(chr2_positions)), collapse = ", "))
#         }
#         stop("Not all chromosomes were found in the chromosome list.",
#             Message_1, Message_2)
#     }
#     if(any(chr1_positions > chr2_positions)){
#         Which_chr1 <- which(chr1_positions > chr2_positions)
#         stop("Not all positions of chr2 were greater than chr1 in",
#             "chromosome list. Problem at coordinate locations",
#             paste(Which_chr1, collapse = ", "))
#     }
#     x_coords_split <- split(x.coords, paste(chr1, chr2, sep = "_"))
#     y_coords_split <- split(y.coords, paste(chr1, chr2, sep = "_"))
#     chr1_split <- split(chr1, paste(chr1, chr2, sep = "_"))
#     chr2_split <- split(chr2, paste(chr1, chr2, sep = "_"))
#     Matrices_list_list <- bplapply(seq_along(x_coords_split), function(x){
#         x_coords <- x.coords[x]
#         y_coords <- y.coords[x]
#         cur_chr1 <- chr1_split[x]
#         cur_chr2 <- chr1_split[x]
#         current_Brick <- BrickContainer_get_path_to_file(Brick, 
#             chr1 = unique(chr1), 
#             chr2 = unique(chr2),
#             resolution = resolution)
#         matrix_subset_list <- lapply(seq_along(x_coords),function(y){
#             Brick_get_matrix_within_coords(current_Brick, 
#                 x.coords = x_coords[y],
#                 y.coords = y_coords[y],
#                 force = force,
#                 FUN = FUN)
#         })
#         names(matrix_subset_list) <- paste(x_coords, y_coords,
#             sep = Reference.object$Ranges.separator)
#     }, BPPARAM = .get_instance_biocparallel(workers = num_cpus))
#     Matrices_list <- do.call(c, Matrices_list_list)
#     Matrices_list <- Matrices_list[paste(x.coords, 
#             y.coords, 
#             sep = Reference.object$Ranges.separator)]
#     return(Matrices_list)
# }