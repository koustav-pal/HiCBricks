configuration_length_check <- function(value, name, length_threshold){
    if(length(value) > length_threshold){
        stop(name, " cannot have more than ", length_threshold," argument.")
    }
}

configuration_null_check <- function(value, name){
    if(any(is.null(value))){
        stop(name," cannot have null values")
    }
}
configuration_na_check <- function(value, name){
    if(any(is.na(value))){
        stop(name," cannot have na values")
    }
}

BrickContainer_class_check <- function(x){
    if(!("BrickContainer" %in% class(x))){
        stop("An object of class BrickContainer is expected!",
            " BrickContainer objects can be created using Create_many_bricks",
            " and loaded to R sessions using load_BrickContainer")
    }
}

BrickContainer_resolution_check <- function(resolution, all_resolutions){
    if(is.na(resolution) & !all_resolutions){
        stop("When Brick is of class BrickContainer resolution must be",
            " provided or all_resolutions must be set to TRUE")
    }
}

configuration_class_check <- function(value, name){
    if(!BrickContainer_class_check(value)){
        stop(name," must be an object of class BrickContainer")
    }
}

read_configuration_file <- function(config_file){
    Configuration_list <- fromJSON(minify(read_lines(file = config_file)))
    return(Configuration_list)
}

setGeneric("return_configuration_resolution_info",
    function(config_file) {
    Configuration_list <- read_configuration_file(config_file)
    return(Configuration_list[["resolution_info"]])
})

setGeneric("return_configuration_matrix_info",
    function(config_file) {
    Configuration_list <- read_configuration_file(config_file)
    return(Configuration_list[["matrix_info"]])
})

setGeneric("return_configuration_header",
    function(config_file) {
    Configuration_list <- read_configuration_file(config_file)
    return(Configuration_list[["headers"]])
})


.create_configuration_header <- function(file_prefix, output_directory, 
    experiment_name, resolution, chromosomes, chromosome_lengths){
    Reference.object <- GenomicMatrix$new()
    header_colnames <- Reference.object$Configurator_JSON_headers_names
    headers = list(
        file_prefix,
        output_directory,
        experiment_name,
        resolution,
        chromosomes, 
        as.integer(chromosome_lengths),
        NA)
    names(headers) <- header_colnames
    return(headers)
}

.create_configuration_resolution_info <- function(rangekeys, 
    ranges_filenames, chrominfo_filename){
    Reference.object <- GenomicMatrix$new()
    resolution_colnames <- 
        Reference.object$Configurator_JSON_resolution_names
    resolution_list <- list(
        rangekeys,
        ranges_filenames,
        chrominfo_filename)
    names(resolution_list) <- resolution_colnames
    return(resolution_list)
}

.create_configuration_matrix_info <- function(resolution, chrom1, chrom2,
    chrom1_binned_length, chrom2_binned_length, chrom1_max_size, 
    chrom2_max_size, attributes_list, type, filename){
    Reference.object <- GenomicMatrix$new()
    matrix_colnames <- Reference.object$Configurator_JSON_matrix_names
    chrom1_chrom2_list <- list(
        chrom1,
        chrom2,
        resolution,
        c(chrom1_binned_length, chrom2_binned_length),
        c(chrom1_max_size, chrom2_max_size),
        attributes_list,
        type,
        filename
    )
    names(chrom1_chrom2_list) <- matrix_colnames
    return(chrom1_chrom2_list)
}

.create_file_list <- function(matrix_info){
    Reference.object <- GenomicMatrix$new()
    matrix_colnames <- Reference.object$Configurator_JSON_matrix_names
    Matrix_df_list <- lapply(matrix_info, function(a_row){
        a_tibble <- data.frame(a_row[matrix_colnames[1]], 
        a_row[matrix_colnames[2]],
        a_row[matrix_colnames[3]],
        a_row[matrix_colnames[7]],
        a_row[matrix_colnames[8]], 
        stringsAsFactors = FALSE)
    })
    Matrix_df <- do.call(rbind, Matrix_df_list)
    colnames(Matrix_df) <- matrix_colnames[c(1,2,3,7,8)]
    Matrix_tibble <- as_tibble(Matrix_df)
    return(Matrix_tibble)
}

.prepare_BrickContainer <- function(header, matrix_info, 
    resolution_info, config_path){
    File_tib <- .create_file_list(matrix_info)
    BrickContainer <- new("BrickContainer",
    name = header$experiment_name,
    resolutions = header$resolutions,
    container_path = config_path,
    chromosomes = header$chromosomes,
    chromosome_length = header$lengths,
    file_list = File_tib,
    headers = header,
    resolution_info = resolution_info,
    matrix_info = matrix_info)
    return(BrickContainer)
}

.make_configuration_path <- function(output_directory){
    Reference.object <- GenomicMatrix$new()
    Config_filename <- Reference.object$brick.config.name
    Config_filepath <- file.path(normalizePath(output_directory),
        Config_filename)
    return(Config_filepath)
}


.format_resolution <- function(x){
    return(trimws(format(x, scientific = FALSE)))
}
