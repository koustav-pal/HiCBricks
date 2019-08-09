BrickContainer <- setClass("BrickContainer", slots = list(name = "character",
    resolutions = "character", container_path = "character",
    chromosomes = "character", chromosome_length = "integer",
    file_list = "data.frame", headers = "list", matrix_info = "list",
    metadata_list = "list"
    ))

setMethod("show",
"BrickContainer",
function(object) {
message("Experiment name: ",object@name)
message("Project directory: ", getRelativePath(pathname = 
    object@headers$project_directory, 
    relativeTo = getwd(), 
    caseSensitive = TRUE))
message("Configuration file: ", getRelativePath(pathname = 
    object@container_path, 
    relativeTo = getwd(), 
    caseSensitive = TRUE))
message("Resolutions: ", paste(object@resolutions, collapse = ", "))
message("Chromosomes: ", paste(object@chromosomes, collapse = ", "))
message("Lengths: ", paste(object@chromosome_length, collapse = ", "))
num_files <- nrow(object@file_list)
# type_names <- names(num_type)
message("containing ", num_files, " matrices across ", 
    length(object@resolutions), " resolutions and ", 
    length(object@chromosomes)," chromosomes")
show(object@file_list)
})

return_experiment_name <- function(x){
        return(x@name)
}

return_output_directory <- function(x){
        return(x@headers$project_directory)
}

return_resolutions <- function(object) {
        return(.format_resolution(object@resolutions))
}

return_chromosomes <- function(object) {
    return(object@chromosomes)
}

return_chromosome_lengths <- function(object) {
        return(object@chromosome_length)
}

.return_file_list <- function(object) {
        return(object@file_list)
}


setMethod("return_configuration_header",
    "BrickContainer",
    function(config_file) {
        return(config_file@headers)
})

setMethod("return_configuration_matrix_info",
    "BrickContainer",
    function(config_file) {
        return(config_file@matrix_info)
})


create_configuration_object <- function(object){
        return(list(headers = object@headers, 
            matrix_info = object@matrix_info))
}

.write_configuration_file <- function(object, config_filepath){
        config_object = create_configuration_object(object)
        config_json = prettify(toJSON(config_object), indent = 4)
        write_lines(x = config_json, path = config_filepath)
}

BrickContainer_list_rangekeys <- function(object, resolution = NA, 
    all_resolutions = NA){
    File_list <- BrickContainer_list_files(object, resolution = resolution)
    Reference.object <- GenomicMatrix$new()
    File_list_colnames <- Reference.object$Configurator_JSON_matrix_names
    Brick_path <- File_list$filepaths
    Dataset <- ._Brick_Get_Something_(
        Group.path = Reference.object$hdf.metadata.root,
        Brick = Brick_path, Name = Reference.object$metadata.chrom.dataset,
        return.what = "data")
    return(Dataset)
}

BrickContainer_list_matrices <- function(object, resolution = NA, 
    all_resolutions = FALSE){
    BrickContainer_resolution_check(resolution, all_resolutions)
    File_list <- BrickContainer_list_files(Brick = object, 
        resolution = resolution)
    Reference.object <- GenomicMatrix$new()
    Colnames <- Reference.object$matrices.chrom.attributes
    Brick_df_list <- lapply(seq_along(File_list$filepaths), function(x){
        Brick <- File_list$filepaths[x]
        temp.df <- Brick_list_matrices(Brick = Brick, 
            chr1 = File_list$chrom1[x], 
            chr2 = File_list$chrom2[x])
        temp.df$resolution = File_list$resolution[x]
        return(temp.df)
    })
    Brick_df <- do.call(rbind, Brick_df_list)
    return(Brick_df)
}