ReturnH5FileConnection = function(File = NULL){
    require(rhdf5)
        HDF.File <- File
        private$HDF.Connection = H5Fopen(name=HDF.File)
        return(private$HDF.Connection)
}
ReturnH5Handler = function(Path=NULL,File = NULL){
        Temp.connection = private$ReturnH5FileConnection(File=File)
        # Path.to.file <- ""
        # for (grp in Group) {
        #     Path.to.file <- file(Path.to.file,grp)
        # }
        Handle=Temp.connection&Path
        return(Handle)
}
Create_Path <- function(groups){
   Path.to.file <- ""
    for (x in groups) {
        Path.to.file <- file.path(Path.to.file,x)
    }
    return(Path.to.file)
}
CreateGroups <- function(Group.path = NULL, File = NULL){
    # Attribute must be a list of max level 2 containing values corresponding to each group in named elements
    h5createGroup(file = File, group = Group.path)

}
CreateDataset <- function(Path = NULL, File = NULL, name = NULL, dims = NULL, maxdims = NULL,
 data_type = "H5T_NATIVE_DOUBLE", chunk_size = NULL, compression_level = 6, fill_with = 0, warnings = TRUE){
    Full.path <- Create_Path(Path)
    if(is.null(chunk_size)){
        chunk_size <- ceiling(Dims/100)
    }
    H5Handler <- ReturnH5Handler(Path=Full.path,File = File)
    h5createDataset(file=H5Handler, dataset=name, dims=dims,
        maxdims = dims, H5type=data_type, chunk = chunk_size,
        level = compression_level, fillValue = fill_with, showWarnings = warnings)
}
CreateAttributes <- function(Path = NULL, File = NULL, Attributes = NULL, data_type = NULL){
    Lego.handler <- ReturnH5Handler(Path = Path,File = File)
    dtype_id <- "H5T_STRING"
    if(data_type != "character"){
        stop("For now we are not considering non-character attributes.\n")
    }
    for (An.attribute in Attributes) {
        H5Acreate(h5obj = Lego.handler, name = An.attribute, dtype_id = dtype_id, h5space = H5Screate_simple(dims = c(1),maxdims = c(1)))
    }
    H5Fclose(Lego.handler)
}
