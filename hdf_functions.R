ReturnH5FileConnection = function(File = NULL){
    require(rhdf5)
        HDF.File <- File
        HDF.Connection = H5Fopen(name=HDF.File)
        return(HDF.Connection)
}
ReturnH5Handler = function(Path=NULL,File = NULL){
        Temp.connection = ReturnH5FileConnection(File=File)
        Handle=Temp.connection&Path
        return(Handle)
}
Create_Path <- function(groups){
   Path.to.file <- ""
    for (x in groups) {
        Path.to.file <- file.path(Path.to.file,x)
    }
    # cat(Path.to.file,"\n")
    return(Path.to.file)
}
CreateGroups <- function(Group.path = NULL, File = NULL){
    Temp.connection = ReturnH5FileConnection(File=File)
    h5createGroup(file = Temp.connection, group = Group.path)
    H5Fclose(Temp.connection)
}
CreateDataset <- function(Path = NULL, File = NULL, name = NULL, dims = NULL, maxdims = NULL,
 data_type = "H5T_NATIVE_DOUBLE", chunk_size = NULL, compression_level = 6, fill_with = 0, warnings = TRUE){
    Full.path <- Create_Path(Path)
    if(is.null(chunk_size)){
        chunk_size <- ceiling(dims/100)
    }
    H5Handler <- ReturnH5Handler(Path=Full.path,File = File)
    h5createDataset(file=H5Handler, dataset=name, dims= dims,
        maxdims = dims, H5type=data_type, chunk = chunk_size,
        level = compression_level, fillValue = fill_with, showWarnings = warnings)
    H5Gclose(H5Handler)
}
CreateAttributes <- function(Path = NULL, File = NULL, Attributes = NULL, data_type = "character", on = "group"){
    Lego.handler <- ReturnH5Handler(Path = Path,File = File)
    dtype_id <- "H5T_NATIVE_CHAR"
    if(data_type != "character"){
        stop("For now we are not considering non-character attributes.\n")
    }
    for (An.attribute in Attributes) {
        H5Acreate(h5obj = Lego.handler, name = An.attribute, dtype_id = dtype_id, 
            h5space = H5Screate_simple(dims = c(1),maxdims = c(1)))
    }
    if(on == "group"){
        H5Gclose(Lego.handler)
    }else if(on == "dataset"){
        H5Dclose(Lego.handler)
    }
}
