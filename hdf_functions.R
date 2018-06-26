ReturnH5FileConnection = function(File = NULL){
    require(rhdf5)
        HDF.File <- File
        private$HDF.Connection = H5Fopen(name=HDF.File)
        return(private$HDF.Connection)
}
ReturnH5Handler = function(Path=NULL,File = NULL){
        Temp.connection = private$ReturnH5FileConnection()
        Path.to.file <- ""
        for (grp in Group) {
            Path.to.file <- file(Path.to.file,grp)
        }
        Handle=Temp.connection&Path.to.file
        return(Handle)
}
CreateGroups <- function(Groups = NULL, File = NULL){
    # Attribute must be a list of max level 2 containing values corresponding to each group in named elements
    GroupPath <- file.path(Groups)
    h5createGroup(file = File, group = GroupPath)

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
