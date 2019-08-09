.get_instance_biocparallel <- function(workers = 1){
    OS_type <- tolower(Sys.info()["sysname"])
    if(OS_type %in% c("linux", "osx", "darwin")){
        return(BiocParallel::register(MulticoreParam(workers=workers)))
    }else if(OS_type %in% c("windows")){
        return(BiocParallel::register(SerialParam(workers=workers)))
    }
}