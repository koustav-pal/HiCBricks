Brick_load_data_from_table <- function(Brick, mcool, chr1,
    chr2, binsize = NULL, cooler.batch.size = 1000000,
    matrix.chunk = 2000, dont.look.for.chr2 = FALSE, 
    remove.prior = FALSE, norm.factor = "Iterative-Correction"){
    Reference.object <- GenomicMatrix$new()
    if(is.null(chr1) | is.null(chr2)){
        stop("chr1, chr2 cannot be NULL.\n")
    }
    if(length(chr1) != length(chr2) | length(chr1) != 1){
        stop("chr1, chr2 are expected to be of length 1.\n")  
    }
    if(!Brick_matrix_exists(Brick = Brick, chr1 = chr1,
        chr2 = chr2)){
        stop("Provided chromosomes do not exist in the chrominfo table\n")
    }
    if(Brick_matrix_isdone(Brick = Brick, chr1 = chr1,
        chr2 = chr2) & !remove.prior){
        stop("A matrix was preloaded before. ",
            "Use remove.prior = TRUE to force value replacement\n")
    }
    resolutions <- Brick_list_mcool_resolutions(mcool = mcool)
    if(!is.null(resolutions) & is.null(binsize)){
        stop("binsize must be provided when ",
            "different resolutions are present in an mcool file.\n")
    }
    if(!is.null(resolutions) & !is.null(binsize)){
        if(!(as.character(as.integer(binsize)) %in% resolutions)){
            stop("binsize not found in mcool file. ",
                "Please check available binsizes ",
                "with Brick_list_mcool_resolutions.\n")            
        }
    }
    if(!is.null(norm.factor)){
        All.factors <- Brick_list_mcool_normalisations(names.only=TRUE)
        if(!Brick_mcool_normalisation_exists(mcool = mcool,
            norm_factor = norm.factor, resolution = as.character(
        as.integer(binsize)))){
            Factors.present <- All.factors[vapply(All.factors, 
                function(x){
                    Brick_mcool_normalisation_exists(mcool = mcool,
                    norm_factor = x, resolution = as.character(
                    as.integer(binsize)))
                }, TRUE)]
            stop(norm.factor," was not found in this mcool file.\n",
                "Normalisation factors present in the mcool file are: ",
                paste(Factors.present, collapse = ", "))
        }
        Norm.factors <- Brick_list_mcool_normalisations()
        Norm.factor <- Norm.factors[norm.factor]
        names(Norm.factor) <- NULL
    }else{
        Norm.factor <- NULL
    }
    RetVar <- ._Process_tsv(Brick = Brick, File = mcool,
        cooler.batch.size = cooler.batch.size, binsize = as.character(
        as.integer(binsize)), matrix.chunk = matrix.chunk, chr1 = chr1,
    chr2 = chr2, dont.look.for.chr2 = dont.look.for.chr2,
    norm.factor = Norm.factor, resolution = !is.null(resolutions))
    return(RetVar)
}
