# HDF structure

#' Create the entire HDF5 structure and load the bintable
#' 
#' `CreateLego` creates the complete HDF5 on-disk data structure
#' 
#' This function creates the complete HDF data structure, loads the binning
#' table associated to the Hi-C experiment and creates (for now) a 2D matrix
#' layout for all chromosome pairs. **Please note**, the binning table must
#' be a discontinuous one (first range end != secode range start), as ranges
#' overlaps using the "any" form will routinely identify adjacent ranges with
#' the same end and start to be in the overlap. Therefore, this criteria is 
#' enforced as default behaviour.
#' 
#' 
#' @param ChromNames \strong{Required} 
#' A character vector containing the chromosomes to be considered for the 
#' dataset. This string is used to verify the presence of all chromosomes in 
#' the provided bitable. 
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
#' @param Output.Filename \strong{Required}
#' A string specifying the location and name of the HDF file to create. If path
#' is not provided, it will be created in the Bioc File cache. Otherwise, it
#' will be created in the specified directory and tracked via Bioc File Cache.
#' 
#' @param bin.delim \strong{Optional}. Defaults to tabs.
#' A character vector of length 1 specifying the delimiter used in the file 
#' containing the binning table.
#' 
#' @param col.index \strong{Optional}. Default "c(1,2,3)".
#' A character vector of length 3 containing the indexes of the required 
#' columns in the binning table. the first index, corresponds to the chr 
#' column, the second to the start column and the third to the end column.
#' 
#' @param impose.discontinuity \strong{Optional}. Default TRUE.
#' If TRUE, this parameter ensures a check to make sure that required the end
#' and start coordinates of consecutive entries are not the same per 
#' chromosome.
#' 
#' @param ChunkSize \strong{Optional}.
#' A numeric vector of length 1. If provided, the HDF dataset will use this 
#' value as the chunk size, for all matrices. By default, the ChunkSize is 
#' set to matrix dimensions/100. 
#' 
#' @param exec \strong{Optional}. Default cat.
#' A string specifying the program or expression to use for reading the file.
#' For bz2 files, use bzcat and for gunzipped files use zcat.
#' 
#' @param remove.existing \strong{Optional}. Default FALSE.
#' If TRUE, will remove the HDF file with the same name and create a new one.
#' By default, it will not replace existing files.
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
#'            \item bin.coverage - \strong{dataset} - proportion of cells with
#' values greater than 0
#'            \item row.sums - \strong{dataset} - total sum of all values in a 
#' row
#'            \item sparsity - \strong{dataset} - proportion of non-zero cells 
#' near the diagonal
#'        }
#'    }
#'    \item Base.ranges - \strong{group}, Ranges tables for quick and easy 
#' access. Additional ranges tables are added here under separate group names.
#'    \itemize{
#'        \item Bintable - \strong{group} - The main binning table associated 
#' to a Lego.
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
#' specifying the chromosomes present in this particular Lego file.
#'        \item other metadata tables.
#'    }
#'}
#' @return This function will generate the target Lego file. Upon completion, 
#' the function will provide the path to the created/tracked HDF file.
#' 
#' @examples 
#' Bintable.path <- system.file("extdata", 
#' "Bintable_40kb.txt", package = "HiCLegos")
#' Chromosomes <- "chr19"
#' Path_to_cached_file <- CreateLego(ChromNames = Chromosomes, 
#' BinTable = Bintable.path, bin.delim = " ", 
#' Output.Filename = "test.hdf", exec = "cat", 
#' remove.existing = TRUE)
#' 
#' \dontrun{
#' Bintable.path <- system.file("extdata", 
#' "Bintable_40kb.txt", package = "HiCLegos")
#' Chromosomes <- c("chr19", "chr20", "chr22", "chr21")
#' Path_to_cached_file <- CreateLego(ChromNames = Chromosomes, 
#' BinTable = Bintable.path, impose.discontinuity=TRUE, 
#' col.index = c(1,2,3), Output.Filename = "test.hdf", 
#' exec = "cat", remove.existing = TRUE)
#' 
#' This will cause an error as the file located at Bintable.path,
#' contains coordinates for only chromosome 19. For this code to work, either
#' all other chromosomes need to be removed from the Chromosomes variable or
#' coordinate information for the other chromosomes need to be provided.
#' 
#' Similarly vice-versa is also true. If the Bintable contains data for other
#' chromosomes, but they were not listed in ChromNames, this will cause an 
#' error.
#' 
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
#' }
#' 
CreateLego <- function(ChromNames=NULL, BinTable=NULL, bin.delim="\t",
    col.index=c(1,2,3), impose.discontinuity=TRUE, ChunkSize=NULL, 
    Output.Filename=NULL, exec="cat", remove.existing=FALSE){
    H5close()
    Dir.path <- dirname(Output.Filename)
    Filename <- basename(Output.Filename)
    Output.Filename <- file.path(normalizePath(Dir.path),Filename)
    Reference.object <- GenomicMatrix$new()
    CreateCacheObject <- FALSE
    CreateNewCacheObject <- FALSE
    Cache.dir <- ._Get_cachedir()
    if(Dir.path == "."){
        warning("HDF file will be created and tracked by",
            "BioC cache directory.\n")
        if(Lego_is_tracked(Lego = Output.Filename)){
            Working.File <- Lego_list_tracked_legos()[Output.Filename]
            CreateNewCacheObject <- TRUE
        }else{
            Working.File <- bfcnew(x = Cache.dir, 
            rname = Output.Filename, ext = Reference.object$lego.extension,
            rtype = "local")
            CreateNewCacheObject <- FALSE
        }
    }else{
        Working.File <- Output.Filename
        CreateCacheObject <- TRUE
    }
    Root.folders <- Reference.object$GetRootFolders()
    if(is.null(ChromNames) | length(ChromNames) == 0){
        stop("Variable ChromNames cannot be empty")
    }
    HDF.File <- Working.File
    if(file.exists(HDF.File)){
        if(remove.existing){
            file.remove(HDF.File)
            Lego_untrack_lego(Lego = Output.Filename)
            if(CreateNewCacheObject){
                HDF.File <- bfcnew(x = Cache.dir, 
                rname = Output.Filename, 
                ext = Reference.object$lego.extension,
                rtype = "local")
            }
        }else{
            stop("Provided HDF file already exists. Please provide ",
                "remove.existing = TRUE to overwrite it\n")
        }
    }
    ChromosomeList<-ChromNames
    if(is.null(BinTable)){
        stop("Variable Bintable cannot be empty. Binning information must be ",
            "provided at startup\n")
    }
    # Read in the binning table
    Bintable.list <- Read_bintable(Filename = BinTable, read.delim = bin.delim, 
        exec = exec, col.index = col.index, chromosomes = ChromosomeList,
        impose.discontinuity = impose.discontinuity)
    Bintable <- Bintable.list[['main.tab']]
    # Create the 0 level directories in the HDF file
    h5createFile(HDF.File)

    for (Folder in Root.folders) {
        CreateGroups(Group.path = Create_Path(Folder), File = HDF.File)
    }
    # Add the chromosome information into the metadata column
    if(!all(ChromosomeList %in% Bintable[,'chr'])){
        stop("All Chromosomes were not listed in the binning table!\n")
    }
    Chrom.lengths <- get_chrom_info(bin.table = Bintable, 
        chrom = ChromosomeList, FUN = length, col.name = 'chr')
    Chrom.sizes <- get_chrom_info(bin.table = Bintable, 
        chrom = ChromosomeList, FUN = max, col.name = 'end')
    Chrom.info.df <- data.frame(chr = names(Chrom.lengths),
        nrow = as.vector(Chrom.lengths),
        size = as.vector(Chrom.sizes),stringsAsFactors = FALSE)
    # Create metadata chromosome groups
    ._Lego_WriteDataFrame_(Lego = HDF.File, Path = c(Root.folders['metadata']), 
        name = Reference.object$metadata.chrom.dataset, object = Chrom.info.df)
    ._Lego_Add_Ranges_(Group.path = Create_Path(c(Root.folders['ranges'],
        Reference.object$hdf.bintable.ranges.group)), Lego = HDF.File, 
        ranges.df = Bintable, name = Reference.object$hdf.ranges.dataset.name, 
        mcol.list = NULL)
    # Create matrices groups
    for (chrom1 in ChromosomeList) {
        CreateGroups(Group.path = Create_Path(c(Root.folders['matrices'],
            chrom1)), File = HDF.File)
        for (chrom2 in ChromosomeList) {
            chr2.path <- Create_Path(c(Root.folders['matrices'],chrom1,chrom2))
            # cat(Chrom.info.df$chr,chrom1,chrom2,"\n")
            CreateGroups(Group.path = chr2.path, File = HDF.File)
            CreateAttributes(Path = chr2.path, File = HDF.File, 
                Attributes = Reference.object$matrices.chrom.attributes,
                data_types = Reference.object$matrices.chrom.attributes.dtype,
                dims = Reference.object$matrices.chrom.attributes.dims,
                maxdims = NULL,
                on = "group")
            Dims <- c(Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"], 
                Chrom.info.df[Chrom.info.df$chr == chrom2,"nrow"])
            if(is.null(ChunkSize)){
                ChunkSize <- ceiling(Dims/100)
            }
            Array.dim <-Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"]
            CreateDataset(Path = c(Root.folders['matrices'],chrom1,chrom2), 
                File = HDF.File, name = Reference.object$hdf.matrix.name, 
                dims = Dims, maxdims = Dims)
        }
    }
    if(CreateCacheObject){
        Lego_track_legos(Lego = HDF.File)
    }
    names(HDF.File) <- Output.Filename
    return(HDF.File)
}

#' Create the entire HDF5 structure and load the bintable from a mcool file
#' 
#' `CreateLego_from_mcool` is a wrapper on CreateLego which creates the Lego 
#' data structure from an mcool file.
#' 
#' mcool are a standard 4D nucleome data structure for Hi-C data. Read more
#' about the 4D nucleome project \href{https://data.4dnucleome.org/}{here}. 
#' 
#' @param chrs \strong{Optional}. 
#' If provided will only create a Lego for these 
#' chromosomes (both cis & trans).
#' 
#' @param binsize \strong{Optional}.
#' The binsize to select from an mcool file.
#' 
#' @inheritParams Lego_load_data_from_mcool
#' 
#' @inheritParams CreateLego
#' 
#' @examples
#' 
#' \dontrun{
#' require(curl)
#' curl_download(url = paste("https://data.4dnucleome.org/"
#' "files-processed/4DNFI7JNCNFB/"
#' "@@download/4DNFI7JNCNFB.mcool",sep = ""),
#' destfile = "./H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' 
#' Output.lego <- paste("./H1-hESC-HiC-4DNFI7JNCNFB-10000",
#' "ICE-normalised-chr1.lego",sep = "-")
#' mcool <- "./H1-hESC-HiC-4DNFI7JNCNFB.mcool"
#' 
#' CreateLego_from_mcool(Lego = Output.lego, 
#' mcool = mcool, 
#' binsize = 10000, 
#' chrs = "chr1")
#' 
#' }
#' 
#' @seealso \code{\link{Lego_load_data_from_mcool}} to load data from 
#' the mcool to a Lego store.
#' 
#' 
CreateLego_from_mcool <- function(Lego = NULL, mcool = NULL, binsize = NULL, 
    chrs = NULL, remove.existing = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(is.null(mcool)){
        stop("mcool must be provided as mcool= /path/to/something")
    }
    if(!file.exists(mcool)){
        stop("mcool not found!")    
    }
    resolutions <- Lego_list_mcool_resolutions(mcool = mcool)
    mcool.version <- GetAttributes(Path = NULL, File=mcool, 
        Attributes="format-version", on = "file", 
        ignore.fun.cast = TRUE)[,"format-version"]
    if(!is.null(resolutions)){
        if(is.null(binsize)){
            stop("binsize cannot be NULL when resolutions are present..\n")
        }
        if(length(binsize) > 1){
            stop("binsize cannot have more than one value\n")
        }
        if(!(as.character(binsize) %in% resolutions)){
            stop("all binsizes were not found in this mcool file. See all",
                " resolutions available with Lego_list_mcool_resolutions\n")
        }
    }
    cooler.remap.chrom <- ._mcool_remap_chromosomes(File = mcool, 
        mcool.version = mcool.version, resolution = !is.null(resolutions), 
        binsize = binsize)
    ChromNames <- cooler.remap.chrom[,"chr.name"]
    ChromNames <- ifelse(!is.null(chrs),ChromNames[ChromNames %in% chrs],
        ChromNames)
    if(!is.null(chrs)){
        if(any(!(chrs %in% ChromNames))){
            stop("Some chrs were not found in this mcool file.\n")
        }
        ChromNames <- ChromNames[ChromNames %in% chrs]
    }
    mcool_bintable_ranges <- ._mcool_bintable_ranges(mcool.file = mcool, 
        resolution = !is.null(resolutions), 
        mcool.remap.chrom = cooler.remap.chrom, binsize = binsize, 
        mcool.version = mcool.version)
    mcool_bintable_ranges <- mcool_bintable_ranges[
    mcool_bintable_ranges[,"chr"] %in% ChromNames,]
    RetVar <- CreateLego(ChromNames=ChromNames, BinTable=mcool_bintable_ranges, 
        Output.Filename=Lego, remove.existing = remove.existing)
    return(RetVar)
}

#' Get all available normalisations in an mcool file.
#' 
#' `Lego_list_mcool_resolutions` lists all available resolutions in the mcool
#' file.
#' 
#' @param mcool \strong{Required}.
#' A parameter specifying the name of an mcool file
#' 
#' @return A named vector listing all possible resolutions in the file. 
#' 
Lego_list_mcool_resolutions <- function(mcool = NULL){
    return(mcool_list_resolutions(mcool = mcool))
}

#' Get all available normalisations in an mcool file.
#' 
#' `Lego_list_mcool_normalisations` lists the names available for 
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
#' Lego_list_mcool_normalisations()
#' 
Lego_list_mcool_normalisations <- function(names.only = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(names.only){
        return(names(Reference.object$mcool.available.normalisations()))
    }
    return(Reference.object$mcool.available.normalisations())
}

#' Check if a normalisation exists in an mcool file.
#' 
#' `Lego_mcool_normalisation_exists` checks if a particular normalisation 
#' exists in an mcool file.
#' 
#' @inheritParams Lego_load_data_from_mcool
#' 
#' @return A boolean vector of length 1 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' require(curl)
#' curl_download(url = paste("https://data.4dnucleome.org/"
#' "files-processed/4DNFI7JNCNFB/"
#' "@@download/4DNFI7JNCNFB.mcool",sep = ""),
#' destfile = "./H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' 
#' mcool <- "./H1-hESC-HiC-4DNFI7JNCNFB.mcool"
#' Lego_mcool_normalisation_exists(mcool = mcool, 
#' norm.factor = "Iterative-Correction",
#' binsize = 10000)
#' 
#' }
#' 
Lego_mcool_normalisation_exists <- function(mcool = NULL, norm.factor = NULL, 
    binsize = NULL){
    Reference.object <- GenomicMatrix$new()
    Norm.factors <- Lego_list_mcool_normalisations()
    Norm.factor <- Norm.factors[norm.factor]
    if(length(Norm.factor)!= 1){
        stop("Please check the available norm factors with ",
            "Lego_list_mcool_normalisations.\n")
    }
    names(Norm.factor) <- NULL
    mcool.version <- GetAttributes(Path = NULL, File=mcool, 
        Attributes="format-version", on = "file", 
        ignore.fun.cast = TRUE)[,"format-version"]
    Bintable.keys <- Reference.object$mcool.bintable.keys(
        version = mcool.version)
    Bintable.group <- Bintable.keys[1]
    resolutions <- Lego_list_mcool_resolutions(mcool = mcool)
    if(!is.null(resolutions) & is.null(binsize)){
        stop("binsize must be provided when different resolutions are present",
            " in an mcool file.\n")
    }
    if(!is.null(resolutions) & !(binsize %in% resolutions)){
        stop("binsize not found in mcool file. Please check available binsizes",
            " with Lego_list_mcool_resolutions.\n")
    }
    if(!is.null(resolutions)){
        Bintable.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Bintable.group))  
    }else{
        Bintable.group.path <- Create_Path(Bintable.group)
    }
    Handler <- ._Lego_Get_Something_(Group.path = Bintable.group.path, 
        Lego = mcool, return.what = "group_handle")
    CloseH5Con(Handle = Handler, type = "group")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    return(Norm.factor %in% GroupList)
}

#' Get the chrominfo for the Hi-C experiment.
#' 
#' `Lego_get_chrominfo` fetches the associated chrominfo table for the 
#' Lego it is associated to. 
#' 
#' @param Lego \strong{Required}.
#' A string specifying the path to the Lego store created with CreateLego.
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
#' Lego.file = system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_chrominfo(Lego = Lego.file)
#' 
Lego_get_chrominfo <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Dataset <- ._Lego_Get_Something_(
        Group.path = Reference.object$hdf.metadata.root, 
        Lego = Lego, Name = Reference.object$metadata.chrom.dataset, 
        return.what = "data")
    return(Dataset)
}

#' Creates a ranges object from provided vectors.
#' 
#' `Lego_make_ranges` creates a GRanges object from the provided arguments
#'  
#' @param Chrom \strong{Required}.
#' A 1 dimensional character vector of size N specifying the chromosomes in the
#' ranges.
#' 
#' @param Start \strong{Required}.
#' A 1 dimensional numeric vector of size N specifying the start positions in 
#' the ranges.
#' 
#' @param End \strong{Required}.
#' A 1 dimensional numeric vector of size N specifying the end positions in 
#' the ranges. Must be less than Start.
#' 
#' @param Strand \strong{Optional}.
#' A 1 dimensional character vector of size N specifying the strand of the 
#' ranges. If not provided, this will be set to the default *.
#' 
#' @param Names \strong{Optional}.
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
#' Test_ranges <- Lego_make_ranges(Chrom = Chrom, Start = Start, End = End)
#' 
Lego_make_ranges = function(Chrom=NULL, Start=NULL, End=NULL, Strand=NULL, 
    Names=NULL){
    Reference.object <- GenomicMatrix$new()

    if(is.null(Names)){
        Names<-paste(Chrom,as.integer(Start),as.integer(End),
            sep=Reference.object$Ranges.separator)
    }
    if(is.null(Strand)){
        Strand<-rep("*",length(Chrom))
    }
    Object<-GenomicRanges::GRanges(
        seqnames=Rle(Chrom),
        ranges=IRanges(Start,end=End,names=Names),
        strand=Rle(strand( Strand )))
    return(Object)
}

#' Store a ranges object in the Lego store.
#' 
#' `Lego_add_ranges` loads a GRanges object into the Lego store.
#' 
#' With this function it is possible to associate other ranges objects with the
#' Lego store. If metadata columns are present, the are also loaded into the 
#' Lego store. Although not explicitly asked for, the metadata columns should
#' not be of type list as this may create complications down the line. We 
#' ask for ranges objects, so if the same ranges object is later retrieved
#' two additional columns will be present. These are the strand and width
#' columns that are obtained when a ranges is converted into a data.frame. 
#' Users can ignore these columns.
#'  
#' @param ranges \strong{Required}.
#' An object of class ranges specifying the ranges to store in the Lego.
#' 
#' @param rangekey \strong{Required}.
#' The name to use for the ranges within the Lego store.
#' 
#' @param remove.existing \strong{Optional}. TRUE
#' Will remove an existing Ranges by the same name and introduce the new one.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @examples
#' 
#' All_Legos <- Lego_list_tracked_legos()
#' Lego_file <- file.path(getwd(),"test.hdf")
#' Lego.file <- All_Legos[Lego_file]
#' Chrom <- c("chrS","chrS","chrS","chrS","chrS")
#' Start <- c(10000,20000,40000,50000,60000)
#' End <- c(10001,20001,40001,50001,60001)
#' Test_ranges <- Lego_make_ranges(Chrom = Chrom, Start = Start, End = End)
#' Lego_add_ranges(Lego = Lego.file, ranges = Test_ranges, 
#' rangekey = "test_ranges")
#' 
#' \dontrun{
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Chrom <- c("chrS","chrS","chrS","chrS","chrS")
#' Start <- c(10000,20000,40000,50000,60000)
#' End <- c(10001,20001,40001,50001,60001)
#' Test_ranges <- Lego_make_ranges(Chrom = Chrom, Start = Start, End = End)
#' Lego_add_ranges(Lego = Lego.file, ranges = Test_ranges, 
#' rangekey = "test_ranges", remove.existing = TRUE)
#' 
#' }
#' 
Lego_add_ranges = function(Lego = NULL, ranges = NULL, rangekey = NULL, 
    remove.existing = TRUE){
    Reference.object <- GenomicMatrix$new()
    if(!(class(ranges) %in% "GRanges") | ("list" %in% class(ranges))){
        stop("Object of class Ranges expected")
    }
    Ranges.df <- as.data.frame(ranges)
    Which.factor <- which(vapply(seq_len(ncol(Ranges.df)), function(x){
            is.factor(Ranges.df[,x])
        }, TRUE))
    Ranges.df[,Which.factor] <- vapply(Which.factor,function(x){
        as.character(Ranges.df[,x])
    },rep("a",nrow(Ranges.df)))
    if(Lego_rangekey_exists(Lego = Lego, rangekey = rangekey)){
        # if(!remove.existing){
            stop("rangekey already exists! Cannot proceed further! ",
                "Please read the documentation to understand Why.\n")
        # }
    }
    Metadata.Cols <- names(Ranges.df)[c(4:ncol(Ranges.df))]
    Metadata.list <- lapply(Metadata.Cols,function(x){
        if(is.factor(Ranges.df[,x])){
            as.character(Ranges.df[,x])
        }else{
            Ranges.df[,x]
        } 
    })
    names(Metadata.list) <- Metadata.Cols
    Ranges.df.coords <- Ranges.df[,c(1,2,3)]
    colnames(Ranges.df.coords) <- Reference.object$NonStrandedColNames
    ._Lego_Add_Ranges_(
        Group.path = Create_Path(c(Reference.object$hdf.ranges.root,rangekey)), 
        Lego = Lego, ranges.df = Ranges.df.coords, name = rangekey, 
        mcol.list = Metadata.list)
    return(TRUE)
}

#' List the ranges tables stored within the Lego.
#' 
#' `Lego_list_rangekeys` lists the names of all ranges associated to a Lego.
#'  
#' @return A one dimensional character vector of length x specifying the names
#' of all ranges currently present in the file.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @examples
#' 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_list_rangekeys(Lego = Lego.file)
#' 
Lego_list_rangekeys = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Handler <- ._Lego_Get_Something_(
        Group.path = Create_Path(Reference.object$hdf.ranges.root), 
        Lego = Lego, return.what = "group_handle")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    return(GroupList)
}

#' Check to see if the Lego contains a ranges with a certain name.
#' 
#' `Lego_rangekey_exists` checks for the presence of a particular ranges with
#' a certain name.
#' 
#' @param rangekey \strong{Required}.
#' A string specifying the name of the ranges to check for.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @return A logical vector of length 1 with either TRUE or FALSE values. 
#' 
#' @examples
#' 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_rangekey_exists(Lego = Lego.file, rangekey = "Bintable")
#' 
Lego_rangekey_exists = function(Lego = NULL, rangekey = NULL){
    Keys <- Lego_list_rangekeys(Lego = Lego)
    return(rangekey %in% Keys)
}

#' Find out what metadata columns are associated to a ranges with a certain 
#' name
#' 
#' `Lego_list_ranges_mcols` will list the metadata columns of the specified
#' ranges if it is present in the Lego store. 
#' 
#' @param rangekey \strong{Optional}.
#' A string specifying the name of the ranges. If not present, the metadata
#' columns of all ranges will be listed.
#'  
#' @inheritParams Lego_get_chrominfo
#' 
#' @return if no metadata columns are present, NA. If metadata columns are
#' present, a data.frame object containing the name of the ranges and the 
#' associated metadata column name.
#' 
#' @examples 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_list_ranges_mcols(Lego = Lego.file, rangekey = "test_ranges")
#' 
Lego_list_ranges_mcols = function(Lego = NULL, rangekey = NULL){
    Reference.object <- GenomicMatrix$new()
    RangeKeys <- Lego_list_rangekeys(Lego = Lego)
    if(!is.null(rangekey)){
        if(!Lego_rangekey_exists(Lego = Lego, rangekey = rangekey)){
            stop("rangekey not found!")
        }
        RangeKeys <- RangeKeys[RangeKeys %in% rangekey]
    }
    mcol.list <- lapply(RangeKeys,function(x){
        Handler <- ._Lego_Get_Something_(
            Group.path = Create_Path(c(Reference.object$hdf.ranges.root, x)), 
            Lego = Lego, return.what = "group_handle")
        GroupList <- h5ls(Handler, 
            datasetinfo = FALSE, 
            recursive = FALSE)[,"name"]
        H5Gclose(Handler)
        data.frame(rangekey = x, 
            m.col = GroupList)
    })
    mcol.df <- do.call(rbind,mcol.list)
    hdf.ranges.protected.names <- Reference.object$hdf.ranges.protected.names()
    mcol.df.filter <- !(mcol.df$m.col %in% hdf.ranges.protected.names)
    mcol.df <- mcol.df[mcol.df.filter,]
    if(nrow(mcol.df)==0){
        mcol.df <- NA
    }
    return(mcol.df)
}

#' List the matrix pairs present in the Lego store.
#' 
#' `Lego_list_matrices` will list all chromosomal pair matrices from the Lego
#' store, with their associated filename, value range, done status and sparse 
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @return Returns a data.frame object with columns chr1, chr2 corresponding 
#' to chromosome pairs, and the associated attributes. filename corresponds to
#' the name of the file that was loaded for the pair. min and max specify the
#' minimum and maximum values in the matrix, done is a logical value 
#' specifying if a matrix has been loaded and sparsity specifies if a matrix
#' is defined as a sparse matrix.
#' 
#' @examples 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_list_matrices(Lego = Lego.file)
#' 
Lego_list_matrices = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    chr1.list <- lapply(ChromInfo[,"chr"], function(chr1){
        chr2.list <- lapply(ChromInfo[,"chr"], function(chr2){
            Colnames <- Reference.object$matrices.chrom.attributes
            Values <- GetAttributes(
                Path = Create_Path(
                    c(Reference.object$hdf.matrices.root, 
                        chr1, 
                        chr2)
                    ),
                File = Lego, 
                Attributes = Colnames, 
                on = "group")
            temp.df <- cbind(data.frame("chr1" = chr1,"chr2" = chr2),Values)
        })
        chr2.df <- do.call(rbind,chr2.list)
    })
    Matrix.list.df <- do.call(rbind,chr1.list)
    rownames(Matrix.list.df) <- NULL
    return(Matrix.list.df)
}

#' Fetch the ranges associated to a rangekey or chromosome.
#' 
#' `Lego_get_ranges` will get a ranges object if present in the Lego store and
#' return a GRanges object.
#' 
#' If a rangekey is present, the ranges will be retrieve and a GRanges 
#' constructed. Metadata columns will also be added. If these are rangekeys
#' other than "Bintable", and had been added using Lego_add_ranges the width 
#' and Strand columns may appear as metadata columns. These will most likely 
#' be artifacts from converting the original ranges object to a data.frame.
#' 
#' @inheritParams Lego_get_chrominfo
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
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_ranges(Lego = Lego.file, chr = "chr19", rangekey = "Bintable")
#' 
Lego_get_ranges = function(Lego = NULL, chr = NULL, rangekey = NULL){
    Reference.object <- GenomicMatrix$new()
    transformlist <- list("strand" = strand, "seqlevels" = seqlevels, 
        "seqlengths" = seqlengths, "width" = width)
    if(is.null(rangekey) | is.null(Lego)){
        stop("rangekey and Lego cannot remain empty!\n")
    }
    if(!Lego_rangekey_exists(Lego = Lego, rangekey = rangekey)){
        stop("rangekey not found!")
    }
    if(length(chr)  > 1){
        stop("chr is expected to be of length 1!\n")   
    }
    Start <- NULL
    Stride <- NULL
    Count <- NULL
    Index <- NULL
    if(!is.null(chr)){
        chromosomes <- ._Lego_Get_Something_(Group.path = Create_Path(
            c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, 
        Name = Reference.object$hdf.ranges.chr.name, 
        return.what = "data")
        if(any(!(chr %in% chromosomes))){
            stop("chr not found!")
        }
        Starts <- ._Lego_Get_Something_(
            Group.path = Create_Path(
                c(Reference.object$hdf.ranges.root, 
                    rangekey)
                ),
        Lego = Lego, 
        Name = Reference.object$hdf.ranges.offset.name, 
        return.what = "data")
        Lengths <- ._Lego_Get_Something_(Group.path = Create_Path(
            c(Reference.object$hdf.ranges.root, 
                rangekey)
            ),
        Lego = Lego, 
        Name = Reference.object$hdf.ranges.lengths.name, 
        return.what = "data")
        Which.one <- chromosomes == chr
        Start <- Starts[Which.one]
        Stride <- 1
        Count <- Lengths[Which.one]
    }
    Dataset <- ._Lego_Get_Something_(Group.path = Create_Path(
        c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.dataset.name, 
        Start = Start, Stride = Stride,
        Count = Count, return.what = "data")
    Dataset <- Lego_make_ranges(Chrom = Dataset[,'chr'], 
        Start = Dataset[,'start'], End = Dataset[,'end'])
    MCols <- Lego_list_ranges_mcols(Lego = Lego, rangekey = rangekey)
    if(is.data.frame(MCols)){
        Fltr <- MCols$m.col %in% 
        Reference.object$genomic.ranges.protected.names
        GRangesCols <- MCols[Fltr,]
        MCols.col <- as.character(MCols[!Fltr,"m.col"])
        m.start <- Start[1]
        m.stride <- Stride[1]
        m.count <- Count[1]
        if(length(GRangesCols) > 0){
            genomic.ranges.FUN.names <- names(transformlist)
            FUN.names <- genomic.ranges.FUN.names[genomic.ranges.FUN.names 
            %in% GRangesCols]
            for(x in FUN.names){
                FUN <- transformlist[[x]]
                FUN(Dataset) <- ._Lego_Get_Something_(Group.path = Create_Path(
                    c(Reference.object$hdf.ranges.root, rangekey)),
                Lego = Lego, Name = x, Start = m.start, Stride = m.stride,
                Count = m.count, return.what = "data")
            }
        }
        if(length(MCols.col) > 0){
            MCols.DF.list <- lapply(MCols.col,function(x){
                Dataset <- ._Lego_Get_Something_(
                    Group.path = Create_Path(c(
                        Reference.object$hdf.ranges.root, 
                        rangekey)), Lego = Lego, Name = x, 
                    Start = m.start, Stride = m.stride,
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
#' `Lego_get_bintable` makes a call to \code{\link{Lego_get_ranges}} to 
#' retrieve the binning table of the associated Lego store. This is equivalent
#' to passing the argument rangekey = "bintable" in 
#' \code{\link{Lego_get_ranges}}
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @param chr \strong{Optional}.
#' A chr string specifying the chromosome to select from the ranges.
#' 
#' @examples 
#' 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_bintable(Lego = Lego.file)
#' 
#' @seealso Lego_get_ranges
Lego_get_bintable = function(Lego = NULL, chr = NULL){
    Reference.object <- GenomicMatrix$new()
    Table <- Lego_get_ranges(Lego = Lego, chr = chr, 
        rangekey = Reference.object$hdf.bintable.ranges.group)
    return(Table)
}

#' Returns the position of the supplied ranges in the binning table 
#' associated to the Hi-C experiment.
#' 
#' `Lego_fetch_range_index` constructs a ranges object using 
#' \code{\link{Lego_make_ranges}}, creates an overlap operation using 
#' \code{\link[GenomicRanges]{findOverlaps}}, where the constructed ranges is
#' the \emph{subject} and the Hi-C experiment associated binning table is the 
#' \emph{query}. The return of this object is a list of ranges with their 
#' corresponding indices in the binning table.
#'  
#' @inheritParams Lego_get_chrominfo
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
#' Indexes is a column of class \code{\link[IRanges]{IntegerList}}, which is 
#' part of the larger \code{\link[IRanges]{AtomicList}} superset. This 
#' "Indexes" column can be accessed like a normal GRanges column with the 
#' additional list accessor [[]] in place of the normal vector accessor [].
#' 
#' @examples 
#' 
#' Chrom <- c("chr19","chr19")
#' Start <- c(1,40000)
#' End <- c(1000000,2000000)
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Test_Run <- Lego_fetch_range_index(Lego = Lego.file, chr = Chrom, 
#' start = Start, end = End)
#' Test_Run$Indexes[[1]]
#' 
Lego_fetch_range_index = function(Lego = NULL, chr = NULL, start = NULL, 
    end = NULL, names = NULL, type = "any"){
    AllTypes<-c("any","within")
    if( any(!(type %in% AllTypes)) ){
        stop("type takes one of two arguments: c(\"any\",\"within\")")
    }
    if(is.null(chr) | is.null(start) | is.null(end) | is.null(Lego)){
        stop("Chrom, start, end and Lego cannot be empty\n")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    if(any(!(chr %in% ChromInfo[,'chr']))){
        stop("Provided chr does not exist in chromosome list.\n")
    }
    if(any(!is.character(chr) | !._Check_numeric(start) | 
        !._Check_numeric(end))){
        stop("Provided chr, start, end do not match expected class ",
            "definitions of character, numeric, numeric\n")
    }
    Unique.chromosomes <- unique(chr)
    OverlapByChromosome.list <- lapply(Unique.chromosomes,function(cur.chr){
        Filter <- chr==cur.chr
        Cur.Chrom <- chr[Filter]
        Cur.Start <- start[Filter]
        Cur.end <- end[Filter]
        Cur.Names <- names[Filter]
        SubjectRanges <- Lego_get_bintable(Lego = Lego, chr = cur.chr)
        if( any(!(Cur.end <= max(end(SubjectRanges)) & 
            Cur.Start >= min(start(SubjectRanges)))) ){
            stop("Start or end is out of ranges for Bintable\n")
        }
        QueryRanges <- Lego_make_ranges(Chrom=Cur.Chrom, Start=Cur.Start, 
            End=Cur.end, Names=Cur.Names)
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
#' `Lego_return_region_position` takes as input a human-readable coordinate
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
#' @inheritParams Lego_get_chrominfo
#' 
#' @param region \strong{Required}.
#' A character vector of length 1 specifying the region to overlap. It must 
#' take the form chr:start:end. 
#' 
#' @return Returns a 1 dimensional vector containing the position of the
#' overlapping regions in the bintable associated the Lego store.
#' 
#' @examples
#' 
#' Coordinate <- "chr19:1:1000000"
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Test_Run <- Lego_return_region_position(Lego = Lego.file, 
#' region = Coordinate)
#' 
#' \dontrun{
#' 
#' Coordinate <- c("chr19:1:1000000","chr19:40000:2000000")
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Test_Run <- Lego_return_region_position(Lego = Lego.file, 
#' region = Coordinate)
#' 
#' This will generate an error because the module itself expects as input a
#' vector of length 1. 
#' 
#' }
Lego_return_region_position = function(Lego = NULL, region=NULL){
    if(!is.character(region) | length(region) > 1){
        stop("region must be a character vector of length 1")
    }
    Coord.Split<- Split_genomic_coordinates(Coordinate=region)
    region.chr<-Coord.Split[[1]][1]
    region.start<-as.numeric(Coord.Split[[1]][2])
    region.stop<-as.numeric(Coord.Split[[1]][3])
    region.ranges<-Lego_fetch_range_index(Lego = Lego, chr=region.chr, 
        start=region.start, end=region.stop, type="within")
    Vector.coordinates <- region.ranges$Indexes[[1]]
    return(Vector.coordinates)
}

#' Load a NxM dimensional matrix into the Lego store.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @param chr1 \strong{Required}.
#' A character vector of length 1 specifying the chromosome corresponding to 
#' the rows of the matrix
#' 
#' @param chr2 \strong{Required}.
#' A character vector of length 1 specifying the chromosome corresponding to 
#' the columns of the matrix
#' 
#' @param matrix.file \strong{Required}.
#' A character vector of length 1 specifying the name of the file to load as a
#' matrix into the Lego store. 
#' 
#' @param exec \strong{Required}.
#' A string specifying the program to use for reading the file. Use cat for txt
#' files, for bz2 files use bzcat and for gz files zcat.
#' 
#' @param delim \strong{Optional}. Default " "
#' The delimiter of the matrix file.
#' 
#' @param remove.prior \strong{Optional}. Default FALSE
#' If a matrix was loaded before, it will not be replaced. Use remove.prior to
#' override and replace the existing matrix.
#' 
#' @param num.rows \strong{Optional}. Default 2000
#' Number of rows to read, in each chunk.
#' 
#' @param is.sparse \strong{Optional}. Default FALSE
#' If true, designates the matrix as being a sparse matrix, and computes the
#' sparsity.index. The sparsity index measures the proportion of non-zero rows
#' or columns at a certain distance from the diagonal (100) in cis interaction
#' matrices.
#'
#' @param sparsity.bins \strong{Optional}. Default 100
#' With regards to computing the sparsity.index, this parameter decides the 
#' number of bins to scan from the diagonal. 
#' 
#' @examples
#' 
#' Test.mat <- matrix(NA,nrow = 800, ncol = 800)
#' Row <- row(Test.mat)
#' Col <- col(Test.mat)
#' Dist <- Col - Row
#' Matrix.file <- "Test_matrix.txt"
#' write.table(x = Dist, file = Matrix.file, sep = " ", quote = FALSE, 
#' row.names = FALSE, col.names = FALSE)
#' All_Legos <- Lego_list_tracked_legos()
#' Lego_file <- file.path(getwd(),"test.hdf")
#' Lego.file <- All_Legos[Lego_file]
#' Lego_load_matrix(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19",
#' matrix.file = Matrix.file, delim = " ", exec = "cat", 
#' remove.prior = TRUE)
#' 
Lego_load_matrix = function(Lego = NULL, chr1 = NULL, chr2 = NULL, 
    matrix.file = NULL, delim = " ", exec = NULL, remove.prior = FALSE, 
    num.rows = 2000, is.sparse = FALSE, sparsity.bins = 100){

    Reference.object <- GenomicMatrix$new()
    ListVars <- list(Lego = Lego, chr1 = chr1, chr2 = chr2, 
        matrix.file = matrix.file, is.sparse = is.sparse, 
        sparsity.bins = sparsity.bins, exec = exec, delim = delim, 
        distance = distance, remove.prior = remove.prior)
    lapply(seq_along(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    lapply(seq_along(ListVars[c("Lego","chr1","chr2","matrix.file","exec")]),
        function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    if(Lego_matrix_isdone(Lego = Lego, chr1 = chr1,
        chr2 = chr2) && !remove.prior){
        stop("A matrix was preloaded before. ",
            "Use remove.prior = TRUE to force value replacement\n")
    }
    Chrom.info.df <- Lego_get_chrominfo(Lego = Lego)
    chr1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr1,"nrow"]
    chr2.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr2,"nrow"]
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,chr1,chr2))
    compute.sparsity <- FALSE
    if(is.sparse && chr1 == chr2){
        compute.sparsity <- TRUE
    }
    RetVar <- ._ProcessMatrix_(Lego = Lego, Matrix.file = matrix.file, 
        delim = delim, exec = exec, Group.path = Group.path, 
        chr1.len = chr1.len, chr2.len = chr2.len, num.rows = num.rows,
        is.sparse = is.sparse, compute.sparsity = compute.sparsity, 
        sparsity.bins = sparsity.bins)
    return(RetVar)
}


#' Load a NxN dimensional sub-distance \emph{cis} matrix into 
#' the Lego store.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
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
#' @param num.rows \strong{Optional}. Default 2000
#' Number of rows to insert per write operation in the HDF file.
#' 
#' @examples
#' 
#' Test.mat <- matrix(NA,nrow = 800, ncol = 800)
#' Row <- row(Test.mat)
#' Col <- col(Test.mat)
#' Dist <- Col - Row
#' Matrix.file <- "Test_matrix.txt"
#' write.table(x = Dist, file = Matrix.file, sep = " ", quote = FALSE, 
#' row.names = FALSE, col.names = FALSE)
#' All_Legos <- Lego_list_tracked_legos()
#' Lego_file <- file.path(getwd(),"test.hdf")
#' Lego.file <- All_Legos[Lego_file]
#' Lego_load_cis_matrix_till_distance(Lego = Lego.file, chr = "chr19", 
#' matrix.file = Matrix.file, delim = " ", distance = 200, remove.prior = TRUE)
#' 
Lego_load_cis_matrix_till_distance = function(Lego = NULL, chr = NULL, 
    matrix.file = NULL, delim = " ", distance = NULL, remove.prior = FALSE, 
    num.rows = 2000, is.sparse = FALSE, sparsity.bins = 100){

    Reference.object <- GenomicMatrix$new()
    ListVars <- list(Lego = Lego, chr = chr, matrix.file = matrix.file, 
        is.sparse = is.sparse, sparsity.bins = sparsity.bins, delim = delim, 
        distance = distance, remove.prior = remove.prior)
    lapply(seq_along(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    lapply(seq_along(ListVars[c("Lego","chr","matrix.file","distance")]),
        function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr, 
        chr2 = chr)){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    if(Lego_matrix_isdone(Lego = Lego, chr1 = chr, 
        chr2 = chr) && !remove.prior){
        stop("A matrix was preloaded before. Use remove.prior = TRUE to ",
            "force value replacement\n")
    }
    Chrom.info.df <- Lego_get_chrominfo(Lego = Lego)
    chr1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr,"nrow"]
    chr2.len <- chr1.len
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,chr,chr))
    compute.sparsity <- FALSE
    if(is.sparse){
        compute.sparsity <- TRUE
    }
    RetVar <- ._Process_matrix_by_distance(Lego = Lego, 
        Matrix.file = matrix.file, delim = delim, Group.path = Group.path, 
        chr1.len = chr1.len, chr2.len = chr2.len, num.rows = num.rows, 
        distance = distance, is.sparse = is.sparse, 
        compute.sparsity = compute.sparsity, sparsity.bins = sparsity.bins)
    return(RetVar)
}


#' Load a NxN dimensional matrix into the Lego store from an mcool file.
#' 
#' Read an mcool contact matrix coming out of 4D nucleome projects into a 
#' Lego store. 
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @inheritParams CreateLego_from_mcool
#' @param mcool \strong{Required}.
#' Path to an mcool file. 
#' 
#' @param dont.look.for.chr2 \strong{Required}.
#' At startup, the function will attempt to search for the first occurence 
#' of a chr2 contact value. This is done to avoid the reading of all chr1
#' values for every chunk processed. If chr1 and chr2 are equivalent, consider
#' setting it to FALSE.
#' 
#' @param norm.factor \strong{Optional}. Default "Iterative-Correction".
#' The normalization factor to use for normalization from an mcool file. 
#' 
#' @param cooler.batch.size \strong{Optional}. Default 1000000.
#' The number of values to read per iteration through a mcool file.
#' 
#' @param matrix.chunk \strong{Optional}. Default 2000.
#' The nxn matrix square to fill per iteration in a mcool file.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' require(curl)
#' curl_download(url = paste("https://data.4dnucleome.org/"
#' "files-processed/4DNFI7JNCNFB/"
#' "@@download/4DNFI7JNCNFB.mcool",sep = ""),
#' destfile = "./H1-hESC-HiC-4DNFI7JNCNFB.mcool")
#' 
#' Output.lego <- paste("./H1-hESC-HiC-4DNFI7JNCNFB-10000",
#' "ICE-normalised-chr1.lego",sep = "-")
#' mcool <- "./H1-hESC-HiC-4DNFI7JNCNFB.mcool"
#' 
#' CreateLego_from_mcool(Lego = Output.lego, 
#' mcool = mcool, 
#' binsize = 10000, 
#' chrs = "chr1")
#' 
#' Lego_load_data_from_mcool(Lego = Output.lego, mcool = mcool, 
#' chr1 = "chr1", chr2 = "chr1", binsize = 10000, 
#' cooler.batch.size = 1000000, matrix.chunk = 2000, 
#' dont.look.for.chr2 = TRUE, remove.prior = TRUE, 
#' norm.factor = "Iterative-Correction")
#' 
#' }
#' 
#' 
#' @seealso \code{\link{CreateLego_from_mcool}} to create matrix from an mcool
#' file, \code{\link{Lego_list_mcool_resolutions}} to list available 
#' resolutions in an mcool file, \code{\link{Lego_list_mcool_normalisations}} 
#' to list available normalisation factors in the mcool file. 
#'
Lego_load_data_from_mcool <- function(Lego = NULL, mcool = NULL, chr1 = NULL, 
    chr2 = NULL, binsize = NULL, cooler.batch.size = 1000000, 
    matrix.chunk = 2000, dont.look.for.chr2 = FALSE, remove.prior = FALSE,
    norm.factor = "Iterative-Correction"){
    Reference.object <- GenomicMatrix$new()
    if(is.null(chr1) | is.null(chr2)){
        stop("chr1, chr2 cannot be NULL.\n")
    }
    if(length(chr1) != length(chr2) | length(chr1) != 1){
        stop("chr1, chr2 are expected to be of length 1.\n")   
    }
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, 
        chr2 = chr2)){
        stop("Provided chromosomes do not exist in the chrominfo table\n")
    }
    if(Lego_matrix_isdone(Lego = Lego, chr1 = chr1, 
        chr2 = chr2) & !remove.prior){
        stop("A matrix was preloaded before. ",
            "Use remove.prior = TRUE to force value replacement\n")
    }
    resolutions <- Lego_list_mcool_resolutions(mcool = mcool)
    if(!is.null(resolutions) & is.null(binsize)){
        stop("binsize must be provided when ",
            "different resolutions are present in an mcool file.\n")
    }
    if(!is.null(resolutions) & !(binsize %in% resolutions)){
        stop("binsize not found in mcool file. ",
            "Please check available binsizes ",
            "with Lego_list_mcool_resolutions.\n")
    }
    if(!is.null(Norm.factor)){
        if(!Lego_mcool_normalisation_exists(mcool = mcool, 
            norm.factor = norm.factor, binsize = binsize)){
            stop(norm.factor," was not found in this mcool file.\n")            
        }
        Norm.factors <- Lego_list_mcool_normalisations()
        Norm.factor <- Norm.factors[norm.factor]
        names(Norm.factor) <- NULL
    }else{
        Norm.factor <- NULL
    }
    RetVar <- ._Process_mcool(Lego = Lego, File = mcool, 
        cooler.batch.size = cooler.batch.size, binsize = binsize,
    matrix.chunk = matrix.chunk, chr1 = chr1, 
    chr2 = chr2, dont.look.for.chr2 = dont.look.for.chr2, 
    norm.factor = Norm.factor, resolution = !is.null(resolutions))
    return(RetVar)
}

#' Check if a matrix has been loaded for a chromosome pair.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns a logical vector of length 1, specifying if a matrix has
#' been loaded or not.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_isdone(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19")
#' 
Lego_matrix_isdone = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_list_matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    return(Matrix.list[Matrix.list$chr1 == chr1 & 
        Matrix.list$chr2 == chr2, "done"])
}

#' Check if a matrix for a chromosome pair is sparse.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns a logical vector of length 1, specifying if a matrix was
#' loaded as a sparse matrix.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_issparse(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19")
#' 
Lego_matrix_issparse = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_list_matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    return(Matrix.list[Matrix.list$chr1 == chr1 & 
        Matrix.list$chr2 == chr2, "sparsity"])
}


#' Get the maximum loaded distance from the diagonal of any matrix.
#' 
#' If values beyond a certain distance were not loaded in the matrix, this
#' distance parameter is useful. This package by default will check this param
#' to make sure that it is not returning non-existent data.
#'  
#' `Lego_matrix_maxdist` will return this parameter.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns an integer vector of length 1, specifying the maximum 
#' distance loaded for that matrix
#'
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_maxdist(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19") 
#' 
Lego_matrix_maxdist = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_list_matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not loaded\n")
    }
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
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns a logical vector of length 1, specifying if the matrix 
#' exists or not.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_exists(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19") 
#' 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_exists(Lego = Lego.file, chr1 = "chr19", chr2 = "chr20")
#' 
Lego_matrix_exists = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    ChromInfo.df <- Lego_get_chrominfo(Lego = Lego)
    return(all(c(chr1,chr2) %in% ChromInfo.df[,"chr"]))
}

#' Return the value range of the matrix
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns a numeric vector of length 2, specifying the minimum and
#' maximum finite real values in the matrix.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_minmax(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19")
#' 
Lego_matrix_minmax = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_list_matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- c(Matrix.list[Filter, "min"],Matrix.list[Filter, "max"])
    return(Extent)
}

#' Return the dimensions of a matrix
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns the dimensions of a Hi-C matrix for any given 
#' chromosome pair.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_dimensions(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19")
#'
Lego_matrix_dimensions = function(Lego=NULL, chr1=NULL, chr2=NULL){
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Reference.object <- GenomicMatrix$new()
    Extents <- ._GetDimensions(group.path = Create_Path(
        c(Reference.object$hdf.matrices.root,chr1,chr2)), 
        dataset.path =Reference.object$hdf.matrix.name, 
        File = Lego, return.what = "size")
    return(Extents)
}


#' Return the filename of the loaded matrix
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @return Returns a character vector of length 1 specifying the filename of 
#' the currently loaded matrix.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_matrix_filename(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19")
#' 
Lego_matrix_filename = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_list_matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- Matrix.list[Filter, "filename"]
    return(Extent)
}

#' Return values separated by a certain distance.
#' 
#' `Lego_get_values_by_distance` can fetch values with or without 
#' transformation or subsetted by a certain distance. Please note, 
#' this module is not an iterable module.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @param chr \strong{Required}.
#' A string specifying the chromosome for the cis Hi-C matrix from which values
#' will be retrieved at a certain distance.
#' 
#' @param distance \strong{Required}. 0 based.
#' Fetch values separated by distance.
#'
#' @param constrain.region \strong{Required}.
#' A character vector of length 1 with the form chr:start:end specifying the 
#' region for which the distance values must be retrieved.
#' 
#' @param batch.size \strong{Optional}. Default 500
#' A numeric vector of length 1 specifying the size of the chunk to retrieve
#' for diagonal selection.
#' 
#' @param FUN \strong{Optional}.
#' If provided a data transformation with FUN will be applied before values
#' are returned.
#' 
#' @return Returns a numeric vector of length N depending on the presence of 
#' constrain.region, FUN and distance from the main diagonal.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_values_by_distance(Lego = Lego.file, chr = "chr19", distance = 0)
#' 
#' Failsafe_median <- function(x){
#'      x[is.nan(x) | is.infinite(x) | is.na(x)] <- 0
#'      return(median(x))
#' }
#' 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_values_by_distance(Lego = Lego.file, chr = "chr19", distance = 4,
#' FUN = Failsafe_median)
#' 
#' @seealso \code{\link{Lego_get_matrix_within_coords}} to get matrix by using 
#' matrix coordinates, \code{\link{Lego_fetch_row_vector}} to get values in a 
#' certain row/col and subset them, \code{\link{Lego_get_vector_values}} to get 
#' values using matrix coordinates, \code{\link{Lego_get_matrix}} to get matrix 
#' by using matrix coordinates.
#' 
Lego_get_values_by_distance = function(Lego = NULL, chr = NULL, 
    distance  = NULL, constrain.region=NULL,batch.size=500,FUN=NULL){
    if(any(vapply(list(Lego,chr,distance),is.null,TRUE))) {
        stop("Lego, chr, distance cannot be NULL.\n")
    }
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr, chr2 = chr) & 
        !Lego_matrix_isdone(Lego = Lego, chr1 = chr, chr2 = chr)){
        stop("Chromosome is not listed or has not been ",
            "loaded in this HDF file.\n")
    }
    if(any(vapply(list(Lego,chr,distance),length,1) > 1)){
        stop("Lego, chr and distance can only have values of length 1.\n")
    }
    Nrow <- (ChromInfo[ChromInfo$chr == chr,"nrow"]-1)
    if(any(distance > Nrow, distance < 0)){
        stop("distance must range between 0 and",Nrow,"\n")
    }
    Max.dist <- Lego_matrix_maxdist(Lego = Lego, chr1 = chr, chr2 = chr)
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
    if(!is.null(constrain.region)){
        Vector.coordinates <- Lego_return_region_position(Lego = Lego, 
            region=constrain.region)
        if(is.null(Vector.coordinates)){
            stop("Overlap operation was unsuccessful! ",
                "Please check coordinates "
                ,constrain.region)
        }
        Vector.start <- min(Vector.coordinates)
        Vector.stop <- max(Vector.coordinates)
    }
    Starting.col <- Vector.start + distance
    Count <- ((Vector.stop - Vector.start) + 1) - distance 
    tot.mat.extracted <- Count * Count
    # cat("Boo ",Count,"\n")
    Start <- c(Vector.start,Starting.col)
    Stride <- c(1,1)
    Counts <- Count
    CumSums <- 0
    Groups <- c(Reference.object$hdf.matrices.root,chr)
    if(Count > batch.size){
        repeated <- floor(Count/batch.size)
        Counts <- rep(batch.size,repeated)
        if(repeated != ceiling(Count/batch.size)){
            cumulative <- sum(Counts)
            Counts <- c(Counts,(Count-cumulative))
        }
        CumSums <- cumsum(c(0,Counts[seq_len(length(Counts)-1)]))
    }
    DistancesVector.list <- lapply(seq_len(length(Counts)),function(x){
        Count <- Counts[x]
        Offset <- CumSums[x]
        cur.start <- Start+Offset
        diag(._Lego_Get_Something_(Group.path = Path, Lego = Lego, 
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
#' `Lego_get_matrix_within_coords` will fetch a matrix subset after 
#' creating an overlap operation between both regions and the bintable 
#' associated to the Lego store. 
#' This function calls \code{\link{Lego_get_matrix}}.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @param x.coords \strong{Required}.
#' A string specifying the region to subset on the rows. It takes the form
#' chr:start:end. An overlap operation with the associated bintable will be 
#' done to identify the bins to subset on the row
#' 
#' @param y.coords \strong{Required}.
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
#' @return Returns a matrix of dimension x.coords binned length by y.coords
#' binned length. This may differ based on FUN.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_matrix_within_coords(Lego = Lego.file, 
#' x.coords = "chr19:40000:2000000", 
#' y.coords = "chr19:40000:2000000")
#'
#' Lego_get_matrix_within_coords(Lego = Lego.file, 
#' x.coords = "chr19:40000:2000000", 
#' y.coords = "chr19:40000:2000000", 
#' FUN = mean)
#' 
#' Lego_get_matrix_within_coords(Lego = Lego.file, 
#' x.coords = "chr19:40000:2000000", 
#' y.coords = "chr19:40000:2000000", 
#' FUN = median)
#'  
#' @seealso \code{\link{Lego_get_matrix}} to get matrix by using matrix 
#' coordinates, \code{\link{Lego_get_values_by_distance}} to get values 
#' separated at a certain distance, \code{\link{Lego_fetch_row_vector}} to get
#' values in a certain row/col and subset them, 
#' \code{\link{Lego_get_vector_values}} to get values using matrix coordinates.
Lego_get_matrix_within_coords = function(Lego = NULL, x.coords=NULL, 
    y.coords=NULL, force = FALSE, FUN=NULL){
    type <- "within"
    if( (is.null(x.coords)) | (is.null(y.coords)) ){
        stop("x.coords, y.coords and Lego cannot be NULL")
    }
    if(!(length(x.coords)==1) | !(length(y.coords)==1)){
        stop("This function processes single process calls at a time. ",
            "Setup an Iterator for more functionality")
    }
    if( !is.character(x.coords) | !is.character(y.coords) ){
        stop("Two string variables were expected for x.coords & y.coords, ",
            "found x.coords class ", class(x.coords), 
            " and y.coords class ", 
            class(y.coords))
    }
    xcoords.split <- Split_genomic_coordinates(Coordinate=x.coords)
    chr1 <- xcoords.split[[1]][1]
    chr1.start <- as.numeric(xcoords.split[[1]][2])
    chr1.stop <- as.numeric(xcoords.split[[1]][3])
    chr2.split <- Split_genomic_coordinates(Coordinate=y.coords)
    chr2 <- chr2.split[[1]][1]
    chr2.start <- as.numeric(chr2.split[[1]][2])
    chr2.stop <- as.numeric(chr2.split[[1]][3])
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    if( any(!(c(chr1,chr2) %in% ChromosomeList)) ){
        stop("Provided chromosomes were not found in chromosome list.")
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop(chr1,chr2," matrix is yet to be loaded into the class.")
    }
    chr1.ranges <- Lego_fetch_range_index(Lego = Lego, chr = chr1,
    start = chr1.start, end = chr1.stop, names=NULL, type=type)
    chr2.ranges <- Lego_fetch_range_index(Lego = Lego, chr = chr2, 
        start = chr2.start, end = chr2.stop, names=NULL, type=type)
    if(is.null(chr1.ranges$Indexes[[1]]) | is.null(chr2.ranges$Indexes[[1]])){
        stop("Overlap operation was unsuccessful! Please check coordinates ",
            x.coords," & ",y.coords)
    }
    x.vector <- chr1.ranges$Indexes[[1]]
    y.vector <- chr2.ranges$Indexes[[1]]
    Matrix <- Lego_get_matrix(Lego = Lego, chr1=chr1, chr2=chr2, 
        x.vector=x.vector, y.vector=y.vector, force = force, FUN=FUN)
    return(Matrix)
}

#' Return a matrix subset.
#' 
#' `Lego_get_matrix` will fetch a matrix subset between row values 
#' ranging from min(x.vector) to max(x.vector) and column values ranging from 
#' min(x.vector) to max(x.vector)
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @param x.vector \strong{Required}.
#' A one-dimensional numeric vector specifying the rows to subset. 
#' 
#' @param y.vector \strong{Required}.
#' A one-dimensional numeric vector specifying the columns to subset. 
#' 
#' @param FUN \strong{Optional}.
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#' 
#' @inheritParams Lego_get_matrix_within_coords
#' 
#' @return Returns a matrix of dimension x.vector length by y.vector length. 
#' This may differ based on the operations with FUN.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_matrix(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19", 
#' x.vector = c(1:10), y.vector = c(1:10))
#' 
#' @seealso \code{\link{Lego_get_matrix_within_coords}} to get matrix by using 
#' matrix genomic coordinates, \code{\link{Lego_get_values_by_distance}} to get 
#' values separated at a certain distance, \code{\link{Lego_fetch_row_vector}} 
#' to getvalues in a certain row/col and subset them, 
#' \code{\link{Lego_get_vector_values}} to get values using matrix coordinates.
Lego_get_matrix = function(Lego = NULL, chr1=NULL, chr2=NULL, x.vector=NULL,
    y.vector=NULL, force = FALSE, FUN=NULL){
    # cat(" Rows: ",x.vector," Cols: ",y.vector,"\n")
    if(any(!._Check_numeric(x.vector) | !._Check_numeric(y.vector))){
        stop("x.vector and y.vector must be numeric.\n")
    }
    if(is.null(chr1) | is.null(chr2) | is.null(x.vector) | is.null(y.vector)){
        stop("Either of chr1, chr2, x.vector or y.vector ",
            "were provided as NULL values.\n")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    if( any(!(c(chr1,chr2) %in% ChromosomeList)) ){
        stop("Provided chromosomes were not found in chromosome list.\n")
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop(chr1,chr2," matrix is yet to be loaded into the class.\n")
    }
    chr1.len <- ChromInfo$nrow[ChromInfo$chr == chr1]
    chr2.len <- ChromInfo$nrow[ChromInfo$chr == chr2]
    if(any(x.vector > chr1.len) | 
        any(y.vector > chr2.len) | 
        min(x.vector,y.vector) < 1 ) {
        stop("x.vector or y.vector falls outside ",
            "the bounds of loaded Bintables") 
    }
    Matrix <- Lego_get_vector_values(Lego = Lego, chr1=chr1, chr2=chr2, 
        xaxis=x.vector, yaxis=y.vector, force = force)
    if(is.null(FUN)){
        return(Matrix)              
    }else{
        return(FUN(Matrix))
    }
}

#' Return row or col vectors.
#' 
#' `Lego_fetch_row_vector` will fetch any given rows from a matrix. If 
#' required, the rows can be subsetted on the columns and transformations 
#' applied. Vice versa is also true, wherein columns can be retrieved and 
#' rows subsetted.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
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
#' @inheritParams Lego_get_matrix_within_coords
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
#' Coordinate <- c("chr19:1:40000","chr19:40001:80000")
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Test_Run <- Lego_fetch_row_vector(Lego = Lego.file, 
#' chr1 = "chr19", chr2 = "chr19", by = "ranges", vector = Coordinate,
#' regions = c("chr19:1:1000000", "chr19:40001:2000000"))
#' 
#' @seealso \code{\link{Lego_get_matrix_within_coords}} to get matrix by using 
#' matrix genomic coordinates, \code{\link{Lego_get_values_by_distance}} to get 
#' values separated at a certain distance, \code{\link{Lego_fetch_row_vector}} 
#' to get values in a certain row/col and subset them, 
#' \code{\link{Lego_get_matrix}} to get matrix by using matrix coordinates
#' 
Lego_fetch_row_vector = function(Lego = NULL, chr1=NULL, chr2=NULL, by=NULL, 
    vector=NULL, regions=NULL, force = FALSE, flip = FALSE, FUN=NULL){
    Chrom.all <- c(chr1,chr2)
    if(is.null(chr1) | is.null(chr2) | is.null(by) | is.null(vector)) {
        stop("Either of chr1, chr2, by or vector were provided as NULL values")
    }
    if(!(by %in% c("position","ranges")) | length(by) != 1){
        stop("by expects a vector of type character, length 1 and ",
            "takes either one of position or ranges as values")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    max.dist <- Lego_matrix_maxdist(Lego = Lego, chr1 = chr1, chr2 = chr2)
    if(!is.character(chr1) | !is.character(chr2)){
        stop("Provided Chromosomes does not appear to be of class character")
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
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
        stop("vector must be of class",
            ifelse(length(Class.exp)>1,paste(Class.exp,collapse=" or "),
                paste(Class.exp)),"when by has value ",by)
    }
    if(length(Chrom.all)>2){
        stop("This module is not iterable")
    }
    chr1.ranges <- Lego_get_bintable(Lego = Lego, chr=chr1)
    chr2.ranges <- Lego_get_bintable(Lego = Lego, chr=chr2)
    if(flip){
        chr1.ranges <- Lego_get_bintable(Lego = Lego, chr=chr2)
        chr2.ranges <- Lego_get_bintable(Lego = Lego, chr=chr1)
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
            y <- Lego_return_region_position(Lego = Lego, region = region)
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
        Values <- Lego_get_vector_values(Lego = Lego, chr1=chr1, chr2=chr2, 
            xaxis=x, yaxis=y, FUN=FUN, force = force)
        return(Values)
    }) 
    return(Vector.values)
}
#' Return a N dimensional vector selection.
#' 
#' `Lego_get_vector_values` is the base function being used by all 
#' other matrix retrieval functions.
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
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
#' @inheritParams Lego_get_matrix_within_coords
#' 
#' @return Returns a vector of length yaxis if length of xaxis is 1. Else
#' returns a matrix of dimension xaxis length by yaxis length.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_matrix(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19", 
#' x.vector = c(1:10), y.vector = c(1:10))
#' 
#' @section Note: Whatever the length of xaxis or yaxis may be, the coordinates 
#' under consideration will range from min(xaxis) to max(xaxis) on the rows or
#' min(yaxis) to max(yaxis) on the columns.
#' 
Lego_get_vector_values = function(Lego = NULL, chr1=NULL, chr2=NULL, xaxis=NULL,
    yaxis=NULL, FUN=NULL, force = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(is.null(chr1) | is.null(chr2)){
        stop("chr1 and chr2 keys cannot be empty!")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    if(!is.character(chr1) | !is.character(chr2)){
        stop("Provided Chromosomes does not appear to be of class character")
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop(chr1,chr2," matrix is yet to be loaded.")
    }
    if(is.null(xaxis) & is.null(yaxis)){
        stop("Both xaxis and yaxis cannot be null")
    }
    Start <- c(min(xaxis),min(yaxis))
    Stride <- c(1,1)
    Count <- c(length(xaxis),length(yaxis))
    Max.dist <- Lego_matrix_maxdist(Lego = Lego, chr1 = chr1, chr2 = chr2)
    if(max(min(xaxis)-max(yaxis),(max(xaxis) - min(yaxis))) > Max.dist &
    !force){
        stop(paste("The farthest pixel loaded for ",
            "this matrix was at a distance of ",
            Max.dist,"bins from the diagonal. ",
            "The current selection subsets out-of-bounds data.\n"))
    }
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root, chr1, chr2))
    Vector <- ._Lego_Get_Something_(Group.path = Group.path, Lego = Lego, 
        Name = Reference.object$hdf.matrix.name, Start = Start, Stride = Stride,
        Count = Count, return.what = "data")
    if(is.null(FUN)){
        return(Vector)
    }else{
        return(FUN(Vector))
    }
}


#' Get the matrix metadata columns in the Lego store.
#' 
#' `Lego_get_matrix_mcols` will get the specified matrix metadata column. 
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#' 
#' @param what \strong{Required}
#' A character vector of length 1 specifying the matrix metric to retrieve
#'  
#' @return Returns a 1xN dimensional vector containing the specified matrix 
#' metric 
#' 
#' @examples 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_get_matrix_mcols(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19",
#' what = "bin.coverage")
Lego_get_matrix_mcols = function(Lego = NULL, chr1 = NULL, chr2 = NULL, 
    what = NULL){
    Reference.object <- GenomicMatrix$new()
    Meta.cols <- Reference.object$hdf.matrix.meta.cols()
    if(any(is.null(c(Lego,chr1,chr2,what)))){
        stop("Lego, chr1, chr2, what cannot be NULL.\n")
    }
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("Matrix for this chromsome pair does not exist.\n")   
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("Matrix for this chromsome pair is yet to be loaded.\n")   
    }
    if(length(what) >1){
        stop("What must be a character vector of length 1\n")      
    }
    if(!Lego_matrix_issparse(Lego = Lego, chr1 = chr1, chr2 = chr2) & what
        == Meta.cols["sparse"]){
        stop("This matrix is not a sparse matrix.",
            " So sparsity.index was not calculated\n")
    }
    if(what == Meta.cols["sparse"] & chr1 != chr2){
        stop("sparsity.index only applies to cis matrices (chr1 == chr2).\n")
    }

    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root, chr1, 
        chr2))
    Vector <- ._Lego_Get_Something_(Group.path = Group.path, Lego = Lego,
    Name = what, return.what = "data")
    return(Vector)
}


#' List the matrix metadata columns in the Lego store.
#' 
#' `Lego_get_matrix_mcols` will list the names of all matrix metadata columns. 
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @inheritParams Lego_load_matrix
#'  
#' @return Returns a vector containing the names of all matrix metadata columns 
#' 
#' @examples 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos")
#' Lego_list_matrix_mcols(Lego = Lego.file, chr1 = "chr19", chr2 = "chr19")
Lego_list_matrix_mcols = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Meta.cols <- Reference.object$hdf.matrix.meta.cols() 
    return(Meta.cols)
}


#' List all Legos that have been created with .
#' 
#' `Lego_list_tracked_legos` will list all Lego objects that are being tracked. 
#' 
#' @param detailed \strong{Optional}
#' If FALSE, produces a named vector of file paths pointing to the 
#' Lego objects. The names correspond to the original Filename provided during
#' Lego creation. If TRUE, produces a tibble containing additional information
#' such as the bioc cache tracking id, and the creation, last accession and 
#' modification times. 
#' 
#' @param preserve.col.names \strong{Optional}
#' If TRUE, will preserve the original colnames from the bioc cache tibble
#' object. If FALSE, will attempt to humanize the colnames to improve
#' readability.
#'  
#' @return Returns a vector or tibble containing the path to tracked 
#' Lego objects. If tibble will contain additional information related to 
#' the tracking. For details see parameter detailed. 
#' 
#' @examples 
#' Lego_list_tracked_legos()
Lego_list_tracked_legos = function(detailed = FALSE,
    preserve.col.names = FALSE){
    Reference.object <- GenomicMatrix$new()
    Cache.dir <- ._Get_cachedir()
    Info.tib <- bfcquery(x = Cache.dir, 
        query = paste(".",Reference.object$lego.extension,sep = ""),
        field = "rpath")
    if(nrow(Info.tib) == 0){
        return(NULL)
    }
    File.names <- Info.tib[["rname"]]
    File.paths <- Info.tib[["rpath"]]
    if(any(duplicated(File.names))){
        File.names <- paste(File.names,Info.tib[["rid"]],sep = "_")
    }
    names(File.paths) <- File.names
    if(!detailed){
        names(File.paths) <- File.names
        return(File.paths)
    }else{
        Info.tib$path_to_file <- File.paths
        Info.tib <- Info.tib[,c(1,2,5,9,3,4,8)]
        if(!preserve.col.names){
            colnames(Info.tib) <- c("Unique_db_identifier",
                "name_of_origin_file", "name_in_cache_dir", "path_to_file",
                "create_time", "access_time", "last_modified_time")
        }
        return(Info.tib)
    }
}


#' Check if a Lego is being tracked.
#' 
#' `Lego_is_tracked` checks if the provided Lego file is being tracked or not. 
#' 
#' @param Lego \strong{required}
#' Path to the Lego file to check. 
#'  
#' @return Returns TRUE or FALSE indicating if the file path is being tracked
#' or not.
#' 
#' @examples
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos") 
#' Lego_is_tracked(Lego = Lego.file)
Lego_is_tracked = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Cache.dir <- ._Get_cachedir()
    Info.tib <- bfcquery(x = Cache.dir, query = Lego, field = "rname")
    nrow(Info.tib) > 0 
}


#' Check if a Lego is being tracked.
#' 
#' `Lego_track_legos` will start tracking the provided Lego file. 
#' 
#' @param Lego \strong{required}
#' Path to the Lego file to track. 
#'  
#' @return Returns a named vector containing the path to the file along with 
#' the bioc Cache id as its name.
#' 
#' @examples 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos") 
#' Lego_track_legos(Lego = Lego.file)
Lego_track_legos = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Cache.dir <- ._Get_cachedir()
    Dir.path <- dirname(Lego)
    Filename <- basename(Lego)
    Lego <- file.path(normalizePath(Dir.path),Filename)
    if(Lego_is_tracked(Lego = Lego)){
        stop("This Lego object is already being tracked!")
    }
    return(bfcadd(x = Cache.dir, rname = Lego, 
            ext = Reference.object$lego.extension,
        rtype = "local", action = "asis"))
}


#' Check if a Lego is being tracked.
#' 
#' `Lego_untrack_lego` will stop tracking the provided Lego file, if it exists. 
#' 
#' @param Lego \strong{required}
#' Path to the Lego file to untrack. 
#' 
#' @examples 
#' Lego.file <- system.file("extdata", "test.hdf", package = "HiCLegos") 
#' Lego_untrack_lego(Lego = Lego.file)
Lego_untrack_lego = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Cache.dir <- ._Get_cachedir()
    Dir.path <- dirname(Lego)
    Filename <- basename(Lego)
    Lego <- file.path(normalizePath(Dir.path),Filename)
    Info.tib <- bfcquery(x = Cache.dir, query = Lego, field = "rname")
    if(Lego_is_tracked(Lego = Lego)){
       DB_id <- Info.tib$rid
       bfcremove(Cache.dir, DB_id)
    }else{
        warning("Lego is not being tracked!")
    }
}
