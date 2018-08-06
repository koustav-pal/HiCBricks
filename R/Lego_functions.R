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
#' @param ChromNames **Required**. 
#' A character vector containing the chromosomes to be considered for the 
#' dataset. This string is used to verify the presence of all chromosomes in the
#' provided bitable. 
#' 
#' @param BinTable **Required**. 
#' A string containing the path to the file to load as the binning table for the
#' Hi-C experiment. The number of entries per chromosome defines the dimension 
#' of the associated Hi-C data matrices. For example, if chr1 contains 250 
#' entries in the binning table, the _cis_ Hi-C data matrix for chr1 will be
#' expected to contain 250 rows and 250 cols. Similary, if the same binning
#' table contained 150 entries for chr2, the _trans_ Hi-C matrices for chr1,chr2
#' will be a matrix with dimension 250 rows and 150 cols.
#' 
#' There are no constraints on the bintable format. As long as the table is in a
#' delimited format, the corresponding table columns can be outlined with the
#' associated parameters. The columns of importance are chr, start and end.
#' 
#' It is recommended to always use binning tables where the end and start of 
#' consecutive ranges are not the same. If they are the same, this may lead to
#' **unexpected behaviour** when using the GenomicRanges "any" overlap function.
#'
#' @param Output.Filename **Required**.
#' A string specifying the location and name of the HDF file to create. If path
#' is not provided, it will be created in the current working directory.
#' 
#' @param bin.delim **Optional**. Default "\t".
#' A character vector of length 1 specifying the delimiter used in the file 
#' containing the binning table.
#' 
#' @param col.index **Optional**. Default "c(1,2,3)".
#' A character vector of length 3 containing the indexes of the required columns
#' in the binning table. the first index, corresponds to the chr column, the
#' second to the start column and the third to the end column.
#' 
#' @param impose.discontinuity **Optional**. Default TRUE.
#' If TRUE, this parameter ensures a check to make sure that required the end
#' and start coordinates of consecutive entries are not the same per chromosome.
#' 
#' @param ChunkSize **Optional**.
#' A numeric vector of length 1. If provided, the HDF dataset will use this 
#' value as the chunk size, for all matrices. By default, the ChunkSize is 
#' set to matrix dimensions/100. 
#' 
#' @param exec **Optional**. Default cat.
#' A string specifying the program or expression to use for reading the file.
#' For bz2 files, use bzcat and for gunzipped files use zcat.
#' 
#' @param remove.existing **Optional**. Default FALSE.
#' If TRUE, will remove the HDF file with the same name and create a new one.
#' By default, it will not replace existing files.
#' 
#' @details The structure of the HDF file is as follows: 
#' The structure contains three major groups which are then hierarchically
#' nested with other groups to finally lead to the corresponding datasets.  
#' _Base.matrices_ **group**
#'      _chromosome_ **group**
#'          _chromosome_ **group**, **Attributes:** Filename, Min, Max, Done
#'              matrix **dataset**
#'              bin.coverage **dataset**
#'              row.sums **dataset**
#'              sparsity **dataset**
#' _Base.ranges_ **group**
#'      _Bintable_ **group**
#'          ranges **dataset**
#'          offsets **dataset**
#'          lengths **dataset**
#'          chr.names **dataset**
#'      _YourRangesTable_ **group**
#'          ranges **dataset**
#'          offsets **dataset**
#'          lengths **dataset**
#'          chr.names **dataset**
#' _Base.metadata_ **group**
#'      chromosomes **dataset**
#'      YourMetadataTable **dataset**
CreateLego <- function(ChromNames=NULL, BinTable=NULL, bin.delim="\t",
    col.index=c(1,2,3), impose.discontinuity=TRUE, ChunkSize=NULL, 
    Output.Filename=NULL, exec="cat", remove.existing=FALSE){
    require(rhdf5)
    H5close()
    Dir.path <- dirname(Output.Filename)
    Filename <- basename(Output.Filename)
    Working.File <- file.path(normalizePath(Dir.path),Filename)
	Reference.object <- GenomicMatrix$new()
    Root.folders <- Reference.object$GetRootFolders()
    if(is.null(ChromNames) | length(ChromNames) == 0){
        stop("Variable ChromNames cannot be empty")   
    }
    HDF.File <- file.path(Working.File)
    if(file.exists(HDF.File)){
        if(remove.existing){
            file.remove(HDF.File)
        }else{
            stop("Provided HDF file already exists. Please provide remove.existing = TRUE to overwrite it\n")
        }
    }
    ChromosomeList<-ChromNames
    if(is.null(BinTable)){
        stop("Variable Bintable cannot be empty. Binning information must be provided at startup")
    }
    # Read in the binning table
    cat("Reading Bintable:",BinTable,"\n")
    Bintable.list <- Read_bintable(Filename = BinTable, read.delim = bin.delim, exec = exec,
                col.index = col.index, chromosomes = ChromosomeList,
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
    Chrom.lengths <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = length, col.name = 'chr')
    Chrom.sizes <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = max, col.name = 'end')
    Chrom.info.df <- data.frame(chr = names(Chrom.lengths),
        nrow = as.vector(Chrom.lengths),
        size = as.vector(Chrom.sizes),stringsAsFactors = FALSE)
    # Create metadata chromosome groups
    ._Lego_WriteDataFrame_(Lego = HDF.File, Path = c(Root.folders['metadata']), name = Reference.object$metadata.chrom.dataset, object = Chrom.info.df)
    ._Lego_Add_Ranges_(Group.path = Create_Path(c(Root.folders['ranges'],Reference.object$hdf.bintable.ranges.group)), Lego = HDF.File, 
        ranges.df = Bintable, name = Reference.object$hdf.ranges.dataset.name, mcol.list = NULL)
    # Create matrices groups
    for (chrom1 in ChromosomeList) {
        CreateGroups(Group.path = Create_Path(c(Root.folders['matrices'],chrom1)), File = HDF.File)
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
            Dims <- c(Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"], Chrom.info.df[Chrom.info.df$chr == chrom2,"nrow"])
            if(is.null(ChunkSize)){
                ChunkSize <- ceiling(Dims/100)
            }
            Array.dim <-Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"]
            CreateDataset(Path = c(Root.folders['matrices'],chrom1,chrom2), File = HDF.File, 
                name = Reference.object$hdf.matrix.name, dims = Dims, maxdims = Dims)
        }
    }
}

#' Get the chrominfo for the Hi-C experiment.
#' 
#' `Lego_get_chrominfo` fetches the associated chrominfo table for the 
#' Lego (HDF) it is associated to. 
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @return A three column file containing chromosomes, nrows and length. 
#' chromosomes corresponds to all chromosomes in the provided bintable, nrows
#' corresponds to the number of entries in the bintable or dimension for that 
#' chromosome in the Hi-C matrix. And length is the total bp length of the same
#' chromosome (max value for that chromosome in the bintable).
Lego_get_chrominfo <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Dataset <- ._Lego_Get_Something_(Group.path = Reference.object$hdf.metadata.root, 
        Lego = Lego, Name = Reference.object$metadata.chrom.dataset, return.what = "data")
    return(Dataset)
}


#' Creates a ranges object from provided vectors.
#' 
#' `Lego_make_ranges` creates a GRanges object from the provided arguments
#'  
#' @param Chrom **Required**.
#' A 1 dimensional character vector of size N specifying the chromosomes in the
#' ranges.
#' 
#' @param Start **Required**.
#' A 1 dimensional numeric vector of size N specifying the start positions in 
#' the ranges.
#' 
#' @param End **Required**.
#' A 1 dimensional numeric vector of size N specifying the end positions in 
#' the ranges. Must be less than Start.
#' 
#' @param Strand **Optional**.
#' A 1 dimensional character vector of size N specifying the strand of the 
#' ranges. If not provided, this will be set to the default *.
#' 
#' @param Names **Optional**.
#' A 1 dimensional character vector of size N specifying the names of the 
#' ranges. If not provided, this will be set to the default chr:start:end.
Lego_make_ranges = function(Chrom=NULL, Start=NULL, End=NULL, Strand=NULL, 
    Names=NULL){
    Reference.object <- GenomicMatrix$new()
    require(GenomicRanges)
    if(is.null(Names)){
        Names<-paste(Chrom,as.integer(Start),as.integer(End),sep=Reference.object$Ranges.separator)
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

#' Perpetually store a ranges object in the Lego store.
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
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param ranges **Required**.
#' An object of class ranges specifying the ranges to store in the Lego.
#' 
#' @param rangekey **Required**.
#' The name to use for the object within the Lego store.
Lego_add_ranges = function(Lego = NULL, ranges = NULL, rangekey = NULL){
    Reference.object <- GenomicMatrix$new()
    if(!(class(ranges) %in% "GRanges") | ("list" %in% class(ranges))){
        stop("Object of class Ranges expected")
    }
    Ranges.df <- as.data.frame(ranges)
    if(is.unsorted(Ranges.df$seqnames)){
        stop("Ranges must be sorted by chromosome!")
    }
    if(Lego_rangekey_exists(Lego = Lego, rangekey = rangekey)){
        stop("rangekey already exists! Cannot proceed further! Please read the documentation to understand Why.")
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
    ._Lego_Add_Ranges_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root,rangekey)), Lego = Lego, 
        ranges.df = Ranges.df.coords, name = rangekey, mcol.list = Metadata.list)
}

#' List the ranges tables stored within the Lego.
#' 
#' `Lego_list_rangekeys` lists the names of all ranges associated to a Lego.
#'  
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @return A one dimensional character vector of length x 
Lego_list_rangekeys = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Handler <- ._Lego_Get_Something_(Group.path = Create_Path(Reference.object$hdf.ranges.root), 
        Lego = Lego, return.what = "group_handle")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    return(GroupList)
}

#' Check to see if the Lego contains a ranges with a certain name.
#' 
#' `Lego_rangekey_exists` checks for the presence of a particular ranges with
#' a certain name.
#'  
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param rangekey **Required**.
#' A string specifying the name of the ranges to check for.
#' 
#' @return A logical vector of length 1 with either TRUE or FALSE values. 
Lego_rangekey_exists = function(Lego = NULL, rangekey = NULL){
    Keys <- Lego_list_rangekeys(Lego = Lego)
    return(rangekey %in% Keys)
}

#' Find out what metadata columns are associated to a ranges with a certain name
#' 
#' `Lego_list_ranges_mcols` will list the metadata columns of the specified
#' ranges if it is present in the Lego store. 
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param rangekey **Optional**.
#' A string specifying the name of the ranges. If not present, the metadata
#' columns of all ranges will be listed. 
#' 
#' @return if no metadata columns are present, NA. If metadata columns are
#' present, a data.frame object containing the name of the ranges and the 
#' associated metadata column name.
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
        Handler <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, x)), Lego = Lego, return.what = "group_handle")
        GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
        H5Gclose(Handler)
        data.frame(rangekey = x, m.col = GroupList)
    })
    mcol.df <- do.call(rbind,mcol.list)
    mcol.df <- mcol.df[!(mcol.df$m.col %in% Reference.object$hdf.ranges.protected.names()),]
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
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @return Returns a data.frame object with columns chr1, chr2 corresponding 
#' to chromosome pairs, and the associated attributes. filename corresponds to
#' the name of the file that was loaded for the pair. min and max specify the
#' minimum and maximum values in the matrix, done is a logical value 
#' specifying if a matrix has been loaded and sparsity specifies if a matrix
#' is defined as a sparse matrix.
Lego_list_matrices = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    chr1.list <- lapply(ChromInfo[,"chr"], function(chr1){
        chr2.list <- lapply(ChromInfo[,"chr"], function(chr2){
            Colnames <- Reference.object$matrices.chrom.attributes
            Values <- GetAttributes(Path = Create_Path(c(Reference.object$hdf.matrices.root, chr1, chr2)),
                File = Lego, Attributes = Colnames, on = "group")
            temp.df <- data.frame(chr1,chr2,t(cbind(Values)))
            colnames(temp.df) <- c("chr1","chr2",Colnames)
            temp.df
        })
        chr2.df <- do.call(rbind,chr2.list)
    })
    Matrix.list.df <- do.call(rbind,chr1.list)
    Matrix.list.df[,"done"] <- as.logical(as.integer(as.character(Matrix.list.df[,"done"])))
    Matrix.list.df[,"sparsity"] <- as.numeric(as.character(Matrix.list.df[,"sparsity"]))
    Matrix.list.df[,"min"] <- as.numeric(as.character(Matrix.list.df[,"min"]))
    Matrix.list.df[,"max"] <- as.numeric(as.character(Matrix.list.df[,"max"]))
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
#' other than "Bintable", and had been added using Lego_add_ranges the width and
#' Strand columns may appear as metadata columns. These will most likely be 
#' artifacts from converting the original ranges object to a data.frame.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr **Optional**.
#' A chr string specifying the chromosome to select from the ranges.
#'
#' @param rangekey **Required**.
#' A string specifying the name of the ranges.
#' 
#' @return Returns a GRanges object with the associated metadata columns that
#' may have been present in the Ranges object. 
Lego_get_ranges = function(Lego = NULL, chr = NULL, rangekey = NULL){
    Reference.object <- GenomicMatrix$new()
    if(is.null(rangekey) | is.null(Lego)){
        stop("rangekey and Lego cannot remain empty!\n")
    }
    if(!Lego_rangekey_exists(Lego = Lego, rangekey = rangekey)){
        stop("rangekey not found!")
    }
    Start <- NULL
    Stride <- NULL
    Count <- NULL
    Index <- NULL
    if(!is.null(chr)){
        chromosomes <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.chr.name, return.what = "data")
        if(!(chr %in% chromosomes)){
            stop("chr not found!")
        }        
        Starts <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.offset.name, return.what = "data")
        Lengths <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.lengths.name, return.what = "data")
        Which.one <- chromosomes == chr
        Start <- Starts[Which.one]
        Stride <- 1
        Count <- Lengths[Which.one]
    }
    Dataset <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.dataset.name, Start = Start, Stride = Stride,
        Count = Count, return.what = "data")
    Dataset <- Lego_make_ranges(Chrom = Dataset[,'chr'], Start = Dataset[,'start'], End = Dataset[,'end'])

    MCols <- Lego_list_ranges_mcols(Lego = Lego, rangekey = rangekey)
    if(class(MCols) == "data.frame"){
        MCols.col <- as.character(MCols[,"m.col"])
        m.start <- ifelse(is.null(Start),NULL,Start[1])
        m.stride <- ifelse(is.null(Start),NULL,Stride[1])
        m.count <- ifelse(is.null(Start),NULL,Count[1])
        MCols.DF.list <- lapply(MCols.col,function(x){
            Dataset <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
                Lego = Lego, Name = x, Start = m.start, Stride = m.stride,
                Count = m.count, return.what = "data")
            DF <- DataFrame(Temp = Dataset)
            colnames(DF) <- x
            DF
        })
        MCols.DF <- do.call(cbind,MCols.DF.list)
        mcols(Dataset) <- MCols.DF
    }
    return(Dataset)
}


#' Returns the binning table associated to the Hi-C experiment.
#' 
#' `Lego_get_bintable` makes a call to Lego_get_ranges to retrieve the 
#' binning table of the associated Lego store. This is equivalent to passing
#' the argument rangekey = "bintable" in Lego_get_ranges
#' 
#' @seealso Lego_get_ranges
Lego_get_bintable = function(Lego = NULL, chr = NULL){
    Reference.object <- GenomicMatrix$new()
    Table <- Lego_get_ranges(Lego = Lego, chr = chr, rangekey = Reference.object$hdf.bintable.ranges.group)
    return(Table)
}

#' Returns the position of the supplied ranges in the binning table associated 
#' to the Hi-C experiment.
#' 
#' `Lego_fetch_range_index` constructs a ranges object and creates an overlap
#' operation between the constructed ranges and the Hi-C experiment associated
#' binning table finally returning a list of ranges with their corresponding 
#' indices in the binning table.
#'  
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr **Required**.
#' A character vector of length N specifying the chromosomes to select from 
#' the ranges.
#' 
#' @param start **Required**.
#' A numeric vector of length N specifying the start positions in the chromosome
#' 
#' @param end **Required**.
#' A numeric vector of length N specifying the end positions in the chromosome
#' 
#' @param names **Optional**.
#' A character vector of length N specifying the names of the chromosomes. If
#' absent, names will take the form chr:start:end.
#' 
#' @param type **Optional**. Default any
#' Type of overlap operation to do. It should be one of two, any or within.
#' any considers any overlap (atleast 1 bp) between the provided ranges and the 
#' binning table.
#' 
#' @return Returns a named list corresponding to one entry for each unique
#' chromosome present in the chr character vector. The list contains another
#' list with length equal to the number of ranges present per chromosome.
#' This list contains two named vectors, SubjectInfo and Indexes. 
#' SubjectInfo contains a ranges created out of the vector from any given 
#' position. Indexes contains a numeric vector corresponding to the rows/cols 
#' overlapping the given ranges.  
Lego_fetch_range_index = function(Lego = NULL, chr = NULL, start = NULL, end = NULL,names = NULL,type = "any"){
    AllTypes<-c("any","within")
    if( any(!(type %in% AllTypes)) ){
        stop("type takes one of two arguments: c(\"any\",\"within\")")
    }
    if(is.null(chr) | is.null(start) | is.null(end) | is.null(Lego)){
        stop("Chrom, start, end and Lego cannot be empty")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    if(any(!(chr %in% ChromInfo[,'chr']))){
        stop("Provided chr does not exist in chromosome list.")
    }
    if(any(class(chr) != "character" | !(class(start) %in% c("numeric","integer")) | !(class(end) %in% c("numeric","integer")))){
        stop("Provided chr, start, end do not match expected class definitions of character, numeric, numeric")
    }
    Unique.chromosomes <- unique(chr)
    OverlapByChromosome <- lapply(Unique.chromosomes,function(cur.chr){
        Filter <- chr==cur.chr
        Cur.Chrom <- chr[Filter]
        Cur.Start <- start[Filter]
        Cur.end <- end[Filter]
        Cur.Names <- names[Filter]
        SubjectRanges <- Lego_get_bintable(Lego = Lego, chr = chr)
        if( any(!(Cur.end <= max(end(SubjectRanges)) & Cur.Start >= min(start(SubjectRanges)))) ){
            stop("Start or end is out of ranges for Bintable")
        }
        QueryRanges <- Lego_make_ranges(Chrom=Cur.Chrom, Start=Cur.Start, End=Cur.end, Names=Cur.Names)
        require(GenomicRanges)
        HitsObject <- findOverlaps(SubjectRanges,QueryRanges,type=type)
        UniqueQueries <- unique(subjectHits(HitsObject))
        MatchingIndexes <- lapply(1:length(UniqueQueries),function(x){
            A.Query <- UniqueQueries[x]
            MatchingQueries <- queryHits(HitsObject)[subjectHits(HitsObject)==A.Query]
            ListObj<-list()
            ListObj[["Indexes"]]<-MatchingQueries
            ListObj[["SubjectInfo"]]<-QueryRanges[A.Query]
            ListObj
            })
    })
    names(OverlapByChromosome)<-Unique.chromosomes
    return(OverlapByChromosome)
}

#' Provides the corresponding overlapping position from the bintable.
#' 
#' `Lego_return_region_position` takes as input a human-readable coordinate
#' format of the form chr:start:end and outputs the overlapping bintable 
#' positions. 
#'  
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param region **Required**.
#' A character vector of length 1 specifying the region to overlap. It must take
#' the form chr:start:end.
#' 
#' @return Returns a 1 dimensional vector containing the position of the
#' overlapping regions in the bintable associated the Lego store. 
Lego_return_region_position = function(Lego = NULL, region=NULL){
    if(!is.character(region) | length(region) > 1){
        stop("region must be a character vector of length 1")
    }
    Coord.Split<- Split_genomic_coordinates(Coordinate=region)
    region.chr<-Coord.Split[[1]][1]
    region.start<-as.numeric(Coord.Split[[1]][2])
    region.stop<-as.numeric(Coord.Split[[1]][3])
    region.ranges<-Lego_fetch_range_index(Lego = Lego, chr=region.chr, start=region.start, end=region.stop, type="within")
    Vector.coordinates <- Region.Ranges[[chr]][[1]][["Indexes"]]
    return(Vector.coordinates)
}

#' Load a NxM dimensional matrix into the Lego store.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' rows of the matrix
#' 
#' @param chr2 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' columns of the matrix
#' 
#' @param matrix.file **Required**.
#' A character vector of length 1 specifying the name of the file to load as a
#' matrix into the Lego store. 
#' 
#' @param exec **Required**.
#' A string specifying the program to use for reading the file. Use cat for txt
#' files, for bz2 files use bzcat and for gz files zcat.
#' 
#' @param delim **Optional**. Default " "
#' The delimiter of the matrix file.
#' 
#' @param remove.priori **Optional**. Default FALSE
#' If a matrix was loaded before, it will not be replaced. Use remove.priori to
#' override and replace the existing matrix.
#' 
#' @param num.rows **Optional**. Default 2000
#' Number of rows to read, in each chunk.
#' 
#' @param distance **Optional**. Default NULL. Not implemented yet.
#' For very high-resolution matrices, read times can become extremely slow and
#' it does not make sense to load the entire matrix into the data structure, as
#' after a certain distance, the matrix will become extremely sparse. If
#' provided, only interactions upto a certain distance from the main diagonal
#' will be loaded into the data structure.
#' 
#' @param is.sparse **Optional**. Default FALSE
#' If true, designates the matrix as being a sparse matrix, and computes the
#' sparsity.index. The sparsity index measures the proportion of non-zero rows
#' or columns at a certain distance from the diagonal (100) in cis interaction
#' matrices.
#'
#' @param sparsity.bins **Optional**. Default 100
#' With regards to computing the sparsity.index, this parameter decides the 
#' number of bins to scan from the diagonal. 
Lego_load_matrix = function(Lego = NULL, chr1 = NULL, chr2 = NULL, matrix.file = NULL, delim = " ", exec = NULL,  
    remove.prior = FALSE, num.rows = 2000, is.sparse = FALSE, distance = NULL, sparsity.bins = 100){
    Reference.object <- GenomicMatrix$new()
    ListVars <- list(Lego = Lego, chr1 = chr1, chr2 = chr2, file = file, is.sparse = is.sparse, 
        sparsity.bins = sparsity.bins, exec = exec, delim = delim, distance = distance, remove.prior = remove.prior)
    sapply(1:length(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    sapply(1:length(ListVars[c("Lego","chr1","chr2","file","exec")]),function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    if(Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2) && !remove.prior){
        stop("A matrix was preloaded before. Use remove.priori = TRUE to force value replacement\n")
    }
    Chrom.info.df <- Lego_get_chrominfo(Lego = Lego)
    chr1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr1,"nrow"]
    chr2.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr2,"nrow"]
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,chr1,chr2))
    compute.sparsity <- FALSE
    if(is.sparse && chr1 == chr2){
        compute.sparsity <- TRUE
    }
    ._ProcessMatrix_(Lego = Lego, Matrix.file = matrix.file, delim = delim, exec = exec, Group.path = Group.path, 
        chr1.len = chr1.len, chr2.len = chr2.len, num.rows = num.rows, is.sparse = is.sparse, 
        compute.sparsity = compute.sparsity, distance = distance, sparsity.bins = sparsity.bins)
}

#' Check if a matrix has been loaded for a chromosome pair.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' rows of the matrix
#' 
#' @param chr2 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' columns of the matrix
#' 
#' @return Returns a logical vector of length 1, specifying if a matrix has
#' been loaded or not.   
Lego_matrix_isdone = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_list_matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    return(as.logical((Matrix.list[Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2, "done"])))
}

#' Check if a chromosome pair exists. This helps if the user thinks that they
#' may have made a mistake while creating the initial Lego file.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' rows of the matrix
#' 
#' @param chr2 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' columns of the matrix
#' 
#' @return Returns a logical vector of length 1, specifying if the matrix exists
#' or not.
Lego_matrix_exists = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    ChromInfo.df <- Lego_get_chrominfo(Lego = Lego)
    all(c(chr1,chr2) %in% ChromInfo.df[,"chr"])
}

#' Return the value range of the matrix
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' rows of the matrix
#' 
#' @param chr2 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' columns of the matrix
#' 
#' @return Returns a numeric vector of length 2, specifying the minimum and
#' maximum finite real values in the matrix.
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

#' Return the filename of the loaded matrix
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' rows of the matrix
#' 
#' @param chr2 **Required**.
#' A character vector of length 1 specifying the chromosome corresponding to the
#' columns of the matrix
#' 
#' @return Returns a character vector of length 1 specifying the initial 
#' filename of the matrix.
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
#' `Lego_get_values_by_distance` can fetch values with or without transformation
#' or subsetted by a certain distance.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr **Required**.
#' A string specifying the chromosome for the cis Hi-C matrix from which values
#' will be retrieved at a certain distance.
#' 
#' @param distance **Required**. 0 based.
#' Fetch values separated by distance.
#'
#' @param constrain.region **Optional**.
#' A character vector of length 1 with the form chr:start:end specifying the 
#' region for which the distance values must be retrieved.
#' 
#' @param batch.size **Optional**. Default 500
#' A numeric vector of length 1 specifying the size of the chunk to retrieve
#' for diagonal selection.
#' 
#' @param FUN **Optional**.
#' If provided a data transformation with FUN will be applied before values
#' are returned.
#' 
#' @return Returns a numeric vector of length N depending on the presence of 
#' constrain.region, FUN and distance from the main diagonal.
Lego_get_values_by_distance = function(Lego = NULL, chr = NULL, distance  = NULL,
    constrain.region=NULL,batch.size=500,FUN=NULL){
    if(any(sapply(list(Lego,chr,distance),is.null))) {
        stop("Lego, chr, distance cannot be NULL.\n")
    }
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr, chr2 = chr) & !Lego_matrix_isdone(Lego = Lego, chr1 = chr, chr2 = chr)){
        stop("Chromosome is not listed or has not been loaded in this HDF file.\n")
    }
    if(any(sapply(list(Lego,chr,distance),length) > 1)){
        stop("Lego, chr and distance can only have values of length 1.\n")
    }
    Nrow <- (ChromInfo[ChromInfo$chr == chr,"nrow"]-1)
    if(any(distance > Nrow, distance < 0)){
        stop("distance must range between 0 and",Nrow,"\n")
    }
    Root.folders <- Reference.object$GetRootFolders()
    Path <- Create_Path(c(Root.folders['matrices'],chr,chr))
    Vector.start <- 1
    Vector.stop <- ChromInfo[ChromInfo$chr==chr,"nrow"]
    if(!is.null(constrain.region)){
        Vector.coordinates <- Lego_return_region_position(Lego = Lego, region=constrain.region, chr=chr)
        if(is.null(Vector.coordinates)){
            stop("Overlap operation was unsuccessful! Please check coordinates ",Constrain.region)
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
        CumSums <- cumsum(c(0,Counts[1:(length(Counts)-1)]))
    }
    DistancesVector.list <- lapply(1:length(Counts),function(x){
        Count <- Counts[x]
        Offset <- CumSums[x]
        cur.start <- Start+Offset
        diag(._Lego_Get_Something_(Group.path = Path, Lego = Lego, Name = Reference.object$hdf.matrix.name,
            Start= cur.start, Stride = Stride, Count = c(Count,Count), return.what = "data"))
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
#' `Lego_get_matrix_within_coords` will fetch a matrix subset after creating an
#' overlap operation between both regions and the bintable associated to the
#' Lego store. This function calls Lego_get_matrix.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param x.coords **Required**.
#' A string specifying the region to subset on the rows. It takes the form
#' chr:start:end. An overlap operation with the associated bintable will be done
#' to identify the bins to subset on the row
#' 
#' @param y.coords **Required**.
#' A string specifying the region to subset on the rows. It takes the form
#' chr:start:end. An overlap operation with the associated bintable will be done
#' to identify the bins to subset on the column
#' 
#' @param constrain.region **Optional**.
#' A character vector of length 1 with the form chr:start:end specifying the 
#' region for which the distance values must be retrieved.
#' 
#' @param FUN **Optional**.
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#' 
#' @return Returns a matrix of dimension x.coords binned length by y.coords
#' binned length. This may differ based on FUN.
#' 
#' @seealso Lego_get_matrix
Lego_get_matrix_within_coords = function(Lego = NULL, x.coords=NULL, y.coords=NULL, FUN=NULL){
    type <- "within"
    if( (is.null(x.coords)) | (is.null(y.coords)) ){
        stop("x.coords, y.coords and Lego cannot be NULL")
    }
    if(!(length(x.coords)==1) | !(length(y.coords)==1)){
        stop("This function processes single process calls at a time. ",
            "Setup an Iterator for more functionality")
    }
    if( class(x.coords)!="character" | class(y.coords)!="character" ){
        stop("Two string variables were expected for x.coords & y.coords,\nfound XCoords class ", class(x.coords), " and YCoords class ", class(y.coords))
    }
    XCoord.split <- Split_genomic_coordinates(Coordinate=x.coords)
    XCoord.Chrom <- XCoord.split[[1]][1]
    XCoord.start <- as.numeric(XCoord.split[[1]][2])
    XCoord.stop <- as.numeric(XCoord.split[[1]][3])
    YCoord.split <- Split_genomic_coordinates(Coordinate=y.coords)
    YCoord.Chrom <- YCoord.split[[1]][1]
    YCoord.start <- as.numeric(YCoord.split[[1]][2])
    YCoord.stop <- as.numeric(YCoord.split[[1]][3])
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    if( any(!(c(XCoord.Chrom,YCoord.Chrom) %in% ChromosomeList)) ){
        stop("Provided chromosomes were not found in chromosome list.")
    }
    chr1 = XCoord.Chrom
    chr2 = YCoord.Chrom
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop(chr1,chr2," matrix is yet to be loaded into the class.")
    }
    XCoords.list <- Lego_fetch_range_index(Lego = Lego, chr=chr1, start=XCoord.start, end=XCoord.stop, names=NULL, type=type)
    YCoords.list <- Lego_fetch_range_index(Lego = Lego, chr=chr2, start=YCoord.start, end=YCoord.stop, names=NULL, type=type)
    if( is.null(XCoords.list[[chr1]][[1]][["Indexes"]]) | is.null(YCoords.list[[chr2]][[1]][["Indexes"]]) ){
        stop("Overlap operation was unsuccessful! Please check coordinates ",x.coords," & ",y.coords)
    }
    x.vector <- XCoords.list[[chr1]][[1]][["Indexes"]]
    y.vector <- YCoords.list[[chr2]][[1]][["Indexes"]]
    Matrix <- Lego_get_matrix(Lego = Lego, chr1=chr1, chr2=chr2, x.vector=x.vector, y.vector=y.vector, FUN=FUN)
    return(Matrix)
}

#' Return a matrix subset.
#' 
#' `Lego_get_matrix` will fetch a matrix subset between row values ranging from 
#' min(x.vector) to max(x.vector) and column values ranging from min(x.vector) 
#' to max(x.vector)
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A string specifying the chromosome to subset on the rows
#' 
#' @param chr2 **Required**.
#' A string specifying the chromosome to subset on the columns
#' 
#' @param x.vector **Required**.
#' A one-dimensional numeric vector specifying the the rows to subset. 
#' 
#' @param y.vector **Required**.
#' A one-dimensional numeric vector specifying the the columns to subset. 
#' 
#' @param FUN **Optional**.
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#' 
#' @return Returns a matrix of dimension x.vector length by y.vector length. 
#' This may differ based on the operations with FUN.
Lego_get_matrix = function(Lego = NULL, chr1=NULL, chr2=NULL, x.vector=NULL, y.vector=NULL, FUN=NULL){
    # cat(" Rows: ",x.vector," Cols: ",y.vector,"\n")
    if(any(!(class(x.vector) %in% c("numeric","integer")) | !(class(y.vector) %in% c("numeric","integer")))){
        stop("x.vector and y.vector must be numeric.\n")
    }
    if(is.null(chr1) | is.null(chr2) | is.null(x.vector) | is.null(y.vector)){
        stop("Either of chr1, chr2, x.vector or y.vector were provided as NULL values.\n")
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
    if(any(x.vector > chr1.len) | any(y.vector > chr2.len) | min(x.vector,y.vector) < 1 ) {
        stop("x.vector or y.vector falls outside the bounds of loaded Bintables") 
    }
    Matrix <- Lego_get_vector_values(Lego = Lego, chr1=chr1, chr2=chr2, xaxis=x.vector, yaxis=y.vector)
    if(is.null(FUN)){
        return(Matrix)              
    }else{
        return(FUN(Matrix))
    }
}

#' Return a matrix subset.
#' 
#' `Lego_fetch_row_vector` will fetch any given rows from a matrix. If required,
#' the rows can be subsetted on the columns and transformations applied.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A string specifying the chromosome to subset on the rows
#' 
#' @param chr2 **Required**.
#' A string specifying the chromosome to subset on the columns
#' 
#' @param by **Required**. One of two possible values, "position" or "ranges"
#' A one-dimensional numeric vector of length 1 specifying one of either 
#' position or ranges. 
#' 
#' @param vector **Required**. 
#' If by is position, a 1 dimensional numeric vector containing the rows to be
#' extracted is expected. If by is ranges, a 1 dimensional numeric vector 
#' containing the names of the bintable is expected. 
#' This function does not do overlaps. Rather it returns any given row or column
#' based on their position or names in the bintable.
#' 
#' @param regions **Optional**.
#' A character vector of length vector is expected. Each element must be of the
#' form chr:start:end. These regions will be converted back to their original
#' positions and the corresponding rows will be subsetted by the corresponding 
#' region element. If the length of regions does not match, the subset operation
#' will not be done and all elements from the rows will be returned.
#'
#' @param flip **Optional**. Default FALSE
#' If present, will flip everything. This is equivalent to selecting columns,
#' and subsetting on the rows.
#' 
#' @param FUN **Optional**.
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#' 
#' @return Returns a matrix of dimension x.vector length by y.vector length. 
#' This may differ based on the operations with FUN.
Lego_fetch_row_vector = function(Lego = NULL, chr1=NULL, chr2=NULL, by=NULL, vector=NULL, regions=NULL, flip = FALSE, FUN=NULL){
    Chrom.all <- c(chr1,chr2)
    if(is.null(chr1) | is.null(chr2) | is.null(by) | is.null(vector)) {
        stop("Either of chr1, chr2, by or vector were provided as NULL values")
    }
    if(!(by %in% c("position","ranges")) | length(by) != 1){
        stop("by expects a vector of type character, length 1 and takes either one of position or ranges as values")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    if(class(chr1)!="character" | class(chr2)!="character"){
        stop("Provided Chromosomes does not appear to be of class character")
    }
    if(!Lego_matrix_isdone(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop(chr1,chr2," matrix is yet to be loaded.")
    }
    if(by=="position"){
        Class.type <- c("numeric","integer")
    }
    if(by=="ranges"){
        Class.type <- "character"
    }
    if(!(class(vector) %in% Class.type)){
        stop("vector must be of class",
            ifelse(length(Class.type)>1,paste(Class.type,collapse=" or "),paste(Class.type)),"when by has value",by)
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
            stop("Position vector falls outside the bounds of ",chr1," Bintable") 
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
            y <- Lego_return_region_position(region=region,Chrom=chr2)
        }else{
            y <- 1:length(chr2.ranges)
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
        Values <- Lego_get_vector_values(Lego = Lego, chr1=chr1, chr2=chr2, xaxis=x, yaxis=y, FUN=FUN)
        return(Values)
    }) 
    return(Vector.values)
}

#' Return a N dimensional vector selection.
#' 
#' `Lego_get_vector_values` is the base function being used by all other matrix
#' retrieval functions.
#' 
#' @param Lego **Required**.
#' A string specifying the path to the Lego store created with CreateLego.
#' 
#' @param chr1 **Required**.
#' A string specifying the chromosome to subset on the rows
#' 
#' @param chr2 **Required**.
#' A string specifying the chromosome to subset on the columns
#' 
#' @param xaxis **Required**. 
#' A 1 dimensional vector containing the rows to retrieve. Gaps in this vector 
#' may result in unexpected behaviour as the values which are considered are 
#' min(xaxis) and max(xaxis) for retrieval.
#' 
#' @param yaxis **Required**.
#' A 1 dimensional vector containing the columns to retrieve. Gaps in this 
#' vector may result in unexpected behaviour as the values which are considered
#' are min(yaxis) and max(yaxis) for retrieval.
#' 
#' @param FUN **Optional**.
#' If provided a data transformation with FUN will be applied before the matrix
#' is returned.
#' 
#' @return Returns a vector of length yaxis if length of xaxis is 1. Else
#' returns a matrix of dimension xaxis length by yaxis length.
#'  
#' @section Note: Whatever the length of xaxis or yaxis may be, the coordinates 
#' under consideration will range from min(xaxis) to max(xaxis) on the rows or
#' min(yaxis) to max(yaxis) on the columns.
Lego_get_vector_values = function(Lego = NULL, chr1=NULL, chr2=NULL, xaxis=NULL, yaxis=NULL, FUN=NULL){
    Reference.object <- GenomicMatrix$new()
    if(is.null(chr1) | is.null(chr2)){
        stop("chr1 and chr2 keys cannot be empty!")
    }
    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    ChromosomeList <- ChromInfo[,"chr"]
    if(class(chr1)!="character" | class(chr2)!="character"){
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

    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root, chr1, chr2))
    Vector <- ._Lego_Get_Something_(Group.path = Group.path, Lego = Lego, Name = Reference.object$hdf.matrix.name,
    Start = Start, Stride = Stride, Count = Count, return.what = "data")
    if(is.null(FUN)){
        return(Vector)
    }else{
        return(FUN(Vector))
    }
}