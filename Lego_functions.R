# HDF structure
# Base.matrices <group>
#   chromsome <group>
#       chromosome <group> Attributes FileName, Min, Max, Done
#           matrix <dataset>
#           bin.coverage <dataset>
#           row.sums <dataset>
#           sparsity <dataset>
# Base.ranges <group>
#   Bintable <group> containing binary attributes Strand, Names
#       Names <dataset>
#       ranges <dataset> based on value of attributes 3,4,5 column
#       offsets <dataset>
#       lengths <dataset>
#       chr.names <dataset>
#       BlaBla1 <dataset>
#       BlaBla2 <dataset>
#   Other <group>
# Base.metadata <group>
# chromosomes <dataset>

CreateLego <- function(ChromNames=NULL, BinTable=NULL, bin.delim="\t",
    col.index=c(1,2,3), impose.discontinuity=TRUE, ChunkSize=NULL, 
    Output.Filename=NULL, exec="cat", remove.existing=FALSE,
    sparse=FALSE, sparsity.compute.bins=100){
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
                col.index = col.index, chromosomes = ChromNames,
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
        ranges.df = Bintable, mcol.list = Metadata.list)
    # Create matrices groups
    for (chrom1 in ChromosomeList) {
        CreateGroups(Group.path = Create_Path(c(Root.folders['matrices'],chrom1)), File = HDF.File)
        for (chrom2 in ChromosomeList) {
            Chrom2.path <- Create_Path(c(Root.folders['matrices'],chrom1,chrom2))
            # cat(Chrom.info.df$chr,chrom1,chrom2,"\n")
            CreateGroups(Group.path = Chrom2.path, File = HDF.File)
            CreateAttributes(Path = Chrom2.path, File = HDF.File, 
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

#### It should be done after the fist write 
Lego_List_Matrices <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Lego_Get_ChromInfo(Lego = Lego)
    chr1.list <- lapply(ChromInfo[,"chr"], function(chr1){
        chr2.list <- lapply(ChromInfo[,"chr"], function(chr2){
            Colnames <- Reference.object$matrices.chrom.attributes
            Values <- GetAttributes(Path = Create_Path(c(Reference.object$hdf.matrices.root, chr1, chr2)),
                File = Lego, Attributes = Colnames, on = "group")
            temp.df <- data.frame(chr1,chr2,cbind(Values))
            colnames(temp.df) <- c("chr1","chr2",Colnames)
            temp.df
        })
        chr2.df <- do.call(rbind,chr2.list)
    })
    Matrix.list.df <- do.call(rbind,chr1.list)
    return(Matrix.list.df)
}

Lego_MakeGRangesObject = function(Chrom=NULL, Start=NULL, End=NULL, Strand=NULL, Names=NULL){
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

Lego_add_ranges <- function(Lego = NULL, ranges = NULL, name = NULL){
    Reference.object <- GenomicMatrix$new()
    if(!(class(Ranges) %in% "GRanges") | ("list" %in% class(Ranges))){
        stop("Object of class Ranges expected")
    }
    Ranges.df <- as.data.frame(Ranges)
    if(is.unsorted(Ranges.df$seqnames)){
        stop("Ranges must be sorted by chromosome!")
    }
    Metadata.Cols <- names(Ranges.df)[,c(4:ncol(Ranges.df))]
    Metadata.list <- lapply(Metadata.Cols,function(x){
        Ranges.df[,x]
    })
    names(Metadata.list) <- Metadata.Cols
    ._Lego_Add_Ranges_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root,name)), Lego = Lego, 
        ranges.df = Ranges.df, mcol.list = Metadata.list)
}

Lego_Get_ChromInfo <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Dataset <- ._Lego_Get_Something_(Group.path = Reference.object$hdf.metadata.root, 
        Lego = Lego, Name = Reference.object$metadata.chrom.dataset, handler = FALSE)
    return(Dataset)
}

Lego_Get_Ranges <- function(Lego = NULL, chr = NULL, rangekey = NULL, as.ranges = FALSE, attach_cols = NULL){
    Reference.object <- GenomicMatrix$new()
    if(is.null(rangekey) | is.null(Lego)){
        stop("rangekey and Lego cannot remain empty!\n")
    }
    if(!Lego_RangeKey_exists(Lego = Lego, rangekey = rangekey)){
        stop("rangekey not found!")
    }
    Start <- NULL
    Stride <- NULL
    Count <- NULL
    Index <- NULL
    if(!is.null(chr)){
        chromosomes <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.chr.name, return.what = "data")
        if(!(chr %in% ChromInfo.df$chr)){
            stop("chr not found!")
        }        
        Starts <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.offset.name, return.what = "data")
        Lengths <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.lengths.name, return.what = "data")
        Which.one <- chromosomes[chromosomes == chr]
        Start <- c(Starts[Which.one],1)
        Stride <- c(1,1)
        Count <- c(Lengths[Which.one],3)
    }
    Dataset <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.dataset.name, Start = Start, Stride = Stride,
        Count = Count, return.what = "data")
    if(as.ranges){
        Dataset <- MakeGRangesObject(Chrom = Dataset[,'chr'], Start = Dataset[,'start'], End = Dataset[,'end'])
    }
    return(Dataset)
}

Lego_Get_Bintable <- function(Lego = NULL, chr = NULL, as.ranges = TRUE){
    Table <- Lego_Get_Ranges(Lego = Lego, chr = chr, rangekey = Reference.object$hdf.bintable.ranges.group, as.ranges = as.ranges)
    return(Table)
}

Lego_ListRangeKeys <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Handler <- ._Lego_Get_Something_(Group.path = Create_Path(Reference.object$hdf.ranges.root), 
        Lego = Lego, return.what = "group_handle")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)
    return(GroupList)
}

Lego_RangeKey_exists <- function(Lego = NULL, rangekey = NULL){
    Keys <- Lego_ListRangeKeys(Lego = Lego)
    return(rangekey %in% Keys)
}

Lego_list_Ranges_mcols <- function(Lego = NULL, rangekey = NULL){
    RangeKeys <- Lego_ListRangeKeys(Lego = Lego)
    if(!is.null(rangekey)){
        if(!Lego_RangeKey_exists(Lego = Lego, rangekey = rangekey)){
            stop("rangekey not found!")
        }
        RangeKeys <- RangeKeys[RangeKeys %in% rangekey]
    }
    mcol.list <- lapply(RangeKeys,function(x){
        Handler <- ._Lego_Get_Something_(Group.path = Create_Path(Reference.object$hdf.ranges.root,x),File = Lego, return.what = "group_handle")
        GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)
        data.frame(rangekey = x, m.col = GroupList)
    })
    mcol.df <- do.call(rbind,mcol.list)
    mcol.df <- mcol.df[!(mcol.df$m.col %in% Reference.object$hdf.ranges.protected.names()),]
    if(nrow(mcol.df)==0){
        mcol.df <- NA
    }
    return(mcol.df)
}

Lego_LoadMatrix <- function(Lego = NULL, chr1 = NULL, chr2 = NULL, FormatAs = "mxnMatrix", 
    delim = " ", Dataset = NULL, exec = NULL){
    ListVars <- list(Lego = Lego, FormatAs = FormatAs, chr1 = chr1, chr2 = chr2, exec = NULL, delim = delim, Dataset = Dataset)
    sapply(1:length(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    sapply(1:length(ListVars[c("Lego","chr1","chr2","Dataset")]),function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })
}

Lego_matrix_exists <- function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    return(chr1 %in% Matrix.list$chr1 & chr2 %in% Matrix.list$chr2)
}


Lego_matrix_isDone <- function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    return(Matrix.list[Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2, "done"])
}


Lego_matrix_minmax <- function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- c(Matrix.list[Filter, "min"],Matrix.list[Filter, "max"])
    return(Extent)
}


Lego_matrix_filename <- function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- Matrix.list[Filter, "filename"]
    return(Extent)
}



Lego_LoadMatrix <- function(Lego = NULL, LOCDIR = NULL, chr1 = NULL, chr2 = NULL, Basename = NULL, 
    create.dir = FALSE, create.recursively = FALSE,  exec = NULL, Dataset = NULL, 
    delim = " ", remove.prior = FALSE, num.rows = 2000, is.sparse = FALSE, sparsity.bins = 100){

    ListVars <- list(Lego = Lego, LOCDIR = LOCDIR, create.dir = create.dir, create.recursively = create.recursively, 
        Basename = Basename, chr1 = chr1, chr2 = chr2, is.sparse = is.sparse, sparsity.bins = sparsity.bins, 
        exec = exec, delim = delim, Dataset = Dataset, remove.prior = remove.prior)
    sapply(1:length(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    sapply(1:length(ListVars[c("Lego","LOCDIR","chr1","chr2","Basename","Dataset")]),function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })
    Chrom.info.df <- Lego_GetChromInfo(Lego = Lego)
    if(!(all(c(chr1, chr2) %in% Chrom.info.df[,"chr"]))){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    Chrom1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr1,"nrow"]
    Chrom2.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr2,"nrow"]
    _ProcessMatrix_(Data = Dataset, delim = NULL, Matrix.file = Matrix.file.path, 
        exec = exec, chr1.len = Chrom1.len, chr2.len = Chrom2.len, 
        fix.num.rows.at = num.rows, is.sparse = is.sparse, sparsity.bins = sparsity.bins)
}



Lego_GetValuesByDistance <- function(Lego = NULL, chr = NULL, distance  = NULL,
    constrain.region=NULL,batch.size=500,FUN=NULL){
    if(any(sapply(list(Lego,chr1,chr2,distance),is.null))) {
        stop("Lego, chr, distance cannot be NULL.\n")
    }
    Reference.object <- GenomicMatrix$new()
    ChromInfo <- Lego_Get_ChromInfo(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr, chr2 = chr)){
        stop("Chromosome is not listed in this HDF file.")
    }
    if(any(sapply(list(Lego,chr,distance),length) > 1)){
        stop("Lego, chr and distance can only have values of length 1.\n")
    }
    Nrow <- (ChromInfo[ChromInfo$chr == chr,"nrow"]-1)
    if(any(distance > Nrow,distance < 0)){
        stop("Distance must range between 0 and",Nrow,"\n")   
    }
    Root.folders <- Reference.object$GetRootFolders()
    Path <- Create_Path(c(Root.folders['matrices'],chr,chr))
    DatasetHandle <- ._Lego_Get_Something_(Group.path = Path, Lego = Lego, return.what = "dataset_handle")
    Vector.start <- 1
    Vector.stop <- ChromInfo[ChromInfo$chr==chr,"nrow"]
    if(!is.null(Constrain.region)){
        Vector.coordinates <- ReturnRegionPosition(Constrain.region=Constrain.region,Chrom=Chrom)
        if(is.null(Vector.coordinates)){
            stop("Overlap operation was unsuccessful! Please check coordinates ",Constrain.region)
        }
        Vector.start <- min(Vector.coordinates)
        Vector.stop <- max(Vector.coordinates)
    }
    Starting.col <- Vector.start + Distance
    Count <- ((Vector.stop - Vector.start) + 1) - Distance 
    tot.mat.extracted <- Count * Count
    # cat("Boo ",Count,"\n")
    Start <- c(Vector.start,Starting.col)
    Stride <- c(1,1)
    Counts <- Count
    CumSums <- 0
    Groups <- c(private$hdf.matrices,Chrom)
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
        # cat(Count,Offset,cur.start,"\n")
        diag(private$FetchFromDataset(Groups=Groups,Chrom=Chrom,Start=cur.start,Stride=c(1,1),Count=c(Count,Count)))
    })

    DistancesVector <- do.call(c,DistancesVector.list)

    # Preparing Start stride and
}

# AddAttribute(Key=private$Bintable.Key,Value=normalizePath(BinTable))


# if(sparse){
#     if(length(sparsity.compute.bins)>1 | !is.numeric(sparsity.compute.bins)){
#         stop("sparsity.compute.bins expects numeric value of length 1\n")
#     }
#     private$ComputeSparsity=sparse
#     private$Sparsity.compute.bins=sparsity.compute.bins
# }