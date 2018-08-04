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
        ranges.df = Bintable, name = Reference.object$hdf.ranges.dataset.name, mcol.list = NULL)
    # # Create matrices groups
    # for (chrom1 in ChromosomeList) {
    #     CreateGroups(Group.path = Create_Path(c(Root.folders['matrices'],chrom1)), File = HDF.File)
    #     for (chrom2 in ChromosomeList) {
    #         chr2.path <- Create_Path(c(Root.folders['matrices'],chrom1,chrom2))
    #         # cat(Chrom.info.df$chr,chrom1,chrom2,"\n")
    #         CreateGroups(Group.path = chr2.path, File = HDF.File)
    #         CreateAttributes(Path = chr2.path, File = HDF.File, 
    #             Attributes = Reference.object$matrices.chrom.attributes,
    #             data_types = Reference.object$matrices.chrom.attributes.dtype,
    #             dims = Reference.object$matrices.chrom.attributes.dims,
    #             maxdims = NULL,
    #             on = "group")
    #         Dims <- c(Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"], Chrom.info.df[Chrom.info.df$chr == chrom2,"nrow"])
    #         if(is.null(ChunkSize)){
    #             ChunkSize <- ceiling(Dims/100)
    #         }
    #         Array.dim <-Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"]
    #         CreateDataset(Path = c(Root.folders['matrices'],chrom1,chrom2), File = HDF.File, 
    #             name = Reference.object$hdf.matrix.name, dims = Dims, maxdims = Dims)
    #     }
    # }
}

Lego_get_chrominfo <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Dataset <- ._Lego_Get_Something_(Group.path = Reference.object$hdf.metadata.root, 
        Lego = Lego, Name = Reference.object$metadata.chrom.dataset, handler = FALSE)
    return(Dataset)
}

Lego_make_ranges = function(Chrom=NULL, Start=NULL, End=NULL, Strand=NULL, Names=NULL){
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

Lego_add_ranges = function(Lego = NULL, ranges = NULL, name = NULL){
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
    Metadata.Cols <- names(Ranges.df)[,c(4:ncol(Ranges.df))]
    Metadata.list <- lapply(Metadata.Cols,function(x){
        Ranges.df[,x]
    })
    names(Metadata.list) <- Metadata.Cols
    ._Lego_Add_Ranges_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root,name)), Lego = Lego, 
        ranges.df = Ranges.df, mcol.list = Metadata.list)
}

Lego_list_rangekeys = function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Handler <- ._Lego_Get_Something_(Group.path = Create_Path(Reference.object$hdf.ranges.root), 
        Lego = Lego, return.what = "group_handle")
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    return(GroupList)
}

Lego_rangekey_exists = function(Lego = NULL, rangekey = NULL){
    Keys <- Lego_list_rangekeys(Lego = Lego)
    return(rangekey %in% Keys)
}

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
        cat(Start,"\n")
        cat(Stride,"\n")
        cat(Count,"\n")
    }
    Dataset <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, rangekey)),
        Lego = Lego, Name = Reference.object$hdf.ranges.dataset.name, Start = Start, Stride = Stride,
        Count = Count, return.what = "data")
    Dataset <- Lego_make_ranges(Chrom = Dataset[,'chr'], Start = Dataset[,'start'], End = Dataset[,'end'])

    MCols <- Lego_list_ranges_mcols(Lego = Lego, rangekey = rangekey)
    if(!is.na(MCols)){
        MCols.col <- MCols[,"m.col"]
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

Lego_get_bintable = function(Lego = NULL, chr = NULL){
    Table <- Lego_get_ranges(Lego = Lego, chr = chr, rangekey = Reference.object$hdf.bintable.ranges.group)
    return(Table)
}


Lego_fetch_range_index = function(Lego = NULL, chr = NULL, start = NULL, end = NULL,names = NULL,type = "any"){
    AllTypes<-c("any","within")
    if( any(!(type %in% AllTypes)) ){
        stop("type takes one of two arguments: c(\"any\",\"within\")")
    }
    if(is.null(chr) | is.null(start) | is.null(end)){
        stop("Chrom, start, end cannot be empty")
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
        Cur.Chrom <- Chrom[Filter]
        Cur.Start <- Start[Filter]
        Cur.end <- end[Filter]
        Cur.Names <- Names[Filter]
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
    names(OverlapByChromosome)<-UniqueChromosomes
    return(OverlapByChromosome)
}

Lego_return_region_position = function(region=NULL,chr=NULL){
    if(!is.character(region) | length(region) > 1){
        stop("region must be a character vector of length 1")
    }
    Coord.Split<- Split_genomic_coordinates(Coordinate=region)
    Region.Chrom<-Coord.Split[[1]][1]
    Region.start<-as.numeric(Coord.Split[[1]][2])
    Region.stop<-as.numeric(Coord.Split[[1]][3])
    if(chr!=Region.Chrom){
        stop("region chr and provided chr must be same.") 
    }
    Region.Ranges<-Lego_fetch_range_index(Chrom=chr, Start=Region.start, end=Region.stop, type="within")
    Vector.coordinates <- Region.Ranges[[chr]][[1]][["Indexes"]]
    return(Vector.coordinates)
}

Lego_load_matrix = function(Lego = NULL, chr1 = NULL, chr2 = NULL, matrix.file = NULL, delim = " ", exec = NULL,  
    remove.prior = FALSE, num.rows = 2000, is.sparse = FALSE, distance = NULL, sparsity.bins = 100){
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
    if(!Lego_matrix_exists(Lego = NULL, chr1 = chr1, chr2 = chr2)){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    if(Lego_matrix_isdone(Lego = NULL, chr1 = chr1, chr2 = chr2) && !remove.prior){
        stop("A matrix was preloaded before. Use remove.priori = TRUE to force value replacement\n")
    }
    Chrom.info.df <- Lego_GetChromInfo(Lego = Lego)
    chr1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr1,"nrow"]
    chr2.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr2,"nrow"]
    ._ProcessMatrix_(Data = Dataset, delim = delim, Matrix.file = matrix.file, 
        exec = exec, chr1.len = chr1.len, chr2.len = chr2.len, 
        fix.num.rows.at = num.rows, is.sparse = is.sparse, distance = distance,
        sparsity.bins = sparsity.bins)
}

Lego_matrix_exists = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    return(chr1 %in% Matrix.list$chr1 & chr2 %in% Matrix.list$chr2)
}

Lego_matrix_isdone = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    return(Matrix.list[Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2, "done"])
}

Lego_matrix_minmax = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- c(Matrix.list[Filter, "min"],Matrix.list[Filter, "max"])
    return(Extent)
}

Lego_matrix_filename = function(Lego = NULL, chr1 = NULL, chr2 = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.list <- Lego_List_Matrices(Lego = Lego)
    if(!Lego_matrix_exists(Lego = Lego, chr1 = chr1, chr2 = chr2)){
        stop("chr1 chr2 pairs were not found\n")
    }
    Filter <- Matrix.list$chr1 == chr1 & Matrix.list$chr2 == chr2
    Extent <- Matrix.list[Filter, "filename"]
    return(Extent)
}

Lego_get_values_by_distance = function(Lego = NULL, chr = NULL, distance  = NULL,
    constrain.region=NULL,batch.size=500,FUN=NULL){
    if(any(sapply(list(Lego,chr1,chr2,distance),is.null))) {
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
        stop("Distance must range between 0 and",Nrow,"\n")
    }
    Root.folders <- Reference.object$GetRootFolders()
    Path <- Create_Path(c(Root.folders['matrices'],chr,chr))
    Vector.start <- 1
    Vector.stop <- ChromInfo[ChromInfo$chr==chr,"nrow"]
    if(!is.null(constrain.region)){
        Vector.coordinates <- Lego_return_region_position(region=constrain.region, chr=chr)
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
        diag(._Lego_Get_Something_(Group.path = Path, Lego = Lego, Start= cur.start, Stride = Stride, Count = c(Count,Count), return.what = "data"))
    })

    DistancesVector <- do.call(c,DistancesVector.list)
    if(is.null(FUN)){
        return(DistancesVector)
    }else{
        return(FUN(DistancesVector))
    }
    # Preparing Start stride and
}

FetchMatrixWithinCoords = function(Lego = NULL, x.coords=NULL, y.coords=NULL, FUN=NULL){
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
    chr1.len <- ChromInfo$nrow[ChromInfo$chr == chr1,]
    chr2.len <- ChromInfo$nrow[ChromInfo$chr == chr2,]
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

Fetch_row_or_col_vector = function(Lego = NULL, chr1=NULL, chr2=NULL, by=NULL, vector=NULL, regions=NULL, FUN=NULL){
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
        Values <- Lego_get_vector_values(Lego = Lego, chr1=chr1, chr2=chr2, xaxis=x, yaxis=y, FUN=FUN)
        return(Values)
    }) 
    return(Vector.values)
}

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

    Group.path <- Create_Path(Reference.object$hdf.matrices.root, chr1, chr2)
    Vector <- ._Lego_Get_Something_(Group.path = Group.path, Lego = Lego, Name = Reference.object$hdf.matrix.name,
    Start = Start, Stride = Stride, Count = Count, return.what = "data")
    if(is.null(FUN)){
        return(Vector)
    }else{
        return(FUN(Vector))
    }
}