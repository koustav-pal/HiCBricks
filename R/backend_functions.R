GenomicMatrix <- R6Class("GenomicMatrix",
    public = list (
        initialize = function(){

        },
        hdf.matrices.root = "Base.matrices",
        hdf.ranges.root = "Base.ranges",
        hdf.metadata.root = "Base.metadata",
        metadata.chrom.dataset = "chrominfo",
        hdf.matrix.name = "matrix",
        hdf.matrix.coverage = "bin.coverage",
        hdf.matrix.rowSums = "row.sums",
        hdf.matrix.sparsity = "sparsity",
        hdf.bintable.ranges.group = "Bintable",
        hdf.ranges.dataset.name = "ranges",
        hdf.ranges.offset.name = "offset",
        hdf.ranges.chr.name = "chr.names",
        hdf.ranges.lengths.name = "lengths",
        mcool.resolutions.name = "resolutions",
        cache.basename = "HiCLegos",
        lego.extension = "lego",
        Max.vector.size=104857600,
        mcool.available.normalisations = function(){
            Names <- c("Knight-Ruitz","Vanilla-coverage",
                "Vanilla-coverage-square-root","Iterative-Correction")
            Vector <- c("KR","VC","VC_SQRT","weight")
            names(Vector) <- Names
            return(Vector)
        },
        genomic.ranges.protected.names = c("seqnames", "strand", "seqlevels", 
            "seqlengths", "isCircular", "start", "end", "width", "element"),
        genomic.ranges.metadata.functions = function(name = NULL){
            transformlist <- list(
                "strand" = strand, 
                "seqlevels" = seqlevels, 
                "seqlengths" = seqlengths, 
                "width" = width)
            if(is.null(name)){
                return(names(transformlist))
            }
            return(transformlist[[name]])
        },
        hdf.matrix.meta.cols = function(){
            Temp <- c(self$hdf.matrix.coverage, 
                self$hdf.matrix.rowSums, 
                self$hdf.matrix.sparsity)
            names(Temp) <- c("bin.cov","row.sums","sparse")
            return(Temp)
        },
        hdf.ranges.protected.names = function(){
            Protect <- c(self$hdf.ranges.dataset.name, 
                self$hdf.ranges.lengths.name, 
                self$hdf.ranges.chr.name, self$hdf.ranges.offset.name)
            return(Protect)
        },
        matrices.chrom.attributes = c("filename","min","max","sparsity",
            "distance","done"),
        matrices.chrom.attributes.dtype = c("character","double","double",
            "integer","integer","integer"),
        matrices.chrom.attributes.fun.cast = function(type = NULL){
            as.logical.as.integer = function(x){
                as.logical(as.integer(x))
            }
            TypeList <- list("filename" = as.character, 
                "min" = as.double, 
                "max" = as.double, 
                "sparsity" = as.logical.as.integer, 
                "distance" = as.integer, 
                "done" = as.logical.as.integer)
            return(TypeList[[type]])
        },
        matrices.chrom.attributes.dims = list(1,1,1,1,1,1),
        bintable.attributes = c("stranded","names"),
        ranges.bintable.dataset = "bintable",
        GeneralFileSeparator = "_",
        Ranges.separator=":",
        NonStrandedColNames=c("chr","start","end"),
        TerrificNumberOfHiCFormats = c("NxNMatrix","PAIRIX","Cooler","HOMER",
            "CustomTable","AggregateFragments"),
        GetRootFolders = function() {
            Folders <- c(self$hdf.matrices.root, 
                self$hdf.ranges.root, self$hdf.metadata.root)
            names(Folders) <- c('matrices','ranges','metadata')
            return(Folders)
        },
        GetBaseRangesFolders = function(){
            Folders <- c(self$hdf.bintable.ranges.group)
            names(Folders) <- c('Bintable')
            return(Folders)            
        },
        mcool.bintable.keys = function(version = NULL){
            if(version <= 1){
                return(c("bins","chrom_id","start","end"))
            }
            if(version >= 2){
                return(c("bins","chrom","start","end"))
            }
        },
        mcool.scaffold.keys = function(version = NULL){
            if(version < 1){
                return(c("scaffolds","length","name"))
            }
            if(version >= 2){
                return(c("chroms","length","name"))
            }
        },
        mcool.matrix.keys = function(version = NULL){
            if(version <1){
                return(c("matrix","bin1_id","bin2_id","count"))
            }
            if(version >=2){
                return(c("pixels","bin1_id","bin2_id","count"))
            }
        },
        mcool.index.keys = function(){
            return(c("indexes","bin1_offset","chrom_offset"))
        }
    ),
    private = list(
        Attribute.List=NA,
        File.List=NA,
        RangesObjects=NA,
        Ranges.Keys=NA,
        Ranges.Col.Keys=NA,
        WhichStarter=NA,
        Working.Dir=NA, 
        Working.File=NA,
        HDF.Connection=NA,
        ComputeSparsity=FALSE,
        Sparsity.compute.bins=NA,
        Matrix.range=NA,
        Num.lines=1,
        hdf.root.folders=c("matrices","base.ranges.tables"),
        Protected.Ranges.Keys=c("Bintable"),
        Bintable.Key = "Bintable",
        Matrice.done = NA
        
    )
)

._Check_numeric <- function(x){
    return(is.numeric(x) | is.integer(x))
}


._GenerateRandomName_ <- function(){
    Seed <- format(Sys.time(),"%Y-%m-%d_%H-%M-%S")
    HashString <- digest(Seed,"crc32")
    return(HashString)
}
._Do_on_vector_PercentGTZero_ = function(x){
    x[is.na(x) | is.infinite(x)] <- 0
    LengthOfRow <- length(x)
    LengthOfGTZero <- length(x[x!=0])
    Fraction <- LengthOfGTZero/LengthOfRow
    return(Fraction)
}
._Do_on_vector_ComputeRowSums_ = function(x){
    x[is.na(x) | is.infinite(x)] <- 0
    return(sum(x))
}
._Do_on_vector_ComputeMinMax_ = function(x){
    x[is.na(x) | is.infinite(x)] <- 0
    return(c(min(x),max(x)))
}
._Do_on_vector_SparsityIndex_ = function(x=NULL,index=NULL,
    sparsity.bins = NULL){
    x[is.na(x) | is.infinite(x)] <- 0
    Range <- (index-sparsity.bins):(index+sparsity.bins)
    Range <- Range[Range>0 & Range<length(x)]
    Rows <- x[Range]
    return(length(Rows[Rows!=0])/length(Rows))
}
._Lego_Get_Something_ <- function(Group.path = NULL, Lego = NULL, Name = NULL,
    Index = NULL, Start = NULL, Stride = NULL, Count = NULL, 
    return.what = "group_handle"){
    Reference.object <- GenomicMatrix$new()
    Group.Handle <- ReturnH5Handler(Path = Group.path, File = Lego)
    if(return.what == "group_handle"){
        return(Group.Handle)
    }
    if(is.null(Name)){
        stop("Name cannot be NULL\n")
    }
    if(return.what == "dataset_handle"){
        Dataset.Handle <- H5Dopen(name = Name, h5loc = Group.Handle)
        H5Gclose(Group.Handle)
        return(Dataset.Handle)
    }
    if(return.what == "data"){
        Dataset <- h5read(name = Name, file = Group.Handle, 
            index = Index, start = Start,
            stride = Stride, count = Count, callGeneric = FALSE)
        H5Gclose(Group.Handle)
        return(Dataset)
    }
}
._Lego_Put_Something_ <- function(Group.path = NULL, Lego = NULL,
    Name = NULL, data = NULL, Index = NULL, Start = NULL, Stride = NULL, 
    Count = NULL, Block = NULL){
    Reference.object <- GenomicMatrix$new()
    Group.handler <- ._Lego_Get_Something_(Group.path = Group.path, 
        Lego = Lego, Name = Name, return.what = "group_handle")
    h5writeDataset(obj=data, h5loc=Group.handler, name=Name, 
        index=Index, start = Start, stride = Stride, count = Count)
    H5Gclose(Group.handler)
}
._Lego_do_on_ComplexSelection_ <- function(Group.path = NULL, Lego = NULL, 
    Name = NULL, Start.list = NULL, Stride.list = NULL, Count.list = NULL, 
    Data = NULL, Block.list = NULL, do.what = "fetch"){
    ListOfArgs <- list(Start.list,Stride.list,Count.list,Block.list)
    if(any(vapply(ListOfArgs,!is.list,TRUE))){
        stop("Start, Stride, Count, Block must be of type list.\n")
    }
    if(unique(vapply(ListOfArgs,length,1))!=1){
        stop("Start, Stride, Count, Block must have same length.\n")
    }
    if(!is.vector(Data)){
        stop("Data should be a vector, when working with ",
            "Start, Stride, Count, Block.\n")
    }
}
._Lego_WriteDataFrame_ <- function(Lego = NULL, Path = NULL, name = NULL, 
    object = NULL){
    if(!(length(c(Lego,Path,name,object))>=4)){
        stop("All arguments are required!")
    }
    Lego.handler <- ._Lego_Get_Something_(Group.path = Path, Lego = Lego, 
        Name = name, return.what = "group_handle")
    h5writeDataset.data.frame(h5loc = Lego.handler, 
        obj = object,
        name = name)
    H5Gclose(Lego.handler)
}
._Lego_WriteArray_ <- function(Lego = NULL, Path = NULL, name = NULL, 
    object = NULL){
    if(!(length(c(Lego,Path,name,object))>=4)){
        stop("All arguments are required!")
    }
    Lego.handler <- ._Lego_Get_Something_(Group.path = Path, Lego = Lego, 
        Name = name, return.what = "group_handle")
    h5writeDataset.array(h5loc = Lego.handler, 
            obj = object, 
            name = name)
    H5Gclose(Lego.handler)
}
._Lego_Add_Ranges_ = function(Group.path = NULL, Lego = NULL, name = NULL, 
    ranges.df = NULL, mcol.list = NULL){
    Reference.object <- GenomicMatrix$new()
    ChrOffsetCols <- Reference.object$hdf.ranges.protected.names()
    Fltr <- !(ChrOffsetCols %in% Reference.object$hdf.ranges.dataset.name)
    ChrOffsetCols <- ChrOffsetCols[Fltr]
    if(any(!(ChrOffsetCols %in% names(mcol.list)))) {
        Chrom.lengths <- get_chrom_info(bin.table = ranges.df, 
            FUN = length, col.name = 'chr')
        Chrom.sizes <- get_chrom_info(bin.table = ranges.df, 
            FUN = max, col.name = 'end')
        Chrom.info.df <- data.frame(chr = names(Chrom.lengths),
            nrow = as.vector(Chrom.lengths),
            size = as.vector(Chrom.sizes),stringsAsFactors = FALSE)
        CumSums <- cumsum(Chrom.info.df[,"nrow"])
        Starts <- c(1,CumSums[-length(CumSums)]+1)
        Temp.list <- list(
            Starts,
            Chrom.info.df[,"nrow"],
            Chrom.info.df[,"chr"]
            )
        names(Temp.list) <- c(Reference.object$hdf.ranges.offset.name,
            Reference.object$hdf.ranges.lengths.name,
            Reference.object$hdf.ranges.chr.name)
        mcol.list <- c(mcol.list,Temp.list)
    }
    CreateGroups(Group.path = Group.path, File = Lego)
    ._Lego_WriteDataFrame_(Lego = Lego, Path = Group.path, 
        name = Reference.object$hdf.ranges.dataset.name, 
        object = ranges.df)
    if(is.null(names(mcol.list))){
        stop("mcol.list must be a named list!\n")
    }
    for (i in seq_along(mcol.list)) {
        m.name <- names(mcol.list[i])
        MCol <- mcol.list[[i]]
        ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
            name = m.name, object = MCol)
    }
}
._FindLineNumbers_ = function(Row.len=NULL,Col.len=NULL){
    Reference.object <- GenomicMatrix$new()
    pixel.mem <- 48
    row.pixel.mem <- pixel.mem * Col.len
    Batch.size <-  floor(Reference.object$Max.vector.size/row.pixel.mem) 
    Batch.size[Batch.size<1] <- 1
    return(Batch.size)
}
._Do_rbind_on_matrices_of_different_sizes_ <- function(Matrix.top = NULL, 
    Matrix.bottom = NULL, row.length = NULL, col.length = NULL, 
    top.coords = NULL, bottom.coords = NULL){
    if(is.null(Matrix.top)){
        return(Matrix.bottom)
    }
    if(bottom.coords[2] != top.coords[2]){
        Makeup.col <- bottom.coords[2] - top.coords[2]
        Matrix.top <- cbind(Matrix.top,matrix(NA, 
            nrow = nrow(Matrix.top), ncol = Makeup.col))
    }
    if(bottom.coords[1] != top.coords[1]){
        Makeup.col <- bottom.coords[1] - top.coords[1]
        Matrix.bottom <- cbind(matrix(NA, 
            nrow = nrow(Matrix.bottom), 
            ncol = Makeup.col),Matrix.bottom)   
    }
    Matrix.return <- rbind(Matrix.top,Matrix.bottom)
    return(Matrix.return)
}
._Create_file_connection = function(Filename=NULL,mode="r"){
    if(grepl('$/.gz',Filename)){
        connection=gzfile(Filename,open=mode)
    }else if(grepl('$/.bz2',Filename)){
        connection=bzfile(Filename,open=mode)
        }else {
            connection=file(Filename,mode)
        }
        return(connection)
}
._Compute_various_matrix_metrics <- function(Matrix = NULL, 
    compute.sparsity=FALSE, sparsity.bins = 100, range = NULL, 
    distance = NULL, diag.position.start = NULL){
    sparsity.bins <- ifelse(sparsity.bins < distance, sparsity.bins, distance)
    Bin.coverage <- vapply(seq_len(nrow(Matrix)),function(x){
        Vec.sub <- Matrix[x,]
        ._Do_on_vector_PercentGTZero_(Vec.sub)
    },1)
    Row.sums <- vapply(seq_len(nrow(Matrix)),function(x){
        Vec.sub <- Matrix[x,]
        ._Do_on_vector_ComputeRowSums_(Vec.sub)
    },1)
    if(compute.sparsity){
        Sparsity.Index <- vapply(seq_len(nrow(Matrix)),function(x){
            Vec.sub <- Matrix[x,]
            sparsity.bin.idexes <- sparsity.bins
            ._Do_on_vector_SparsityIndex_(x=Vec.sub,
                index=diag.position.start + (x - 1),
                sparsity.bins = sparsity.bins)
        },1)
    }else{
        Sparsity.Index <- NULL
    }
    Row.extent <- ._Do_on_vector_ComputeMinMax_(Matrix)
    if(range[1] > Row.extent[1] | is.na(range[1])) {
        range[1] <- Row.extent[1]
    }
    if(range[2] < Row.extent[2] | is.na(range[2])){
        range[2] <- Row.extent[2]
    }
    A.list <- list("bin.cov" = Bin.coverage, "row.sum" = Row.sums, 
        "sparsity" = Sparsity.Index, "extent" = range)
    return(A.list)
}
humanize_size <- function(x){
    Size <- x/1024
    if(Size < 1){
        return(paste(x,"bytes"))
    }
    Size <- x/1024/1024
    if(Size < 1){
        return(paste(round(x/1024,2),"KB"))
    }
    Size <- x/1024/1024/1024
    if(Size < 1){
        return(paste(round(x/1024/1024,2),"MB"))
    }
    return(paste(round(x/1024/1024/1024,2),"GB"))
}
._Process_matrix_by_distance <- function(Lego = NULL, Matrix.file = NULL, 
    delim = NULL, Group.path = NULL, chr1.len = NULL, chr2.len = NULL, 
    num.rows = 2000, distance = NULL, is.sparse = NULL, compute.sparsity = NULL,
    sparsity.bins = 100){
    Reference.object <- GenomicMatrix$new()
    if(is.sparse){
        Sparsity.bins = sparsity.bins
    }
    Handler <- ._Create_file_connection(Filename = Matrix.file, mode = "r")
    Start.row <- 1
    Start.col <- 1
    Row.Offset <- 0
    Col.Offset <- 0
    Path.to.file <- Lego
    Matrix.range <- c(NA,NA)
    if(chr1.len <= num.rows){
        num.rows <- chr1.len
    }
    Bin.coverage <- NULL
    Row.sums <- NULL
    Sparsity.Index <- NULL
    Matrix <- NULL
    i <- 1
    while(i<=chr1.len) {
        Col.upper.limit <- ifelse((i + distance) > chr2.len, 
            chr2.len, i + distance)
        Col.lower.limit <- ifelse((i - distance) <= 0, 1, i - distance)
        if(is.null(Matrix)){
            Start.col <- Col.lower.limit
            Col.Offset <- Col.lower.limit - 1
            num.row.upper.limit <- ifelse((i + num.rows + distance) > chr2.len, 
                chr2.len, i + num.rows + distance)
            Matrix <- matrix(data = 0, nrow = num.rows, 
                ncol = (num.row.upper.limit - Col.lower.limit + 1))
        }
        Vector <- scan(file=Handler, what=double(), sep=delim, 
            nlines=1, quiet=TRUE)
        Row.loc <- (i - Row.Offset)
        Col.loc <- c((Col.lower.limit-Col.Offset):(Col.upper.limit-Col.Offset))
        Matrix[Row.loc, Col.loc] <- Vector[Col.lower.limit:Col.upper.limit]
        if((i - Row.Offset) == num.rows){
            Start <- c(Start.row,Start.col)
            Stride <- c(1,1)
            Count <- c(nrow(Matrix),ncol(Matrix))
            Metrics.list <- ._Compute_various_matrix_metrics(Matrix = Matrix, 
                compute.sparsity = compute.sparsity, 
                sparsity.bins = sparsity.bins, range = Matrix.range, 
                distance = distance, diag.position.start = i - Col.Offset)
            Matrix.range <- Metrics.list[["extent"]]
            Bin.coverage <- c(Bin.coverage,Metrics.list[["bin.cov"]])
            Row.sums <- c(Row.sums,Metrics.list[["row.sum"]])
            Sparsity.Index <- c(Sparsity.Index,Metrics.list[["sparsity"]])

            message("Inserting Data at location:",Start[1],Start[2],"\n")
            message("Data length:",Count[1],"\n")
            ._Lego_Put_Something_(Group.path=Group.path, Lego = Lego, 
                Name = Reference.object$hdf.matrix.name,
                data = Matrix, Start = Start, Stride = Stride, Count = Count)
            Start.row <- Start.row + Count[1]
            Row.Offset <- Row.Offset + num.rows 
            Object.size <- object.size(Matrix)
            Matrix <- NULL
            num.rows <- ifelse((i + num.rows) >= chr1.len, 
                chr1.len - i, num.rows)
            message("Loaded",humanize_size(Object.size),"of data...\n")
        }
        i<-i+1
    }
    close(Handler)
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
        name = Reference.object$hdf.matrix.rowSums, object = Row.sums)
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage, object = Bin.coverage)
    if(compute.sparsity){
        ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
            name = Reference.object$hdf.matrix.sparsity, 
            object = Sparsity.Index)
    }
    Attributes <- Reference.object$matrices.chrom.attributes
    Attr.vals <- c(basename(Matrix.file),as.double(Matrix.range),
        as.integer(is.sparse), 
        as.integer(distance), 
        as.integer(TRUE))
    WriteAttributes(Path = Group.path, 
        File = Lego, 
        Attributes = Attributes, 
        values = Attr.vals, on = "group")
    return(TRUE)
}
._ProcessMatrix_ <- function(Lego = NULL, Matrix.file = NULL, delim = NULL, 
    exec = NULL, Group.path = NULL, chr1.len = NULL, chr2.len = NULL, 
    num.rows = 2000, is.sparse = NULL, compute.sparsity = NULL,
    distance = NULL, sparsity.bins = 100){

    Reference.object <- GenomicMatrix$new()
    if(is.sparse){
        Sparsity.bins = sparsity.bins
    }
    Command <- paste(exec,Matrix.file,sep=" ")
    Start.row <- 1
    Start.col <- 1
    Set.col <- TRUE
    Path.to.file <- Lego
    Cumulative.data <- NULL
    Cumulative.indices <- NULL
    Matrix.range <- c(NA,NA)
    NumLines <- ._FindLineNumbers_(Row.len=chr1.len,Col.len=chr2.len)
    if(NumLines <= num.rows){
        NumLines <- num.rows
    }
    if(chr1.len <= NumLines){
        NumLines <- chr1.len
    }
    Bin.coverage <- NULL
    Row.sums <- NULL
    Sparsity.Index <- NULL
    Iterations.number <- chr1.len / NumLines
    Iterations <- rep(NumLines,floor(Iterations.number))
    if(is.null(distance)){
        distance <- chr2.len
    }
    if(floor(Iterations.number)!=ceiling(Iterations.number)){
        cumulative <- sum(Iterations)
        Iterations <- c(Iterations,(chr1.len-cumulative))
    }
    Skippity<-0
    if(length(Iterations)>1){
        Skippity.cumsum <- cumsum(Iterations)
        Skippity <- c(0,Skippity.cumsum[-length(Skippity.cumsum)])
    }
    i<-1
    # top.coords <- NULL
    # bottom.coords <- NULL
    Drop.what <- seq_len(chr2.len)

    while(i<=length(Iterations)) {
        Iter <- Iterations[i]
        Skip <- Skippity[i]
        Matrix <- as.matrix(fread(input=Command, sep=delim, nrows=Iter, 
            na.strings="NA", stringsAsFactors=FALSE, skip=Skip, verbose=FALSE, 
            dec=".", showProgress=TRUE))
        message("Read",Iter,"lines after Skipping",Skip,"lines\n")

        Metrics.list <- ._Compute_various_matrix_metrics(Matrix = Matrix, 
            compute.sparsity = compute.sparsity, sparsity.bins = sparsity.bins, 
            range = Matrix.range, distance = distance, 
            diag.position.start = Skip + 1)
        Matrix.range <- Metrics.list[["extent"]]
        Bin.coverage <- c(Bin.coverage,Metrics.list[["bin.cov"]])
        Row.sums <- c(Row.sums,Metrics.list[["row.sum"]])
        Sparsity.Index <- c(Sparsity.Index,Metrics.list[["sparsity"]])

        Cumulative.data <- rbind(Cumulative.data,Matrix)
        Obj.size <- object.size(Cumulative.data)
        if(Obj.size >= Reference.object$Max.vector.size | 
            i == length(Iterations)){
            Start <- c(Start.row,1)
            Stride <- c(1,1)
            Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
            message("Inserting Data at location:",Start[1],"\n")
            message("Data length:",Count[1],"\n")
            ._Lego_Put_Something_(Group.path=Group.path, Lego = Lego, 
                Name = Reference.object$hdf.matrix.name,
                data = Cumulative.data, Start = Start, Stride = Stride, 
                Count = Count)
            Start.row <- Start.row + Count[1]
            Set.col <- TRUE
            Cumulative.data <- NULL
            message("Loaded ",Obj.size," bytes of data...\n")
        }
        message("Read ",(Skip+Iter),"records...\n")
        i<-i+1
    }
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
        name = Reference.object$hdf.matrix.rowSums, object = Row.sums)
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage, object = Bin.coverage)
    if(compute.sparsity){
        ._Lego_WriteArray_(Lego = Lego, Path = Group.path, 
            name = Reference.object$hdf.matrix.sparsity, 
            object = Sparsity.Index)
    }
    Attributes <- Reference.object$matrices.chrom.attributes
    Attr.vals <- c(basename(Matrix.file),
        as.double(Matrix.range),
        as.integer(is.sparse),
        as.integer(distance),
        as.integer(TRUE))
    WriteAttributes(Path = Group.path, File = Lego, 
        Attributes = Attributes, 
        values = Attr.vals, 
        on = "group")
    return(TRUE)
}

._GetDimensions <- function(group.path = NULL, dataset.path = NULL, 
    File = NULL, return.what = NULL){   
    File.handler <- ._Lego_Get_Something_(Group.path = group.path, Lego = File, 
        Name = dataset.path, return.what = "dataset_handle")
    Dataspace <- H5Dget_space(File.handler)
    Extents <- H5Sget_simple_extent_dims(Dataspace)
    CloseH5Con(Handle = File.handler, type = "dataset")
    if(!is.null(return.what)){
        return(Extents[[return.what]])
    }else{
        return(Extents)
    }
}

._Get_cachedir <- function(){
    Reference.object <- GenomicMatrix$new()
    Cache.dir.basename <- Reference.object$cache.basename
    Cachebasepath <- user_cache_dir()
    Cache.dir <- BiocFileCache(cache = file.path(Cachebasepath,
        Cache.dir.basename), ask = FALSE)
    return(Cache.dir)
}