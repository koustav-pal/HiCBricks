GenomicMatrix <- R6Class("GenomicMatrix",
    public = list (
        initialize = function(){

        },
        hdf.matrices.root = "Base.matrices",
        hdf.ranges.root = "Base.ranges",
        hdf.metadata.root = "Base.metadata",
        metadata.chrom.dataset = "chrominfo",
        hdf.matrix.name = "matrix",
        hdf.matrix.coverage = "chr1_bin_coverage",
        hdf.matrix.coverage.t = "chr2_bin_coverage",
        hdf.matrix.rowSums = "chr1_row_sums",
        hdf.matrix.colSums = "chr2_col_sums",
        hdf.matrix.sparsity = "sparsity",
        hdf.bintable.ranges.group = "Bintable",
        hdf.ranges.dataset.name = "ranges",
        hdf.ranges.offset.name = "offset",
        hdf.ranges.chr.name = "chr.names",
        hdf.ranges.lengths.name = "lengths",
        mcool.resolutions.name = "resolutions",
        cache.basename = "HiCBricks",
        cache.metacol.cols = "HiCBricks_MetaInfo",
        cache.metacol.dbid = "rid",
        cache.metacol.hashname = "hashname",
        cache.metacol.dirpath = "dirpath",
        brick.extension = "brick",
        Max.vector.size=104857600,
        brick.config.name = "HiCBricks_builder_config.json",
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
                self$hdf.matrix.coverage.t,
                self$hdf.matrix.rowSums,
                self$hdf.matrix.colSums, 
                self$hdf.matrix.sparsity)
            names(Temp) <- c("chr1_bin_coverage",
                "chr2_bin_coverage",
                "chr1_row_sums",
                "chr2_col_sums",
                "sparse")
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
        },
        Configurator_JSON_headers_names = c("file_prefix", 
            "project_directory", "experiment_name", "resolutions", 
            "chromosomes", "lengths", "queues"),
        Configurator_JSON_matrix_names = c("chrom1", 
            "chrom2", "resolution", "dimensions", "mat_type", "filename")
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

._format_by_charlen <- function(string = NULL, length = 80){
    Space_locations <- str_locate_all(pattern = " ", paste(string," "))[[1]]
    Space_locations_ends <- Space_locations[,"end"]
    Binned_locations <- rle(ceiling(Space_locations_ends/length))
    Break_loc_end <- Space_locations_ends[cumsum(Binned_locations$lengths)]
    Break_loc_start <- 1
    if(length(Break_loc_end) > 1){
        Break_loc_start <- c(Break_loc_start, 
            Break_loc_end[-length(Break_loc_end)]+1)
    }
    Trimmed_string <- str_trim(str_sub(string, 
        start = Break_loc_start, 
        end = Break_loc_end))
    Finished_message <- paste(Trimmed_string, collapse = "\n")
}


._GenerateRandomName_ <- function(seed = NULL, 
    algo = c("sha1", "crc32", "md5", "sha256", 
        "xxhash32", "xxhash64", "murmur32")){
    if(is.null(seed)){
        seed <- format(Sys.time(),"%Y-%m-%d_%H-%M-%S")
    }
    algo <- match.arg(algo)
    HashString <- digest(seed, algo)
    return(HashString)
}

._Get_Brick_hashname <- function(Brick = NULL){
    Filename <- basename(Brick)
    Dir.path <- normalizePath(dirname(Brick))
    Brick.path <- file.path(Dir.path,Filename)
    Tracked_name <- ._GenerateRandomName_(Brick.path)
    return(Tracked_name)
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

._Brick_Get_Something_ <- function(Group.path = NULL, Brick = NULL, 
    Name = NULL, Index = NULL, Start = NULL, Stride = NULL, Count = NULL, 
    return.what = "group_handle"){
    Reference.object <- GenomicMatrix$new()
    Group.Handle <- ReturnH5Handler(Path = Group.path, File = Brick)
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

._Brick_Put_Something_ <- function(Group.path = NULL, Brick = NULL,
    Name = NULL, data = NULL, Index = NULL, Start = NULL, Stride = NULL, 
    Count = NULL){
    Reference.object <- GenomicMatrix$new()
    Group.handler <- ._Brick_Get_Something_(Group.path = Group.path, 
        Brick = Brick, Name = Name, return.what = "group_handle")
    h5writeDataset(obj=data, h5loc=Group.handler, name=Name, 
        index=Index, start = Start, stride = Stride, count = Count)
    H5Gclose(Group.handler)
}

._Brick_do_on_ComplexSelection_ <- function(Group.path = NULL, Brick = NULL, 
    Name = NULL, Start.list = NULL, Stride.list = NULL, Count.list = NULL, 
    Data = NULL, Brick.list = NULL, do.what = "fetch"){
    ListOfArgs <- list(Start.list,Stride.list,Count.list,Brick.list)
    if(any(vapply(ListOfArgs,!is.list,TRUE))){
        stop("Start, Stride, Count, Brick must be of type list.\n")
    }
    if(unique(vapply(ListOfArgs,length,1))!=1){
        stop("Start, Stride, Count, Brick must have same length.\n")
    }
    if(!is.vector(Data)){
        stop("Data should be a vector, when working with ",
            "Start, Stride, Count, Brick.\n")
    }
}

._Brick_WriteDataFrame_ <- function(Brick = NULL, Path = NULL, name = NULL, 
    object = NULL){
    if(!(length(c(Brick,Path,name,object))>=4)){
        stop("All arguments are required!")
    }
    Brick.handler <- ._Brick_Get_Something_(Group.path = Path, Brick = Brick, 
        Name = name, return.what = "group_handle")
    h5writeDataset.data.frame(h5loc = Brick.handler, 
        obj = object,
        name = name)
    H5Gclose(Brick.handler)
}

._Brick_WriteArray_ <- function(Brick = NULL, Path = NULL, name = NULL, 
    object = NULL){
    if(!(length(c(Brick,Path,name,object))>=4)){
        stop("All arguments are required!")
    }
    Brick.handler <- ._Brick_Get_Something_(Group.path = Path, Brick = Brick, 
        Name = name, return.what = "group_handle")
    h5writeDataset.array(h5loc = Brick.handler, 
            obj = object, 
            name = name)
    H5Gclose(Brick.handler)
}


._Brick_Add_Ranges_ = function(Group.path = NULL, Brick = NULL, name = NULL, 
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
    CreateGroups(Group.path = Group.path, File = Brick)
    ._Brick_WriteDataFrame_(Brick = Brick, Path = Group.path, 
        name = Reference.object$hdf.ranges.dataset.name, 
        object = ranges.df)
    if(is.null(names(mcol.list))){
        stop("mcol.list must be a named list!\n")
    }
    for (i in seq_along(mcol.list)) {
        m.name <- names(mcol.list[i])
        MCol <- mcol.list[[i]]
        ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
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

.create_file_connection = function(Filename=NULL,mode="r"){
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
._Compute_various_col_matrix_metrics <- function(Matrix = NULL, 
    metrics.list){
    Matrix[is.na(Matrix) | is.infinite(Matrix)] <- 0
    metrics.list[["bin.coverage"]] <- 
        metrics.list[["bin.coverage"]] + colSums(Matrix > 0)
    metrics.list[["col.sums"]] <- metrics.list[["col.sums"]] +
    colSums(Matrix)
    return(metrics.list)
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
._Process_matrix_by_distance <- function(Brick = NULL, Matrix.file = NULL, 
    delim = NULL, Group.path = NULL, chr1.len = NULL, chr2.len = NULL, 
    num.rows = 2000, distance = NULL, is.sparse = NULL, 
    compute.sparsity = NULL, sparsity.bins = 100){
    Reference.object <- GenomicMatrix$new()
    if(is.sparse){
        Sparsity.bins = sparsity.bins
    }
    Handler <- .create_file_connection(Filename = Matrix.file, mode = "r")
    Start.row <- 1
    Start.col <- 1
    Row.Offset <- 0
    Col.Offset <- 0
    Path.to.file <- Brick
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

            message(paste("Inserting Data at location:",
                paste(Start[1],", ",Start[2],sep = "")))
            message("Data length: ",Count[1],"\n")
            ._Brick_Put_Something_(Group.path=Group.path, Brick = Brick, 
                Name = Reference.object$hdf.matrix.name,
                data = Matrix, Start = Start, Stride = Stride, Count = Count)
            Start.row <- Start.row + Count[1]
            Row.Offset <- Row.Offset + num.rows 
            Object.size <- object.size(Matrix)
            Matrix <- NULL
            num.rows <- ifelse((i + num.rows) >= chr1.len, 
                chr1.len - i, num.rows)
            message(paste("Loaded",humanize_size(Object.size),"of data..."))
        }
        i<-i+1
    }
    close(Handler)
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.rowSums, 
        object = Row.sums)
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.colSums, 
        object = Row.sums)
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage, 
        object = Bin.coverage)
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage.t, 
        object = Bin.coverage)
    if(compute.sparsity){
        ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
            name = Reference.object$hdf.matrix.sparsity, 
            object = Sparsity.Index)
    }
    Attributes <- Reference.object$matrices.chrom.attributes
    Attr.vals <- c(basename(Matrix.file),as.double(Matrix.range),
        as.integer(is.sparse), 
        as.integer(distance), 
        as.integer(TRUE))
    WriteAttributes(Path = Group.path, 
        File = Brick, 
        Attributes = Attributes, 
        values = Attr.vals, on = "group")
    return(TRUE)
}
._ProcessMatrix_ <- function(Brick = NULL, Matrix.file = NULL, delim = NULL, 
    Group.path = NULL, chr1.len = NULL, chr2.len = NULL, num.rows = 2000, 
    is.sparse = NULL, compute.sparsity = NULL, distance = NULL, 
    sparsity.bins = 100){
    Reference.object <- GenomicMatrix$new()
    if(is.sparse){
        Sparsity.bins = sparsity.bins
    }
    options(datatable.fread.input.cmd.message=FALSE)
    Start.row <- 1
    Start.col <- 1
    Set.col <- TRUE
    Path.to.file <- Brick
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
    Col.sums <- NULL
    Col.bin.coverage <- NULL
    Sparsity.Index <- NULL
    Iterations.number <- chr1.len / NumLines
    Iterations <- rep(NumLines,floor(Iterations.number))
    Col.metrics.list <- list("bin.coverage" = rep(0,chr2.len),
        "col.sums" = rep(0,chr2.len))
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
        Matrix <- as.matrix(fread(file = Matrix.file, sep=delim, nrows=Iter, 
            na.strings="NA", stringsAsFactors=FALSE, skip=Skip, verbose=FALSE, 
            dec=".", showProgress=TRUE))
        message("Read ",Iter," lines after Skipping ",Skip," lines")

        Metrics.list <- ._Compute_various_matrix_metrics(Matrix = Matrix, 
            compute.sparsity = compute.sparsity, sparsity.bins = sparsity.bins, 
            range = Matrix.range, distance = distance, 
            diag.position.start = Skip + 1)
        Col.metrics.list <- ._Compute_various_col_matrix_metrics(
            Matrix = Matrix, 
            metrics.list = Col.metrics.list)
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
            message("Inserting Data at location: ",Start[1])
            message("Data length: ",Count[1])
            ._Brick_Put_Something_(Group.path=Group.path, Brick = Brick, 
                Name = Reference.object$hdf.matrix.name,
                data = Cumulative.data, Start = Start, Stride = Stride, 
                Count = Count)
            Start.row <- Start.row + Count[1]
            Set.col <- TRUE
            Cumulative.data <- NULL
            message("Loaded ",Obj.size," bytes of data...")
        }
        message("Read ",(Skip+Iter)," records...")
        i<-i+1
    }
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.rowSums, 
        object = Row.sums)
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.colSums, 
        object = Col.metrics.list[["col.sums"]])
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage, 
        object = Bin.coverage)
    ._Brick_WriteArray_(Brick = Brick, 
        Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage.t, 
        object = Col.metrics.list[["bin.coverage"]]/chr1.len)
    if(compute.sparsity){
        ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
            name = Reference.object$hdf.matrix.sparsity, 
            object = Sparsity.Index)
    }
    Attributes <- Reference.object$matrices.chrom.attributes
    options(datatable.fread.input.cmd.message=FALSE)
    Attr.vals <- c(basename(Matrix.file),
        as.double(Matrix.range),
        as.integer(is.sparse),
        as.integer(distance),
        as.integer(TRUE))
    WriteAttributes(Path = Group.path, File = Brick, 
        Attributes = Attributes, 
        values = Attr.vals, 
        on = "group")
    return(TRUE)
}

.create_brick <- function(output_directory, filename, 
    chrom1, chrom2, resolution, bintable_df, 
    hdf_chunksize=NULL, remove_existing=FALSE,
    link_existing = FALSE){

    Reference.object <- GenomicMatrix$new()
    Root.folders <- Reference.object$GetRootFolders()
    Chrom_info_df <- return_chrominfo_df(bintable_df = bintable_df, 
        chromosomes = NULL)
    Dims <- c(Chrom_info_df[Chrom_info_df$chr == chrom1,"nrow"],
        Chrom_info_df[Chrom_info_df$chr == chrom2,"nrow"])
    mat_type <- ifelse(chrom1 == chrom2, "cis", "trans")
    hdf_file <- file.path(output_directory, filename)
    if(file.exists(hdf_file)){
        if(link_existing){
            return(.create_configuration_matrix_info(
                resolution = resolution,
                chrom1 = chrom1,
                chrom2 = chrom2,
                chrom1_binned_length = Dims[1],
                chrom2_binned_length = Dims[2],
                type = mat_type,
                filename = filename))
        }
        if(!remove_existing){
            stop("An HDF file by the same name already exists.",
                " Please provide remove_existing = TRUE to overwrite it,",
                " or link_existing = TRUE to add the same file",
                " to the BrickContainer")
        }
        file.remove(hdf_file)
    }
    h5createFile(hdf_file)
    for (Folder in Root.folders) {
        CreateGroups(Group.path = Create_Path(Folder), 
            File = hdf_file)
    }
    # Create metadata chromosome groups
    ._Brick_WriteDataFrame_(Brick = hdf_file, 
        Path = c(Root.folders['metadata']), 
        name = Reference.object$metadata.chrom.dataset, 
        object = Chrom_info_df)
    ._Brick_Add_Ranges_(Group.path = Create_Path(
        c(Root.folders['ranges'],
        Reference.object$hdf.bintable.ranges.group)), 
        Brick = hdf_file,
        ranges.df = bintable_df, 
        name = Reference.object$hdf.ranges.dataset.name,
        mcol.list = NULL)
    CreateGroups(Group.path = Create_Path(c(Root.folders['matrices'],
        chrom1)), File = hdf_file)
    chr2.path <- Create_Path(
        c(Root.folders['matrices'],chrom1,chrom2))
    # cat(Chrom.info.df$chr,chrom1,chrom2,"\n")
    CreateGroups(Group.path = chr2.path, File = hdf_file)
    CreateAttributes(Path = chr2.path, File = hdf_file,
        Attributes = Reference.object$matrices.chrom.attributes,
        data_types = Reference.object$matrices.chrom.attributes.dtype,
        dims = Reference.object$matrices.chrom.attributes.dims,
        maxdims = NULL,
        on = "group")
    if(is.null(hdf_chunksize)){
        hdf_chunksize <- ceiling(Dims/100)
    }
    Array.dim <- Chrom_info_df[Chrom_info_df$chr == chrom1,"nrow"]
    CreateDataset(Path = c(Root.folders['matrices'], chrom1, chrom2),
        File = hdf_file, name = Reference.object$hdf.matrix.name,
        dims = Dims, maxdims = Dims)
    return(.create_configuration_matrix_info(
        resolution = resolution,
        chrom1 = chrom1,
        chrom2 = chrom2,
        chrom1_binned_length = Dims[1],
        chrom2_binned_length = Dims[2],
        type = mat_type,
        filename = filename))
}

._GetDimensions <- function(group.path = NULL, dataset.path = NULL, 
    File = NULL, return.what = NULL){
    File.handler <- ._Brick_Get_Something_(Group.path = group.path, 
        Brick = File, Name = dataset.path, return.what = "dataset_handle")
    Dataspace <- H5Dget_space(File.handler)
    Extents <- H5Sget_simple_extent_dims(Dataspace)
    CloseH5Con(Handle = File.handler, type = "dataset")
    if(!is.null(return.what)){
        return(Extents[[return.what]])
    }else{
        return(Extents)
    }
}

# ._Get_cachedir <- function(){
#     Reference.object <- GenomicMatrix$new()
#     Cache.dir.basename <- Reference.object$cache.basename
#     Cachebasepath <- user_cache_dir()
#     Cache.dir <- BiocFileCache(cache = file.path(Cachebasepath,
#         Cache.dir.basename), ask = FALSE)
#     return(Cache.dir)
# }

# ._add_metadata_to_cache <- function(metadata.df = NULL, cache.dir = NULL){
#     Reference.object <- GenomicMatrix$new()
#     bfcmeta(cache.dir, 
#         name = Reference.object$cache.metacol.cols,
#         append = TRUE) <- metadata.df    
# }
# ._Create_new_cached_file <- function(Cache.dir = NULL, Brick.path = NULL){
#     Reference.object <- GenomicMatrix$new()
#     Filename <- basename(Brick.path)
#     Dir.path <- normalizePath(dirname(Brick.path))
#     Tracked_name <- ._Get_Brick_hashname(Brick.path)
#     Working.File <- bfcnew(x = Cache.dir, 
#         rname = Filename, 
#         ext = Reference.object$brick.extension,
#         rtype = "local")
#     DB_id <- names(Working.File)
#     Metadata.df <- data.frame(dbid = DB_id,
#         hashname = Tracked_name,
#         dirpath = Dir.path)
#     colnames(Metadata.df) <- c(
#         Reference.object$cache.metacol.dbid,
#         Reference.object$cache.metacol.hashname,
#         Reference.object$cache.metacol.dirpath)
#     ._add_metadata_to_cache(metadata.df = Metadata.df, 
#         cache.dir = Cache.dir)
#     return(Working.File)
# }
