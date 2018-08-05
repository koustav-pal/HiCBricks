require(R6)
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
        Max.vector.size=104857600,
        hdf.ranges.protected.names = function(){
            Protect <- c(self$hdf.ranges.dataset.name, self$hdf.ranges.lengths.name, 
                self$hdf.ranges.chr.name, self$hdf.ranges.offset.name)
            return(Protect)
        },
        matrices.chrom.attributes = c("filename","min","max","sparsity","done"),
        matrices.chrom.attributes.dtype = c("character","double","double","integer","integer"),
        matrices.chrom.attributes.dims = list(1,1,1,1,1),
        bintable.attributes = c("stranded","names"),
        ranges.bintable.dataset = "bintable",
        GeneralFileSeparator = "_",
        Ranges.separator=":",
        NonStrandedColNames=c("chr","start","end"),
        TerrificNumberOfHiCFormats = c("NxNMatrix","PAIRIX","Cooler","HOMER","CustomTable","AggregateFragments"),
        GetRootFolders = function() {
            Folders <- c(self$hdf.matrices.root, self$hdf.ranges.root, self$hdf.metadata.root)
            names(Folders) <- c('matrices','ranges','metadata')
            return(Folders)
        },
        GetBaseRangesFolders = function(){
            Folders <- c(self$hdf.bintable.ranges.group)
            names(Folders) <- c('Bintable')
            return(Folders)            
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

._GenerateRandomName_ <- function(){
    library("digest")
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
._Do_on_vector_SparsityIndex_ = function(x=NULL,index=NULL,sparsity.bins = NULL,length=NULL){
    x[is.na(x) | is.infinite(x)] <- 0
    Range <- (index-sparsity.bins):(index+sparsity.bins)
    Range <- Range[Range>0 & Range<length]
    Rows <- x[Range]
    return(length(Rows[Rows!=0])/length(Rows))
}
._Lego_Get_Something_ <- function(Group.path = NULL, Lego = NULL, Name = NULL,
    Index = NULL, Start = NULL, Stride = NULL, Count = NULL, return.what = "group_handle"){
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
        Dataset <- h5read(name = Name, file = Group.Handle, index = Index, start = Start,
            stride = Stride, count = Count)
        H5Gclose(Group.Handle)
        return(Dataset)
    }
}
._Lego_Put_Something_ <- function(Group.path = NULL, Lego = NULL, Name = NULL, data = NULL,
    Index = NULL, Start = NULL, Stride = NULL, Count = NULL, Block = NULL){
    Reference.object <- GenomicMatrix$new()
    Group.handler <- ._Lego_Get_Something_(Group.path = Group.path, Lego = Lego, Name = Name, return.what = "group_handle")
    h5writeDataset(obj=data, h5loc=Group.handler, name=Name, index=Index, start = Start, stride = Stride, count = Count)
    H5Gclose(Group.handler)
}
._Lego_do_on_ComplexSelection_ <- function(Group.path = NULL, Lego = NULL, Name = NULL,
    Start.list = NULL, Stride.list = NULL, Count.list = NULL, 
    Block.list = NULL, do.what = "fetch"){
    ListOfArgs <- list(Start.list,Stride.list,Count.list,Block.list)
    if(any(sapply(ListOfArgs,!is.list))){
        stop("Start, Stride, Count, Block must be of type list.\n")
    }
    if(unique(sapply(ListOfArgs,length))!=1){
        stop("Start, Stride, Count, Block must have same length.\n")
    }
    if(!is.vector(data)){
        stop("data should be a vector, when working with Start, Stride, Count, Block.\n")
    }
}
._Lego_WriteDataFrame_ <- function(Lego = NULL, Path = NULL, name = NULL, object = NULL){
    library(stringr)
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
._Lego_WriteArray_ <- function(Lego = NULL, Path = NULL, name = NULL, object = NULL){
    library(stringr)
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
._Lego_Add_Ranges_ = function(Group.path = NULL, Lego = NULL, name = NULL, ranges.df = NULL, mcol.list = NULL){
    Reference.object <- GenomicMatrix$new()
    ChrOffsetCols <- Reference.object$hdf.ranges.protected.names()
    ChrOffsetCols <- ChrOffsetCols[!(ChrOffsetCols %in% Reference.object$hdf.ranges.dataset.name)]
    if(any(!(ChrOffsetCols %in% names(mcol.list)))) {
        Chrom.lengths <- get_chrom_info(bin.table = ranges.df, FUN = length, col.name = 'chr')
        Chrom.sizes <- get_chrom_info(bin.table = ranges.df, FUN = max, col.name = 'end')
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
    ._Lego_WriteDataFrame_(Lego = Lego, Path = Group.path, name = Reference.object$hdf.ranges.dataset.name, object = ranges.df)
    if(is.null(names(mcol.list))){
        stop("mcol.list must be a named list!\n")
    }
    for (i in seq_along(mcol.list)) {
        m.name <- names(mcol.list[i])
        MCol <- mcol.list[[i]]
        ._Lego_WriteArray_(Lego = Lego, Path = Group.path, name = m.name, object = MCol)
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
._ProcessMatrix_ <- function(Read.file = NULL, delim = NULL, exec = NULL, Group.path = NULL, dataset.name = NULL,
    chr1.len = NULL, chr2.len = NULL, num.rows = 2000, is.sparse = NULL, distance = NULL, sparsity.bins = NULL){
    require(data.table)
    Reference.object <- GenomicMatrix$new()
    if(is.sparse){
        Sparsity.bins = sparsity.bins
    }
    Command <- paste(exec,Read.file,sep=" ")
    Start.row <- 0
    Path.to.file <- file.path(private$Output.Directory,private$Output.Filename)
    Cumulative.data <- NULL
    Cumulative.distances.data <- NULL
    Cumulative.indices <- NULL
    Matrix.range <- c(NA,NA)
    NumLines <- Reference.object$FindLineNumbers(Row.len=chr1.len,Col.len=chr2.len)
    if(NumLines <= fix.num.rows.at){
        NumLines <- fix.num.rows.at
    }
    if(Chrom1.len <= NumLines){
        NumLines <- Chrom1.len
    }
    Bin.coverage <- NULL
    Row.sums <- NULL
    Sparsity.Index <- NULL
    Iterations.number <- Chrom1.len / NumLines
    Iterations <- rep(NumLines,floor(Iterations.number))
    if(floor(Iterations.number)!=ceiling(Iterations.number)){
        cumulative <- sum(Iterations)
        Iterations <- c(Iterations,(Chrom1.len-cumulative))
    }
    Skippity<-0
    if(length(Iterations)>1){
        Skippity.cumsum <- cumsum(Iterations)
        Skippity <- c(0,Skippity.cumsum[1:(length(Skippity.cumsum)-1)])
    }
    i<-1
    while(i<=length(Iterations)) {
        Iter <- Iterations[i]
        Skip <- Skippity[i]
        Matrix <- as.matrix(fread(input=Command, sep=delim, nrows=Iter, na.strings="NA", 
            stringsAsFactors=FALSE, skip=Skip, verbose=FALSE, dec=".", showProgress=FALSE))
        cat("Read",Iter,"lines after Skipping",Skip,"lines\n")
        Bin.coverage <- c(Bin.coverage,sapply(1:nrow(Matrix),function(x){
            Vec.sub <- Matrix[x,]
            ._Do_on_vector_PercentGTZero_(Vec.sub)
        }))
        Row.sums <- c(Row.sums,sapply(1:nrow(Matrix),function(x){
            Vec.sub <- Matrix[x,]
            ._Do_on_vector_ComputeRowSums_(Vec.sub)
        }))
        Sparsity.Index <- c(Sparsity.Index,sapply(1:nrow(Matrix),function(x){
            Vec.sub <- Matrix[x,]
            sparsity.bin.idexes <- x + Skip
            ._Do_on_vector_SparsityIndex_(x=Vec.sub,index=sparsity.bin.idexes,length=Chrom2.len)
        }))
        Row.extent <- ._Do_on_vector_ComputeMinMax_(Matrix)
        if(Matrix.range[1] > Row.extent[1] | is.na(Matrix.range[1])) {
            Matrix.range[1] <- Row.extent[1]
        }
        if(Matrix.range[2] < Row.extent[2] | is.na(Matrix.range[2])){
            Matrix.range[2] <- Row.extent[2]
        }
        Cumulative.data <- rbind(Cumulative.data,Matrix)
        Obj.size <- object.size(Cumulative.data)
        if(Obj.size >= Reference.object$Max.vector.size | i == length(Iterations)){
            Start <- c(Start.row+1,1) 
            Stride <- c(1,1)
            Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
            cat("Inserting Data at location:",Start[1],"\n")
            cat("Data length:",Count[1],"\n")
            ._Lego_Put_Something_(Group.path=Group.path, Lego = Lego, Name = Reference.object$hdf.matrix.name,
                data = Cumulative.data, Start = Start, Stride = Stride, Count = Count)
            Start.row <- Start.row+Count[1]
            Cumulative.data <- NULL
            cat("Loaded ",Obj.size," bytes of data...\n")
        }
        cat("Read ",(Skip+Iter),"records...\n")
        i<-i+1
    }
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, name = Reference.object$hdf.matrix.rowSums, object = Row.sums)
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, name = Reference.object$hdf.matrix.coverage, object = Bin.coverage)
    ._Lego_WriteArray_(Lego = Lego, Path = Group.path, name = Reference.object$hdf.matrix.sparsity, object = Sparsity.Index)
    Attributes <- Reference.object$matrices.chrom.attributes
    Attr.vals <- c(basename(Read.file),Matrix.range,is.sparse,TRUE)
    WriteAttributes(Path = Group.path, File = Lego, Attributes = Attributes, values = Attr.vals, on = "group")
}








# ._ProcessMatrix_ = function(Read.file = NULL, delim=" ", exec= " ", fix.num.rows.at=2000,
#     is.spa){
#     require(data.table)

#     Chromosomes.all <- c(Chromosome1,Chromosome2)
#     if(is.null(exec)){
#         stop("exec takes as input the shell command to be used by data.table for reading")
#     }
#     if((!is.null(Chromosome1) & !is.null(Chromosome2))){
#         private$CheckForChromosomes(Chrom1=Chromosome1,Chrom2=Chromosome2)
#     }
#     if(!is.null(Filename)){
#         Command <- paste(exec,Filename,sep=" ")
#         Start.row <- 0
#         Path.to.file <- file.path(private$Output.Directory,private$Output.Filename)
#         Chrom1.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chromosome1)
#         Chrom2.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chromosome2)
#         Chrom1.len <- length(Chrom1.ranges)
#         Chrom2.len <- length(Chrom2.ranges)
#         Chrom1.Group <- private$ReturnH5GroupHandler(Groups=c(private$hdf.matrices,Chromosome1))
#         Cumulative.data <- NULL
#         Cumulative.distances.data <- NULL
#         Cumulative.indices <- NULL
#         Matrix.range <- c(NA,NA)
#         NumLines <- private$FindLineNumbers(Row.len=Chrom1.len,Col.len=Chrom2.len)
#         if(NumLines <= fix.num.rows.at){
#             NumLines <- fix.num.rows.at
#         }
#         if(Chrom1.len <= NumLines){
#             NumLines <- Chrom1.len
#         }
#         Bin.coverage <- NULL
#         Row.sums <- NULL
#         Sparsity.Index <- NULL
#         Iterations.number <- Chrom1.len / NumLines
#         Iterations <- rep(NumLines,floor(Iterations.number))
#         if(floor(Iterations.number)!=ceiling(Iterations.number)){
#             cumulative <- sum(Iterations)
#             Iterations <- c(Iterations,(Chrom1.len-cumulative))
#         }
#         Skippity<-0
#         if(length(Iterations)>1){
#             Skippity.cumsum <- cumsum(Iterations)
#             Skippity <- c(0,Skippity.cumsum[1:(length(Skippity.cumsum)-1)])
#         }
#         i<-1
#         while(i<=length(Iterations)) {
#             Iter <- Iterations[i]
#             Skip <- Skippity[i]
#             Matrix <- as.matrix(fread(input=Command, sep=delim, nrows=Iter, na.strings="NA", 
#                 stringsAsFactors=FALSE, skip=Skip,verbose=FALSE, dec=".",
#                 showProgress=FALSE))
#             cat("Read",Iter,"lines after Skipping",Skip,"lines\n")
#             Bin.coverage <- c(Bin.coverage,sapply(1:nrow(Matrix),function(x){
#                 Vec.sub <- Matrix[x,]
#                 PercentGTZero(Vec.sub)
#             }))
#             Row.sums <- c(Row.sums,sapply(1:nrow(Matrix),function(x){
#                 Vec.sub <- Matrix[x,]
#                 ComputeRowSums(Vec.sub)
#             }))
#             if(self$is.sparse() & Chromosome1==Chromosome2){
#                 Sparsity.Index <- c(Sparsity.Index,sapply(1:nrow(Matrix),function(x){
#                     Vec.sub <- Matrix[x,]
#                     sparsity.bin.idexes <- x + Skip
#                     SparsityIndex(x=Vec.sub,index=sparsity.bin.idexes,length=Chrom2.len)
#                 }))
#             }
#             Row.extent <- ComputeMinMax(Matrix)
#             if(Matrix.range[1] > Row.extent[1] | is.na(Matrix.range[1])) {
#                 Matrix.range[1] <- Row.extent[1]
#             }
#             if(Matrix.range[2] > Row.extent[2] | is.na(Matrix.range[2])){
#                Matrix.range[2] <- Row.extent[2]
#             }
#             Cumulative.data <- rbind(Cumulative.data,Matrix)
#             Obj.size <- object.size(Cumulative.data)
#             if(Obj.size>=private$Max.vector.size | i==length(Iterations)){
#                 Start <- c(Start.row+1,1) 
#                 Stride <- c(1,1)
#                 Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
#                 cat("Inserting Data at location:",Start[1],"\n")
#                 cat("Data length:",Count[1],"\n")
#                 private$InsertIntoDataset(Connection=Chrom1.Group,
#                     Chrom=Chromosome2,Data=Cumulative.data,Start=Start,Stride=Stride,Count=Count)
#                 Start.row <- Start.row+Count[1]
#                 Cumulative.data <- NULL
#                 cat("Loaded ",Obj.size," bytes of data...\n")
                
#             }
#             cat("Read ",(Skip+Iter),"records...\n")
#             i<-i+1
#         }
#         H5Gclose(Chrom1.Group)
#         private$CloseH5FileConnection()
#         # cat(length(Bin.coverage),"\n")
#         private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Bin.coverage,
#             Col.Name=paste(Chromosome2,"cov",sep="."),Replace=TRUE,na.function=as.numeric)
#         # cat(length(Row.sums),"\n")
#         private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Row.sums,
#             Col.Name=paste(Chromosome2,"sum",sep="."),Replace=TRUE,na.function=as.numeric)
#         if((Chromosome1 == Chromosome2) & self$is.sparse()){
#             private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Sparsity.Index,
#                 Col.Name=paste("sparsity","idx",sep="."),Replace=TRUE,na.function=as.numeric)                    
#         }
#         private$AddFileToList(Key1=Chromosome1,Key2=Chromosome2,Value=normalizePath(Filename))
#         private$AddDoneForChromosome(Key1=Chromosome1,Key2=Chromosome2)
#     }
# }
