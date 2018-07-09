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
        hdf.bintable.ranges.group = "Bintable",
        hdf.ranges.dataset.name = "ranges",
        matrices.chrom.attributes = c("filename","done"),
        ranges.bintable.dataset = "bintable",
        GeneralFileSeparator = "_",
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
        Max.vector.size=104857600,
        Matrix.range=NA,
        Num.lines=1,
        hdf.root.folders=c("matrices","base.ranges.tables"),
        Protected.Ranges.Keys=c("Bintable"),
        Bintable.Key = "Bintable",
        NonStrandedColNames=c("chr","start","end"),
        Matrice.done = NA,
        Ranges.separator=":"
    )
)

._GenerateRandomName_ <- function(){
    library("digest")
    Seed <- format(Sys.time(),"%Y-%m-%d_%H-%M-%S")
    HashString <- digest(Seed,"crc32")
    return(HashString)
}

._ProcessMatrix_ <- function(Read.file = NULL, delim = NULL, exec = NULL, DatasetHandle = NULL,
    chr1.len = NULL, chr2.len = NULL, fix.num.rows.at = NULL, is.sparse = NULL, sparsity.bins = NULL){
    require(data.table)
    Reference.object <- GenomicMatrix$new()
    PercentGTZero = function(x){
        x[is.na(x) | is.infinite(x)] <- 0
        LengthOfRow <- length(x)
        LengthOfGTZero <- length(x[x!=0])
        Fraction <- LengthOfGTZero/LengthOfRow
        return(Fraction)
    }
    ComputeRowSums = function(x){
        x[is.na(x) | is.infinite(x)] <- 0
        return(sum(x))
    }
    ComputeMinMax = function(x){
        x[is.na(x) | is.infinite(x)] <- 0
        return(c(min(x),max(x)))
    }
    if(is.sparse){
        Sparsity.bins = sparsity.bins
        SparsityIndex = function(x=NULL,index=NULL,length=NULL){
            x[is.na(x) | is.infinite(x)] <- 0
            Range <- (index-Sparsity.bins):(index+Sparsity.bins)
            Range <- Range[Range>0 & Range<length]
            Rows <- x[Range]
            return(length(Rows[Rows!=0])/length(Rows))
        }
    }
    Command <- paste(exec,Read.file,sep=" ")
    Start.row <- 0
    Path.to.file <- file.path(private$Output.Directory,private$Output.Filename)
    Chrom1.Group <- private$ReturnH5GroupHandler(Groups=c(private$hdf.matrices,Chromosome1))
    Cumulative.data <- NULL
    Cumulative.distances.data <- NULL
    Cumulative.indices <- NULL
    NumLines <- private$FindLineNumbers(Row.len=chr1.len,Col.len=chr2.len)
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
            stringsAsFactors=FALSE, skip=Skip,verbose=FALSE, dec=".",
            showProgress=FALSE))
        cat("Read",Iter,"lines after Skipping",Skip,"lines\n")
        Bin.coverage <- c(Bin.coverage,sapply(1:nrow(Matrix),function(x){
            Vec.sub <- Matrix[x,]
            PercentGTZero(Vec.sub)
        }))
        Row.sums <- c(Row.sums,sapply(1:nrow(Matrix),function(x){
            Vec.sub <- Matrix[x,]
            ComputeRowSums(Vec.sub)
        }))
        if(self$is.sparse() & Chromosome1==Chromosome2){
            Sparsity.Index <- c(Sparsity.Index,sapply(1:nrow(Matrix),function(x){
                Vec.sub <- Matrix[x,]
                sparsity.bin.idexes <- x + Skip
                SparsityIndex(x=Vec.sub,index=sparsity.bin.idexes,length=Chrom2.len)
            }))
        }
        Row.extent <- ComputeMinMax(Matrix)
        if(private$Matrix.range[[Chromosome1]][[Chromosome2]][1] > Row.extent[1]) {
            private$Matrix.range[[Chromosome1]][[Chromosome2]][1] <- Row.extent[1]
        }
        if(private$Matrix.range[[Chromosome1]][[Chromosome2]][2] < Row.extent[2]){
            private$Matrix.range[[Chromosome1]][[Chromosome2]][2] <- Row.extent[2]  
        }
        Cumulative.data <- rbind(Cumulative.data,Matrix)
        Obj.size <- object.size(Cumulative.data)
        if(Obj.size>=private$Max.vector.size | i==length(Iterations)){
            Start <- c(Start.row+1,1) 
            Stride <- c(1,1)
            Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
            cat("Inserting Data at location:",Start[1],"\n")
            cat("Data length:",Count[1],"\n")
            private$InsertIntoDataset(Connection=Chrom1.Group,
                Chrom=Chromosome2,Data=Cumulative.data,Start=Start,Stride=Stride,Count=Count)
            Start.row <- Start.row+Count[1]
            Cumulative.data <- NULL
            cat("Loaded ",Obj.size," bytes of data...\n")
        }
        cat("Read ",(Skip+Iter),"records...\n")
        i<-i+1
        H5Gclose(Chrom1.Group)
        private$CloseH5FileConnection()
        # cat(length(Bin.coverage),"\n")
        private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Bin.coverage,
            Col.Name=paste(Chromosome2,"cov",sep="."),Replace=TRUE,na.function=as.numeric)
        # cat(length(Row.sums),"\n")
        private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Row.sums,
            Col.Name=paste(Chromosome2,"sum",sep="."),Replace=TRUE,na.function=as.numeric)
        if((Chromosome1 == Chromosome2) & self$is.sparse()){
            private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Sparsity.Index,
                Col.Name=paste("sparsity","idx",sep="."),Replace=TRUE,na.function=as.numeric)                    
        }
        private$AddFileToList(Key1=Chromosome1,Key2=Chromosome2,Value=normalizePath(Filename))
        private$AddDoneForChromosome(Key1=Chromosome1,Key2=Chromosome2)
    }
}

._ProcessMatrix_ = function(Read.file = NULL, delim=" ", exec= " ", fix.num.rows.at=2000,
    is.spa){
    require(data.table)

    Chromosomes.all <- c(Chromosome1,Chromosome2)
    if(is.null(exec)){
        stop("exec takes as input the shell command to be used by data.table for reading")
    }
    if((!is.null(Chromosome1) & !is.null(Chromosome2))){
        private$CheckForChromosomes(Chrom1=Chromosome1,Chrom2=Chromosome2)
    }
    if(!is.null(Filename)){
        Command <- paste(exec,Filename,sep=" ")
        Start.row <- 0
        Path.to.file <- file.path(private$Output.Directory,private$Output.Filename)
        Chrom1.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chromosome1)
        Chrom2.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chromosome2)
        Chrom1.len <- length(Chrom1.ranges)
        Chrom2.len <- length(Chrom2.ranges)
        Chrom1.Group <- private$ReturnH5GroupHandler(Groups=c(private$hdf.matrices,Chromosome1))
        Cumulative.data <- NULL
        Cumulative.distances.data <- NULL
        Cumulative.indices <- NULL
        NumLines <- private$FindLineNumbers(Row.len=Chrom1.len,Col.len=Chrom2.len)
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
                stringsAsFactors=FALSE, skip=Skip,verbose=FALSE, dec=".",
                showProgress=FALSE))
            cat("Read",Iter,"lines after Skipping",Skip,"lines\n")
            Bin.coverage <- c(Bin.coverage,sapply(1:nrow(Matrix),function(x){
                Vec.sub <- Matrix[x,]
                PercentGTZero(Vec.sub)
            }))
            Row.sums <- c(Row.sums,sapply(1:nrow(Matrix),function(x){
                Vec.sub <- Matrix[x,]
                ComputeRowSums(Vec.sub)
            }))
            if(self$is.sparse() & Chromosome1==Chromosome2){
                Sparsity.Index <- c(Sparsity.Index,sapply(1:nrow(Matrix),function(x){
                    Vec.sub <- Matrix[x,]
                    sparsity.bin.idexes <- x + Skip
                    SparsityIndex(x=Vec.sub,index=sparsity.bin.idexes,length=Chrom2.len)
                }))
            }
            Row.extent <- ComputeMinMax(Matrix)
            if(private$Matrix.range[[Chromosome1]][[Chromosome2]][1] > Row.extent[1]) {
                private$Matrix.range[[Chromosome1]][[Chromosome2]][1] <- Row.extent[1]
            }
            if(private$Matrix.range[[Chromosome1]][[Chromosome2]][2] < Row.extent[2]){
                private$Matrix.range[[Chromosome1]][[Chromosome2]][2] <- Row.extent[2]  
            }
            Cumulative.data <- rbind(Cumulative.data,Matrix)
            Obj.size <- object.size(Cumulative.data)
            if(Obj.size>=private$Max.vector.size | i==length(Iterations)){
                Start <- c(Start.row+1,1) 
                Stride <- c(1,1)
                Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
                cat("Inserting Data at location:",Start[1],"\n")
                cat("Data length:",Count[1],"\n")
                private$InsertIntoDataset(Connection=Chrom1.Group,
                    Chrom=Chromosome2,Data=Cumulative.data,Start=Start,Stride=Stride,Count=Count)
                Start.row <- Start.row+Count[1]
                Cumulative.data <- NULL
                cat("Loaded ",Obj.size," bytes of data...\n")
                
            }
            cat("Read ",(Skip+Iter),"records...\n")
            i<-i+1
        }
        H5Gclose(Chrom1.Group)
        private$CloseH5FileConnection()
        # cat(length(Bin.coverage),"\n")
        private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Bin.coverage,
            Col.Name=paste(Chromosome2,"cov",sep="."),Replace=TRUE,na.function=as.numeric)
        # cat(length(Row.sums),"\n")
        private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Row.sums,
            Col.Name=paste(Chromosome2,"sum",sep="."),Replace=TRUE,na.function=as.numeric)
        if((Chromosome1 == Chromosome2) & self$is.sparse()){
            private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Sparsity.Index,
                Col.Name=paste("sparsity","idx",sep="."),Replace=TRUE,na.function=as.numeric)                    
        }
        private$AddFileToList(Key1=Chromosome1,Key2=Chromosome2,Value=normalizePath(Filename))
        private$AddDoneForChromosome(Key1=Chromosome1,Key2=Chromosome2)
    }
}
