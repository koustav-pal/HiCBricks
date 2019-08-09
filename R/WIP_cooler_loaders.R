# Read 1000000 rows at a time
# Create_Indices

Create_read_indexes <- function(table.handle, delim, big.seek = 1000000,
    chr1s, chr2s, chr1.col, chr2.col, value.col, Starts, Ends){
    message("Indexing the file...")
    i == 1
    while(i == 1) {
        Table.df <- read.table(file = table.handle, sep = sep, 
            nrows = big.seek)
        Table.df <- Table.df[,c(chr1.col, chr2.col, value.col)]
        colnames(Table.df) <- c("chr1", "chr2", "value")
        if(any(Table.df[,"chr2"] < Table.df[,"chr1"])){
            stop(paste("Provided chr2.col was not larger",
                "than chr1.col!","This file is not in",
                "upper.tri sparse format!", sep = " "))
        }
        chr1.pos <- NULL
        chr2.pos <- NULL
        for (i in seq_along(chr1s)) {
            chr1 <- chr1s[i]
            chr2 <- chr2s[i]
            chr1.start <- Starts[chr1]
            chr2.start <- Starts[chr2]
            chr1.end <- Ends[chr1]
            chr2.end <- Ends[chr2]
            chr1.filter <- Table.df[,chr1.col] >= chr1.start & 
                                Table.df[,chr1.col] < chr1.end
            chr2.filter <- Table.df[,chr2.col] >= chr2.start & 
                                Table.df[,chr2.col] < chr2.end
            chr1.pos <- min(which(chr1.filter))
            if(!any(chr1.filter)){
                break
            }
            chr2.pos <- min()
        }
    }
}


._ProcessTable_ <- function(Brick, table.handle, delim, 
    Group.path, exec, num.rows = 1000000, is.sparse, 
    chr1.col, chr2.col, value.col, chr1.start, chr2.start, 
    chr1.end, chr2.end, compute.sparsity, distance,
    sep, sparsity.bins = 100){

    Reference.object <- GenomicMatrix$new()

    if(is.sparse){
        Sparsity.bins = sparsity.bins
    }
    
    options(datatable.fread.input.cmd.message=FALSE)

    Command <- paste(exec, table.file)    
    Path.to.file <- Brick
    Cumulative.data <- NULL
    Cumulative.indices <- NULL
    Matrix.range <- c(NA,NA)
    Bin.coverage <- NULL
    Row.sums <- NULL
    Sparsity.Index <- NULL

    while(chr1.start <= chr1.end & 
        chr2.start <= chr2.end) {
        Table.df <- read.table(file = table.handle, sep = sep, nrows = num.rows)
        Table.df <- Table.df[,c(chr1.col, chr2.col, value.col)]
        colnames(Table.df) <- c("chr1", "chr2", "value")
        if(!all(Table.df[,"chr1"] %in% seq(chr1.start, chr1.end))){
            Which <- min(which(!(Table.df[,"chr1"] %in% seq(chr1.start, chr1.end))))

        }
    }

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

    ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
        name = Reference.object$hdf.matrix.rowSums, object = Row.sums)
    ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
        name = Reference.object$hdf.matrix.coverage, object = Bin.coverage)
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