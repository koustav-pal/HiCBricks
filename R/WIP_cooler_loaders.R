.create_read_indexes <- function(file, delim, big_seek = 1000000,
    chrs, chr1_col, chr2_col, value_col, starts, ends){
    i == 1
    indexes_list <- list()
    records_read <- 0
    Test_table_df <- fread(file = file, sep = sep, nrows = 10, 
                    skip = records_read, data.table = FALSE)
    if(!(ncol(Test_table_df) >= 3)){
        stop("A table in the upper triangle sparse format is",
            " expected with atleast 3 columns")
    }
    Test_table_df <- Test_table_df[,c(chr1_col, chr2_col, value_col)]
    chr1_class <- class(Test_table_df[,chr1_col])
    chr2_class <- class(Test_table_df[,chr2_col])
    value_class <- class(Test_table_df[,value_col])
    if(!all(c("numeric", "numeric", "numeric") %in% 
            c(chr1_class, chr2_class, value_class))){
        stop("Expected numeric values for chr1, chr2 bins and values.")
    }
    message("Indexing the file...")
    while(i == 1) {
        Table_df <- try(fread(file = file, sep = sep, nrows = big.seek, 
                    skip = records_read, data.table = FALSE), silent = TRUE)
        if(!is.data.frame(Table_df)){
            break
        }
        Table_df <- Table_df[,c(chr1_col, chr2_col, value_col)]
        colnames(Table_df) <- c("chr1", "chr2", "value")
        if(any(Table_df[,"chr2"] < Table_df[,"chr1"])){
            stop(paste("Provided chr2_col was not larger",
                "than chr1_col!","This file is not in",
                "upper.tri sparse format!", sep = " "))
        }
        for (i in seq_along(chrs)) {
            chr1 <- chrs[i]
            chr1_start <- starts[i]
            chr1_end <- ends[i]
            chr1_filter <- Table_df[,"chr1"] >= chr1_start & 
                                Table_df[,"chr1"] <= chr1_end
            if(any(!chr1_filter)){
                break
            }
            chr1_min <- min(which(chr1_filter)) + records_read
            chr1_max <- max(which(chr1_filter)) + records_read
            for(j in seq(from = i, to = length(chrs))){
                chr2 <- chrs[j]
                chr2_start <- starts[j]
                chr2_end <- ends[j]
                chr1_chr2_pair_name <- paste(chr1, chr2, sep = "_")
                chr2_filter <- Table_df[,"chr2"] >= chr2_start & 
                                    Table_df[,"chr2"] <= chr2_end
                chr1_and_chr2 <- (chr1_filter & chr2_filter)
                if(any(!chr1_and_chr2)){
                    break
                }
                chr1_chr2_min <- min(which(chr1_and_chr2))
                chr1_chr2_max <- max(which(chr1_and_chr2))
                if(chr1_chr2_pair_name %in% names(indexes_list)){
                    temp_df <- indexes_list[[chr1_chr2_pair_name]]
                    if(chr1_max > temp_df$chr1_end){
                        temp_df$chr1_end <- chr1_max
                    }
                    if(chr2_max > temp_df$chr2_end){
                        temp_df$chr2_end <- chr1_chr2_max
                    }
                }else{
                    temp_df <- data.frame(
                        chr1 = chr1,
                        chr1_start = chr1_min,
                        chr1_end = chr1_max,
                        chr2 = chr2,
                        chr2_start = chr1_chr2_min,
                        chr2_end = chr1_chr2_max)
                }
                indexes_list[[chr1_chr2_pair_name]] <- temp_df
            }
        }
        records_read <- records_read + nrow(Table_df)
    }
    indexes_df <- do.call(rbind, indexes_list)
    row.names(indexes_df) <- NULL
    return(indexes_df)
}


._ProcessTable_ <- function(Brick, table_handle, delim, 
    group_path, exec, num_rows = 1000000, is_sparse, 
    chr1_col, chr2_col, value_col, chr1_start, chr2_start, 
    chr1_end, chr2_end, compute_sparsity, distance,
    sep, sparsity_bins = 100){

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