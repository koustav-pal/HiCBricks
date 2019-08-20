.return_table_df <- function(file, delim, num_rows, skip_rows, col_index, 
    chr1_extent, chr2_extent){
    Table_df <- try(fread(file = file, sep = delim, nrows = num_rows, 
                skip = skip_rows, data.table = FALSE), silent = TRUE)
    if(!is.data.frame(Table_df)){
        return(NA)
    }
    Table_df <- Table_df[,col_index]
    colnames(Table_df) <- c("chr1", "chr2", "value")
    return(Table_df)
}

.replace_chr1_row_indexes <- function(a_dataframe, end_position_vector, 
    id_vector){
    Which <- which(a_dataframe$chr1_row %in% id_vector)
    a_dataframe$chr1_end[Which] <- vapply(Which, function(x){
        my_id <- a_dataframe$chr1_row[x]
        end_position_vector[id_vector == my_id]
    },1)
    return(a_dataframe)
}

.check_column_classes <- function(a_dataframe){
    if(!(ncol(a_dataframe) >= 3)){
        stop("A table in the upper triangle sparse format is",
            " expected with atleast 3 columns")
    }
    chr1_class <- class(a_dataframe[,1])
    chr2_class <- class(a_dataframe[,2])
    value_class <- class(a_dataframe[,3])
    if(!all(c("numeric", "numeric", "numeric") %in% 
            c(chr1_class, chr2_class, value_class))){
        stop("Expected numeric values for chr1, chr2 bins and values.")
    }    
}

.check_upper_triangle <- function(a_dataframe){
    if(any(a_dataframe[,"chr2"] < a_dataframe[,"chr1"])){
    stop(paste("Provided chr2_col was not larger",
        "than chr1_col!","This file is not in",
        "upper.tri sparse format!", sep = " "))
    }
}

.create_chr1_indexes <- function(file, delim, big_seek = 1000000,
    chrs, col_index, starts, ends){
    Test_table_df <- .return_table_df(file = file, delim = sep, num_rows = 10,
            skip_rows = 0, col_index = col_index)
    .check_column_classes(a_dataframe = Test_table_df)
    message("Indexing the file...")
    indexes_list <- list()
    records_read <- 0
    chr_start_shift <- 0
    i == 1
    while(i == 1) {
        Table_df <- .return_table_df(file = file, delim = sep, 
            num_rows = big_seek, skip_rows = records_read, 
            col_index = col_index)
        if(is.na(Table_df)){
            break
        }
        .check_upper_triangle(a_dataframe = Table_df)
        for (i in seq_along(chrs)) {
            chr1 <- chrs[i]
            chr1_start <- starts[i]
            chr1_end <- ends[i]
            chr1_filter <- Table_df[,"chr1"] >= chr1_start & 
                Table_df[,"chr1"] <= chr1_end
            start_pos <- min(which(chr1_filter))
            Table_df_subset <- Table_df[chr1_filter,]
            chr1_run_lengths <- rle(Table_df[,"chr1"])
            chr1_end_positions <- cumsum(chr1_run_lengths$lengths) + 
                chr_start_shift
            chr1_start_positions <- c(start_pos + chr_start_shift, 
                chr1_end_positions[length(chr1_end_positions)] - 1)
            if(all(!chr1_filter)){
                break
            }
            if(chr1 %in% names(indexes_list)){
                temp_df <- indexes_list[[chr1]]
                if(any(chr1_run_lengths$values %in% temp_df$chr1_row)){
                    temp_df <- .replace_chr1_row_indexes(
                        a_dataframe = temp_df,
                        end_position_vector = chr1_end_positions,
                        # start_position_vector = chr1_start_positions,
                        id_vector = chr1_run_lengths$values)
                }
            }else{
                temp_df <- data.frame(
                    chr1 = chr1,
                    chr1_row = chr1_run_lengths$values,
                    chr1_start = chr1_start_positions,
                    chr1_end = chr1_end_positions)
            }
            indexes_list[[chr1]] <- temp_df
        }
        records_read <- records_read + nrow(Table_df)
        chr_start_shift <- chr_start_shift + records_read
    }
    indexes_df <- do.call(rbind, indexes_list)
    row.names(indexes_df) <- NULL
    return(indexes_df)
}

.process_tsv <- function(Brick, table_file, delim, resolution, matrix_chunk, 
    batch_size, col_index, remove_prior = TRUE, is_sparse = FALSE, 
    sparsity_bins){
    Reference.object <- GenomicMatrix$new()
    if(is_sparse){
        sparsity_bins = sparsity_bins
    }
    if(length(col_index) != 3){
        stop("col_index must be of length 3, defining columns",
            " from_bin, to_bin, value.")
    }

    Brick_files_tib <- BrickContainer_list_files(Brick, 
        resolution = resolution)
    Chrominfo_df <- Brick_get_chrominfo(Brick)
    end_positions <- cumsum(Chrominfo_df[,"nrow"])
    start_positions <- c(1, end_positions[-length(end_positions)] + 1)
    position_split_chr <- split(cbind(start_positions, end_positions), 
        Chrominfo_df[,"chr"])
    chr_positions_list <- lapply(position_split_chr, function(x){
        make_mcool_iterations(Start.pos = x[1], End.pos = x[2], 
            step = matrix_chunk)
    })
    Indexes_df <- .create_read_indexes(file = table_file, delim = delim, 
        big_seek = 5000000, chrs = names(chr_positions_list), 
        col_index = col_index, starts = start_positions, ends = end_positions)


    names(start_positions) <- Chrominfo_df$chr

    skip_rows <- 0
    for (i in seq_along(nrow(Indexes_df))) {

        current_index_df <- Indexes_df[i,]
        rows_to_read <- current_index_df$chr2_end - 
            current_index_df$chr2_start + 1
        Table_df <- fread(file = file, sep = sep, nrows = rows_to_read, 
                    skip = skip_rows, data.table = FALSE)
        colnames(Table_df) <- c("from_bin", "to_bin", "value")
        skip_rows <- skip_rows + rows_to_read

        chr1_positions <- chr_positions_list[[current_index_df$chr1]]
        chr1_offset <- start_positions[current_index_df$chr1] - 1
        chr1_starts <- chr1_positions$start
        chr1_ends <- chr1_positions$end
        chr2_positions <- chr_positions_list[[current_index_df$chr2]]
        chr2_offset <- start_positions[current_index_df$chr2] - 1
        chr2_starts <- chr2_positions$start
        chr2_ends <- chr2_positions$end

        k <- 0
        chr1_previous_starts <- 0
        for (i in seq_along(chr1_starts)) {
            chr1_start <- chr1_starts[i] - chr1_offset - chr1_previous_starts
            chr1_end <- chr1_ends[i] - chr1_offset - chr1_previous_starts
            chr2_previous_starts <- 0
            for (j in seq_along(chr2_starts)) {
                Table_df_subset <- Table_df[
                    Table_df$from_bin >= chr1_starts[i] & 
                    Table_df$from_bin <= chr1_end[i] & 
                    Table_df$to_bin >= chr2_starts[j] & 
                    Table_df$to_bin <= chr2_end[j]]
                chr2_start <- chr2_starts[j] - chr2_offset - 
                    chr2_previous_starts
                chr2_end <- chr2_ends[j] - chr2_offset - 
                    chr2_previous_starts
                Table_df_subset$from_bin <- Table_df_subset$from_bin - 
                    chr1_offset - chr1_previous_starts
                Table_df_subset$to_bin <- Table_df_subset$to_bin - 
                    chr2_offset - chr2_previous_starts
                Matrix <- matrix(0, nrow = (chr1_end - chr1_start + 1),
                    ncol = (chr2_end - chr2_start + 1))
                Matrix[cbind(Table_df_subset$from_bin, 
                        Table_df_subset$to_bin)] <- Table_df_subset$value
                Matrix[is.na(Matrix) | 
                    is.infinite(Matrix) | 
                    is.nan(Matrix)] <- 0
                Start <- c(chr1_starts[i] - chr1_offset, 
                    chr2_ends[j] - chr2_offset)
                Stride <- c(1,1)
                Count <- dim(Matrix)
                group_path <- c(Reference.object$hdf.matrices.root, 
                    current_index_df$chr1, 
                    current_index_df$chr2)
                ._Brick_Put_Something_(
                    Group.path = group_path, 
                    Brick = Brick,
                    Name = Reference.object$hdf.matrix.name, 
                    data = Matrix, 
                    Start = Start,
                    Stride = Stride, 
                    Count = Count)
                ._Brick_Put_Something_(
                    Group.path = group_path, 
                    Brick = Brick,
                    Name = Reference.object$hdf.matrix.name, 
                    data = t(Matrix), 
                    Start = rev(Start),
                    Stride = Stride, 
                    Count = rev(Count))
                chr2_previous_starts <- chr2_previous_starts + chr2_start[i]
            }
            chr1_previous_starts <- chr1_previous_starts + chr1_starts[i]
        }
    }

    # k <- 1
    # skip_rows <- 0
    # while(i == 1) {
    #     Table_df <- try(fread(file = table_handle, sep = delim, 
    #         skip = skip_rows, nrows = batch_size, data.table = FALSE))
    #     skip_rows <- skip_rows + nrow(Table_df)
    #     skip_next <- 
    #     Table_df <- Table_df[,col_index]
    #     colnames(Table_df) <- c("from_bin", "to_bin", "value")
    #     chr1_positions <- position_split_chr[sapply(position_split_chr, 
    #         function(x){ any(x[2] >= Table_df$from_bin & 
    #             x[1] <= Table_df$from_bin)
    #     })]
    #     chr2_positions <- position_split_chr[sapply(position_split_chr, 
    #         function(x){ any(x[2] >= Table_df$to_bin & 
    #             x[1] <= Table_df$to_bin)
    #     })]
    #     for (chr1 in names(chr1_positions)) {
    #         chr1_position_list <- chr_positions_list[[chr1]]
    #         for (chr2 in names(chr1_positions)) {
    #             chr2_position_list <- chr_positions_list[[chr2]]
    #             chr1
    #         }
    #     }

    #     k <- k + 1
    # }


    # Metrics.list <- ._Compute_various_matrix_metrics(Matrix = Matrix, 
    #     compute.sparsity = compute.sparsity, sparsity.bins = sparsity_bins, 
    #     range = Matrix.range, distance = distance, 
    #     diag.position.start = Skip + 1)
    #     Matrix.range <- Metrics.list[["extent"]]
    #     Bin.coverage <- c(Bin.coverage,Metrics.list[["bin.cov"]])
    #     Row.sums <- c(Row.sums,Metrics.list[["row.sum"]])
    #     Sparsity.Index <- c(Sparsity.Index,Metrics.list[["sparsity"]])

    #     Cumulative.data <- rbind(Cumulative.data,Matrix)
    #     Obj.size <- object.size(Cumulative.data)
    #     if(Obj.size >= Reference.object$Max.vector.size | 
    #         i == length(Iterations)){
    #         Start <- c(Start.row,1)
    #         Stride <- c(1,1)
    #         Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
    #         message("Inserting Data at location: ",Start[1])
    #         message("Data length: ",Count[1])
    #         ._Brick_Put_Something_(Group.path=Group.path, Brick = Brick, 
    #             Name = Reference.object$hdf.matrix.name,
    #             data = Cumulative.data, Start = Start, Stride = Stride, 
    #             Count = Count)
    #         Start.row <- Start.row + Count[1]
    #         Set.col <- TRUE
    #         Cumulative.data <- NULL
    #         message("Loaded ",Obj.size," bytes of data...")
    #     }
    #     message("Read ",(Skip+Iter)," records...")
    #     i<-i+1

    # ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
    #     name = Reference.object$hdf.matrix.rowSums, object = Row.sums)
    # ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
    #     name = Reference.object$hdf.matrix.coverage, object = Bin.coverage)
    # if(compute.sparsity){
    #     ._Brick_WriteArray_(Brick = Brick, Path = Group.path, 
    #         name = Reference.object$hdf.matrix.sparsity, 
    #         object = Sparsity.Index)
    # }
    # Attributes <- Reference.object$matrices.chrom.attributes
    # Attr.vals <- c(basename(Matrix.file),
    #     as.double(Matrix.range),
    #     as.integer(is.sparse),
    #     as.integer(distance),
    #     as.integer(TRUE))
    # WriteAttributes(Path = Group.path, File = Brick, 
    #     Attributes = Attributes, 
    #     values = Attr.vals, 
    #     on = "group")
    # return(TRUE)
}


Brick_load_data_from_table <- function(Brick, table_file, delim = " ", 
    resolution = NULL, batch_size = 1000000, matrix_chunk = 2000, 
    col_index = c(1, 2, 3), remove_prior = FALSE, is_sparse = FALSE, 
    sparsity_bins = 100) {
    Reference.object <- GenomicMatrix$new()

    BrickContainer_class_check(Brick)
    Resolutions <- BrickContainer_list_resolutions(Brick)
    resolution <- .format_resolution(resolution)

    if(!(resolution %in% Resolutions)) {
        stop("resolution does not exist in the BrickContainer")
    }

    RetVar <- .process_tsv(Brick = Brick, table_file = table_file, 
        delim = delim, col_index = col_index, remove_prior = remove_prior,
        is_sparse = is_sparse, resolution = resolution, 
        sparsity_bins = sparsity_bins, matrix_chunk = matrix_chunk, 
        batch_size = batch_size, has_header = has_header)
}

#
# TODO
#
# => rewrite this part 
