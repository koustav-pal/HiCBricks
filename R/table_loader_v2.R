add_to_data <- function(Vector = NULL, start = NULL, end = NULL,  
    data = NULL){ 
    Vector[start:end] <- Vector[start:end] + data 
    return(Vector)    
}

prepare_empty_metrics_list <- function(starts.1 = NULL, ends.1 = NULL,    
    starts.2 = NULL, ends.2 = NULL, chrom1 = NULL, chrom2 = NULL){    
    if(chrom1 != chrom2){ 
        row.sums <- list(rep(0,(max(ends.1) - min(starts.1))+1),
            rep(0,(max(ends.2) - min(starts.2))+1))
        names(row.sums) <- c(chrom1,chrom2)
        bin.coverage <- list(rep(0,(max(ends.1) - min(starts.1))+1),  
            rep(0,(max(ends.2) - min(starts.2))+1))
        names(bin.coverage) <- c(chrom1,chrom2)
        extent <- c(0,0)  
    }else{    
        row.sums <- rep(0,(max(ends.1) - min(starts.1))+1)    
        bin.coverage <- rep(0,(max(ends.1) - min(starts.1))+1)    
        extent <- c(0,0)  
        # sparsity.index <- list("values" = rep(0,(   
        # max(ends.1) - min(starts.1))+1),    
        #     "remaining.val" = rep(0,(max(ends.1) - min(starts.1))+1),   
        #     "remaining.val.coords" = NA)    
    } 
    A.list <- list("row.sums" = row.sums, "bin.coverage" = bin.coverage,  
        "extent" = extent) #"sparsity" = sparsity.index)  
    return(A.list)    
}

insert_data_and_computemetrics_both_matrices <- function(Brick = NULL,  
    Matrix = NULL, group.path = NULL, chrom1 = NULL, chrom2 = NULL,   
    row.offset = NULL, col.offset = NULL, row.pos = NULL, col.pos = NULL,
    metrics.list = NULL){ 
    Reference.object <- GenomicMatrix$new()   
    real.row.coords <- seq(1,nrow(Matrix),by = 1) + row.offset
    real.col.coords <- seq(1,ncol(Matrix),by = 1) + col.offset
    Values <- Matrix[cbind(row.pos,col.pos)]
    dataset.name <- as.character(Reference.object$hdf.matrix.name)    
    if(chrom1 == chrom2){ 
        if(all(real.col.coords %in% real.row.coords)){    
            Values <- Matrix[cbind(row.pos,col.pos)]  
            Matrix[cbind(col.pos,row.pos)] <- Values  
        }else{    
            Start <- c(min(real.col.coords), min(real.row.coords))    
            Stride <- c(1,1)  
            Count <- dim(t(Matrix))
            ._Brick_Put_Something_(Group.path = group.path, Brick = Brick,    
                Name = dataset.name, data = t(Matrix), Start = Start, 
                Stride = Stride, Count = Count)
            Matrix[is.na(Matrix) | is.infinite(Matrix) | is.nan(Matrix)] <- 0 
            metrics.list[["row.sums"]] <- add_to_data(    
                Vector = metrics.list[["row.sums"]],  
                start = min(real.col.coords), 
                end = max(real.col.coords),   
                data = rowSums(t(Matrix)))
            metrics.list[["bin.coverage"]] <-add_to_data( 
                Vector = metrics.list[["bin.coverage"]],  
                start = min(real.col.coords), 
                end = max(real.col.coords),   
                data = rowSums(t(Matrix) > 0))    
        } 
        Start <- c(min(real.row.coords), min(real.col.coords))
        Stride <- c(1,1)  
        Count <- dim(Matrix)
        ._Brick_Put_Something_(Group.path = group.path, Brick = Brick,    
            Name = dataset.name, data = Matrix, Start = Start,    
            Stride = Stride, Count = Count)   
        Matrix[is.na(Matrix) | is.infinite(Matrix) | is.nan(Matrix)] <- 0 
        metrics.list[["row.sums"]] <- add_to_data(    
            Vector = metrics.list[["row.sums"]],  
            start = min(real.row.coords), 
            end = max(real.row.coords),   
            data = rowSums(Matrix))   
        metrics.list[["bin.coverage"]] <-add_to_data( 
            Vector = metrics.list[["bin.coverage"]],  
            start = min(real.row.coords), 
            end = max(real.row.coords),   
            data = rowSums(Matrix > 0))   
    }else{    
        Start <- c(min(real.row.coords), min(real.col.coords))    
        Stride <- c(1,1)  
        Count <- dim(Matrix)
        ._Brick_Put_Something_(Group.path = group.path, Brick = Brick,    
            Name = dataset.name, data = Matrix, Start = Start,    
            Stride = Stride, Count = Count)
        Matrix[is.na(Matrix) | is.infinite(Matrix) | is.nan(Matrix)] <- 0 
        metrics.list[["row.sums"]][[chrom2]] <- add_to_data(
                Vector = metrics.list[["row.sums"]][[chrom2]],
                start = min(real.col.coords),
                end = max(real.col.coords),
                data = rowSums(t(Matrix)))
        metrics.list[["row.sums"]][[chrom1]] <- add_to_data(
                Vector = metrics.list[["row.sums"]][[chrom1]],
                start = min(real.row.coords),
                end = max(real.row.coords),
                data = rowSums(Matrix))
        metrics.list[["bin.coverage"]][[chrom2]] <- add_to_data(
                Vector = metrics.list[["bin.coverage"]][[chrom2]],
                start = min(real.col.coords),
                end = max(real.col.coords),
                data = rowSums(t(Matrix) > 0))    
        metrics.list[["bin.coverage"]][[chrom1]] <- add_to_data(
                Vector = metrics.list[["bin.coverage"]][[chrom1]],    
                start = min(real.row.coords),
                end = max(real.row.coords),
                data = rowSums(Matrix > 0))
    } 
    Min <- metrics.list[["extent"]][1]    
    Max <- metrics.list[["extent"]][2]    
    if(min(Matrix) < Min | Min == 0){ 
        metrics.list[["extent"]][1] <- min(Matrix)
    } 
    if(max(Matrix) > Max){    
        metrics.list[["extent"]][2] <- max(Matrix)
    } 
    return(metrics.list)  
}

.make_iterations <- function(Start.pos = NULL, End.pos = NULL,
    step = NULL){
    if((End.pos - Start.pos) < step){
        return(list(start = Start.pos, end = End.pos))    
    }
    Starts <- seq(from = Start.pos, to = End.pos, by = step)
    Starts <- Starts[Starts != End.pos]
    Ends <- c(Starts[-1] - 1, End.pos)
    return(list(start = Starts, end = Ends))
}

.return_table_df <- function(file, delim, num_rows, skip_rows, col_index, 
    chr1_extent, chr2_extent){
    Table_df <- try(fread(file = file, sep = delim, nrows = num_rows, 
                skip = skip_rows, data.table = FALSE), silent = TRUE)
    if(!is.data.frame(Table_df)){
        return(data.frame(chr1 = c(), chr2 = c(), value = c()))
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
    Test_table_df <- .return_table_df(file = file, delim = delim, 
        num_rows = 10, skip_rows = 0, col_index = col_index)
    # .check_column_classes(a_dataframe = Test_table_df)
    message("Indexing the file...")
    indexes_list <- list()
    records_read <- 0
    chr_start_shift <- 0
    # i = 1
    Table_df <- .return_table_df(file = file, delim = delim, 
        num_rows = big_seek, skip_rows = records_read, 
        col_index = col_index)
    while(nrow(Table_df) > 0) {
        .check_upper_triangle(a_dataframe = Table_df)
        for (i in seq_along(chrs)) {
            chr1 <- chrs[i]
            chr1_start <- starts[i]
            chr1_end <- ends[i]
            chr1_filter <- Table_df[,"chr1"] >= chr1_start & 
                Table_df[,"chr1"] <= chr1_end
            if(all(!chr1_filter)){
                next
            }
            # start_pos <- min(which(chr1_filter))
            Table_df_subset <- Table_df[chr1_filter,]
            chr1_run_lengths <- rle(Table_df_subset[,"chr1"])
            chr1_end_positions <- cumsum(chr1_run_lengths$lengths) + 
                chr_start_shift
            chr1_start_positions <- chr1_end_positions - 
                chr1_run_lengths$lengths + 1 
            if(chr1 %in% names(indexes_list)){
                temp_df <- indexes_list[[chr1]]
                temp_df_1 <- NULL
                temp_df_2 <- NULL
                Filter <- (chr1_run_lengths$values %in% temp_df$chr1_row)
                if(any(Filter)){
                    temp_df_1 <- .replace_chr1_row_indexes(
                        a_dataframe = temp_df,
                        end_position_vector = chr1_end_positions,
                        id_vector = chr1_run_lengths$values)
                }
                if(any(!Filter)){
                    temp_df_2 <- data.frame(
                        chr1 = chr1,
                        chr1_row = chr1_run_lengths$values[!Filter],
                        chr1_start = chr1_start_positions[!Filter],
                        chr1_end = chr1_end_positions[!Filter])
                }
                temp_df <- rbind(temp_df_1, temp_df_2)
            }else{
                temp_df <- data.frame(
                    chr1 = chr1,
                    chr1_row = chr1_run_lengths$values,
                    chr1_start = chr1_start_positions,
                    chr1_end = chr1_end_positions)
            }
            chr_start_shift <- chr_start_shift + nrow(Table_df_subset)
            indexes_list[[chr1]] <- temp_df
        }
        # message("Read ", chr_start_shift,"...")
        Table_df <- .return_table_df(file = file, delim = delim, 
            num_rows = big_seek, skip_rows = chr_start_shift, 
            col_index = col_index)
    }
    indexes_df <- do.call(rbind, indexes_list)
    row.names(indexes_df) <- NULL
    return(indexes_df)
}

# ==========================================================================
# Process a delimited file and load it into HDF files
# ==========================================================================
#
# Parameter definitions
# --------------------------------------------------------------------------
#
# Brick: A S4 object of class BrickContainer
# 
# table_file: A vector of class character specifying the path to a tsv 
# file to reac.
# 
# delim: The delimiter of the tsv file.
# 
# resolution: The resolution of the Hi-C matrix.
# 
# matrix_chunk: The chunk of the matrix to process.
# 
# batch_size: The number of lines to read from the tsv file per iteration.
# 
# col_index: A vector specifying the column index of interacting bins. The 
# index description is as follows: 
# - The originating bin of the interaction 
# - The destination bin of the interaction
# - The interaction signal itself
# 
# remove_prior: A binary vector of length 1 specifying if a previously loaded 
# matrix should be removed or not.
# 
# is_sparse: A binary vector of length 1 specifying if the matrix being loaded
# is sparse. Not to be confused with a sparse matrix format, sparse here means
# if it contains a lot of zeros. In this case a metric called the sparsity 
# index is computed which defines how many zeros are present per 100 bins from
# the diagonal of any given bin.
#
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
    Chrominfo_df <- Brick_get_chrominfo(Brick, resolution = resolution)
    end_positions <- cumsum(Chrominfo_df[,"nrow"])
    start_positions <- c(1, end_positions[-length(end_positions)] + 1)
    chr1_indexes_df <- .create_chr1_indexes(file = table_file, 
        delim = delim, 
        big_seek = batch_size, 
        chrs = Chrominfo_df$chr, 
        col_index = col_index, 
        starts = start_positions, 
        ends = end_positions)
    position_split_chr <- split(cbind(start_positions, end_positions), 
        Chrominfo_df[,"chr"])
    position_split_chr <- 
    chr_positions_list <- lapply(position_split_chr, function(x){
        iter_list <- .make_iterations(Start.pos = x[1], 
            End.pos = x[2], 
            step = matrix_chunk)
        iter_list
    })
    names(start_positions) <- Chrominfo_df$chr
    for (i in seq_len(nrow(Brick_files_tib))){
        chr1 <- Brick_files_tib$chrom1[i]
        chr1_starts <- chr_positions_list[[chr1]]$start
        chr1_ends <- chr_positions_list[[chr1]]$end
        chr1_widths <- cumsum(chr1_ends - chr1_starts + 1)
        chr1_hdf_offsets <- c(0, chr1_widths[-length(chr1_widths)])
        chr1_offset <- start_positions[chr1] - 1

        Indexes_chr1_filter <- chr1_indexes_df$chr1 == chr1

        chr2 <- Brick_files_tib$chrom2[i]
        chr2_starts <- chr_positions_list[[chr2]]$start
        chr2_ends <- chr_positions_list[[chr2]]$end
        chr2_widths <- cumsum(chr2_ends - chr2_starts + 1)
        chr2_hdf_offsets <- c(0, chr2_widths[-length(chr2_widths)])
        chr2_offset <- start_positions[chr2] - 1
        message(chr1, " vs ", chr2)
        message("chr1 starts: ", paste(chr1_starts, collapse = ", "))
        message("chr1 ends: ", paste(chr1_ends, collapse = ", "))
        message("chr2 starts: ", paste(chr2_starts, collapse = ", "))
        message("chr2 ends: ", paste(chr2_ends, collapse = ", "))
        metrics.list <- prepare_empty_metrics_list(
            starts.1 = chr1_starts,
            ends.1 = chr1_ends,
            starts.2 = chr2_starts,
            ends.2 = chr2_ends,
            chrom1 = chr1,
            chrom2 = chr2)

        Brick_filepath <- Brick_files_tib$filepaths[i]
        group_path <- Create_Path(c(
            Reference.object$hdf.matrices.root,
            chr1, chr2))

        chunk_pairs <- cbind(rep(seq_along(chr1_starts), 
                each = length(chr2_starts)),
            rep(seq_along(chr2_starts), 
                times = length(chr1_starts)))
        for(j in seq_len(nrow(chunk_pairs))) {
            chr1_start <- chr1_starts[chunk_pairs[j,1]]
            chr1_end <- chr1_ends[chunk_pairs[j,1]]
            chr2_start <- chr2_starts[chunk_pairs[j,2]]
            chr2_end <- chr2_ends[chunk_pairs[j,2]]
            chr1_hdf_offset <- chr1_hdf_offsets[chunk_pairs[j,1]]
            chr2_hdf_offset <- chr2_hdf_offsets[chunk_pairs[j,2]]

            chr1_rows_filter <- chr1_indexes_df$chr1_row >= chr1_start & 
                chr1_indexes_df$chr1_row <= chr1_end
            chr1_skip_rows <- min(chr1_indexes_df$chr1_start[
                            Indexes_chr1_filter & chr1_rows_filter]) - 1
            chr1_read_rows <- max(chr1_indexes_df$chr1_end[
                            Indexes_chr1_filter & 
                            chr1_rows_filter]) - chr1_skip_rows
            Table_df <- .return_table_df(
                file = table_file, 
                delim = delim, 
                num_rows = chr1_read_rows, 
                skip_rows = chr1_skip_rows, 
                col_index = col_index)

            Matrix <- matrix(data = 0, 
                nrow = (chr1_end - chr1_start) + 1,
                ncol = (chr2_end - chr2_start) + 1)
            Filter <- Table_df$chr2 >= chr2_start & Table_df$chr2 <= chr2_end
            if(!any(Filter)){
                next
            }
            Temp_table_df <- Table_df[Filter,]
            Temp_table_df$chr1 <- Temp_table_df$chr1 - 
                chr1_offset - ((chr1_start - chr1_offset) - 1)
            Temp_table_df$chr2 <- Temp_table_df$chr2 - 
                chr2_offset - ((chr2_start - chr2_offset) - 1)
                # message("Row segment: ", 
                #     paste(chr1, chr1_start - chr1_offset, 
                #         chr1_end - chr1_offset, sep = ":"), 
                #     "; Col segment: ", 
                #     paste(chr2, chr2_start - chr2_offset, 
                #         chr2_end - chr2_offset, sep = ":"))
                # message("Row range: ", min(Temp_table_df$chr1), " ", 
                #     max(Temp_table_df$chr1),"; Col range:", 
                #     min(Temp_table_df$chr2), " ", max(Temp_table_df$chr2))
            Matrix[cbind(Temp_table_df$chr1, 
                Temp_table_df$chr2)] <- Temp_table_df$value
            metrics.list <- insert_data_and_computemetrics_both_matrices(
                Brick = Brick_filepath,
                Matrix = Matrix,
                group.path = group_path,
                chrom1 = chr1,
                chrom2 = chr2,
                row.offset = chr1_hdf_offset,
                col.offset = chr2_hdf_offset,
                row.pos = Temp_table_df$chr1,
                col.pos = Temp_table_df$chr2,
                metrics.list = metrics.list)
        }
        distance <- max(chr2_ends) - chr2_offset - 1
        matrix_range <- metrics.list[["extent"]]
        Attributes <- Reference.object$matrices.chrom.attributes
        attr_vals <- c(basename(table_file),
            as.double(matrix_range),
            as.integer(is_sparse),
            as.integer(distance),
            as.integer(TRUE))

        if(is.list(metrics.list[["row.sums"]])){
            chr1_length <- length(metrics.list[["row.sums"]][[chr1]])
            chr2_length <- length(metrics.list[["row.sums"]][[chr2]])
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.rowSums,
                object = metrics.list[["row.sums"]][[chr1]])
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.coverage,
                object = metrics.list[["bin.coverage"]][[chr1]]/chr1_length)
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.colSums,
                object = metrics.list[["row.sums"]][[chr2]])
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.coverage.t,
                object = metrics.list[["bin.coverage"]][[chr2]]/chr2_length)
            WriteAttributes(Path = group_path, File = Brick_filepath,
                Attributes = Attributes, values = attr_vals, on = "group")
        }else{
            chr1_length <- length(metrics.list[["row.sums"]])
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.rowSums,
                object = metrics.list[["row.sums"]])
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.colSums,
                object = metrics.list[["row.sums"]])
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.coverage,
                object = metrics.list[["bin.coverage"]]/chr1_length)
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
                name = Reference.object$hdf.matrix.coverage.t,
                object = metrics.list[["bin.coverage"]]/chr1_length)
            WriteAttributes(Path = group_path, File = Brick_filepath,
                Attributes = Attributes, values = attr_vals, on = "group")
        }
    }
    return(TRUE)
}