.make_coords_list <- function(Brick, resolution){
    Reference.object <- GenomicMatrix$new()
    Brick_filelist <- BrickContainer_list_files(Brick = Brick, 
        resolution = resolution)
    filepath_split <- split(Brick_filelist$filepaths, 
        Brick_filelist$chrom1)
    chr2_split <- split(Brick_filelist$chrom2, 
        Brick_filelist$chrom1)
    Brick_filepath_list <- lapply(names(filepath_split), function(chr1){
        paths <- filepath_split[[chr1]]
        done_vector <- vapply(chr2_split[[chr1]], function(chr2){
            Brick_matrix_isdone(Brick = Brick, resolution = resolution,
                chr1 = chr1, chr2 = chr2)
        }, TRUE)

        chr2s <- chr2_split[[chr1]][done_vector]
        paths <- paths[done_vector]
        if(length(chr2s) == 0){
            return(NULL)
        }
        names(paths) <- chr2s
        return(paths)
    })
    names(Brick_filepath_list) <- names(filepath_split)
    Brick_filepath_list <- Brick_filepath_list[
        !vapply(Brick_filepath_list, is.null, TRUE)]
    if(length(Brick_filepath_list) == 0){
        stop("It looks like no files were loaded in this container")
    }


    Coords_df_list <- lapply(names(Brick_filepath_list), function(chr1){
        chr2s <- names(Brick_filepath_list[[chr1]])

        chr1_chr2_pairs_df_list <- lapply(chr2s, function(chr2){
            max_dist <- Brick_matrix_maxdist(Brick = Brick, 
                chr1 = chr1, chr2 = chr2, resolution = resolution)
            Dimensions <- Brick_matrix_dimensions(Brick = Brick, 
                chr1 = chr1, chr2 = chr2, resolution = resolution)
            col_length <- min(max_dist, Dimensions[2])
            chr1_row_sums <- Brick_get_matrix_mcols(
                Brick = Brick, 
                chr1 = chr1, 
                chr2 = chr2,
                resolution = resolution,
                what = Reference.object$hdf.matrix.rowSums)
            pairs_df <- .make_chr1_chr2_pairs(chr1 = chr1, chr2 = chr2, 
                chr1_start = 1, chr1_end = Dimensions[1], 
                chr1_row_sums = chr1_row_sums, chr2_row = 1, 
                chr2_length = col_length)
            pairs_df$filepath <- Brick_filepath_list[[chr1]][chr2]
            return(pairs_df)
        })

        chr1_chr2_pairs_df <- do.call(rbind, chr1_chr2_pairs_df_list)
        chr1_chr2_pairs_df <- arrange(chr1_chr2_pairs_df, "chr1_start")
        return(chr1_chr2_pairs_df)

    })

    Coords_df <- do.call(rbind, Coords_df_list)
    return(Coords_df)
}

.make_chr1_chr2_pairs <- function(chr1, chr2, chr1_start, chr1_end, 
    chr1_row_sums, chr2_row, chr2_length){
    Row_seq <- seq(from = chr1_start, to = chr1_end)
    chr2_diag_offset <- TRUE
    if(chr1 != chr2){
        chr2_diag_offset <- FALSE
    }
    Col_rep <- rep(chr2_row, times = length(Row_seq))
    length_rep <- rep(chr2_length, times = length(Row_seq))
    if(chr2_diag_offset){
        length_rep <- length_rep - (Row_seq - 1)
        Col_rep <- Row_seq
    }
    temp_df <- data.frame(chr1 = chr1, chr2 = chr2, chr1_start = Row_seq,
            chr2_start = Col_rep, chr2_length = length_rep, 
            stringsAsFactors = FALSE)
    if(any(chr1_row_sums == 0)){
        Which <- which(chr1_row_sums == 0)
        temp_df <- temp_df[!(temp_df$chr1_start %in% Which),]
    }
    return(temp_df)
}

.fetch_upper_tri_value_by_row <- function(Brick_filepath, chr1, chr2, 
    row, col, chr2_length){
    Reference.object <- GenomicMatrix$new()
    Group.path <- Create_Path(c(Reference.object$hdf.matrices.root,
        chr1, chr2))
    Start <- c(row, col)
    Stride <- c(1,1)
    Count <- c(1, chr2_length)
    Vector <- ._Brick_Get_Something_(Group.path = Group.path, 
        Brick = Brick_filepath, Name = Reference.object$hdf.matrix.name, 
        Start = Start, Stride = Stride, Count = Count, return.what = "data")
    return(Vector)
}