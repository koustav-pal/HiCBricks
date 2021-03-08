.create_matrix_group_path <- function(has_resolution = FALSE, 
    mcool_version, resolution){
    Reference_object <- GenomicMatrix$new()
    Keys <- Reference_object$mcool.matrix.keys(version = mcool_version)
    if(has_resolution){
        group_path <- Create_Path(
            c(Reference_object$mcool.resolutions.name, 
                resolution,
                Keys[1]))
    }else{
        group_path <- Create_Path(c(Keys[1]))
    }
    return(group_path)
}

.create_bintable_group_path <- function(has_resolution = FALSE, 
    mcool_version, resolution){
    Reference_object <- GenomicMatrix$new()
    Keys <- Reference_object$mcool.bintable.keys(
        version = mcool_version)
    if(has_resolution){
        group_path <- Create_Path(
            c(Reference_object$mcool.resolutions.name, 
                resolution,
                Keys[1]))
    }else{
        group_path <- Create_Path(c(Keys[1]))
    }
    return(group_path)
}

.create_indexes_group_path <- function(has_resolution = FALSE, 
    mcool_version, resolution){
    Reference_object <- GenomicMatrix$new()
    Keys <- Reference_object$mcool.index.keys()
    indexes_group <- Keys[1]
    if(has_resolution){
        group_path <- Create_Path(
            c(Reference_object$mcool.resolutions.name, 
                resolution,
                indexes_group[1]))
    }else{
        group_path <- Create_Path(c(indexes_group[1]))
    }
    return(group_path)
}

.return_bias_vector <- function(mcool_path, mcool_version, has_resolution, 
    resolution, norm_factor = NULL){
    bintable_group_path <- .create_bintable_group_path(
        has_resolution = has_resolution, 
        mcool_version = mcool_version, 
        resolution = resolution)
    bias_vector <- NULL
    if(!is.null(norm_factor)){
        bias_vector <- ._Brick_Get_Something_(
            Group.path = bintable_group_path,
            Brick = mcool_path, 
            Name = norm_factor, 
            return.what = "data")
    }
    return(bias_vector)
}

.return_offsets <- function(mcool_path, mcool_version, has_resolution, 
    resolution){
    Reference_object <- GenomicMatrix$new()
    indexes_keys <- Reference_object$mcool.index.keys()
    indexes_bin <- indexes_keys[2]
    indexes_chrom <- indexes_keys[3]
    index_group_path <- .create_indexes_group_path(
        has_resolution = has_resolution, 
        mcool_version = mcool_version, 
        resolution = resolution)
    chrom_offset <- ._Brick_Get_Something_(
        Group.path = index_group_path,
        Brick = mcool_path, 
        Name = indexes_chrom, 
        return.what = "data")
    bin_offset <- ._Brick_Get_Something_(
        Group.path = index_group_path,
        Brick = mcool_path, 
        Name = indexes_bin, 
        return.what = "data")
    return(list("chrom_offset" = chrom_offset,
            "bin_offset" = bin_offset))
}

.return_chromosomes_list <- function(mcool_path, mcool_version, 
    has_resolution, resolution){
    cooler_remap_chrom <- ._mcool_remap_chromosomes(File = mcool_path,
        resolution = has_resolution, 
        binsize = resolution,
        mcool.version = mcool_version)
    return(cooler_remap_chrom[,"chr.name"])
}

.return_values_list <- function(mcool_path, mcool_version, read_from, read_to,
    has_resolution, resolution){
    Reference_object <- GenomicMatrix$new()
    matrix_keys <- Reference_object$mcool.matrix.keys(version = mcool_version)
    matrix_group_path <- .create_matrix_group_path(resolution = resolution,
        has_resolution = has_resolution, mcool_version = mcool_version)
    Bin1_id <- .fetch_data(file = mcool_path, group = matrix_group_path,
        name = matrix_keys[2], start = read_from, num_rows = read_to)
    Bin2_id <- .fetch_data(file = mcool_path, group = matrix_group_path,
        name = matrix_keys[3], start = read_from, num_rows = read_to)
    Counts <- .fetch_data(file = mcool_path, group = matrix_group_path,
        name = matrix_keys[4], start = read_from, num_rows = read_to)
    return(list("bin1_id" = Bin1_id, "bin2_id" = Bin2_id, "counts" = Counts))
}

.create_mcool_chr1_indexes <- function(chromosomes, ignore_chrs, matrix_chunk, 
    chr_offsets, bin_offsets, num_records_limit = 10000000){
    chr1_index_iterations_list <- .convert_chr_to_chunks(
        chr_offsets = chr_offsets, chromosomes = chromosomes,
        ignore_chrs = ignore_chrs, matrix_chunk = matrix_chunk)
    chromosomes <- chromosomes[!(seq_along(chromosomes) %in% ignore_chrs)]
    chr1_indexes_df_list <- lapply(seq_along(chromosomes), function(x){
        chr1_index_iterations <- chr1_index_iterations_list[[x]]
        records_read_per_iter <- chr1_index_iterations$read_to - 
            chr1_index_iterations$read_from
        if(any(records_read_per_iter > num_records_limit)){
            chr1_index_iterations <- .make_chunks_from_bin_offset(
                bin_offset = bin_offsets,
                num_records_limit = num_records_limit,
                chr1 = chromosomes[x],
                chr1_start = chr1_index_iterations$start[1],
                chr1_end = chr1_index_iterations$end[
                    length(chr1_index_iterations$end)])
        }
        names(chr1_index_iterations$end) <- NULL
        Bin_iter_df_list <- lapply(seq_along(chr1_index_iterations$start), 
            function(y){
                chr1_start <- chr1_index_iterations$start[y]
                chr1_end <- chr1_index_iterations$end[y]
                data.frame(
                    chr1 = chromosomes[x],
                    chr1_start = chr1_start,
                    chr1_end =  chr1_end,
                    read_from = bin_offsets[chr1_start] + 1,
                    read_to = bin_offsets[chr1_end+1],
                    stringsAsFactors = FALSE)
        })
        Bin_iter_df_list <- Bin_iter_df_list[!vapply(Bin_iter_df_list, 
            is.null, TRUE)]
        Bin_iter_df <- do.call(rbind, Bin_iter_df_list)
    })
    chr1_indexes_df <- do.call(rbind, chr1_indexes_df_list)
    chr1_indexes_df <- chr1_indexes_df[(chr1_indexes_df$read_to - 
        chr1_indexes_df$read_from) >= 0,]
    rownames(chr1_indexes_df) <- NULL
    return(chr1_indexes_df)
}

.convert_chr_to_chunks <- function(chr_offsets, chromosomes, ignore_chrs, 
    matrix_chunk){
    names(chr_offsets) <- NULL
    chunk_list <- lapply(seq_along(chromosomes), function(x){
        if(x %in% ignore_chrs){
            return(NULL)
        }
        .make_mcool_iterations(Start.pos = chr_offsets[x]+1,
        End.pos = chr_offsets[x+1],
        step = matrix_chunk)
    })
    return(chunk_list[!vapply(chunk_list, is.null, TRUE)])
}

.make_chunks_from_bin_offset <- function(bin_offset, num_records_limit, 
    chr1, chr1_start, chr1_end){
    chr1_seq <- seq(from = chr1_start, to = chr1_end+1, by = 1)
    num_records_per_bin <- diff(bin_offset[chr1_seq])
    chr1_seq <- chr1_seq[-length(chr1_seq)]
    Total_sums <- cumsum(num_records_per_bin)
    Divide_to_bins <- floor(Total_sums/num_records_limit) + 1
    chr1_chunk_list <- split(chr1_seq, Divide_to_bins)
    chr1_indexes_df_list <- lapply(chr1_chunk_list, function(chr1_chunk){
                data.frame(
                    chr1 = chr1,
                    chr1_start = min(chr1_chunk),
                    chr1_end =  max(chr1_chunk),
                    read_from = bin_offset[min(chr1_chunk)] + 1,
                    read_to = bin_offset[max(chr1_chunk)+1],
                    stringsAsFactors = FALSE)
    })
    chr1_indexes_df <- do.call(rbind, chr1_indexes_df_list)
    return(chr1_indexes_df)
}

.return_chr1_chr2_pairs <- function(Brick, chromosomes, matrix_chunk, 
    chr_offsets, ignore_chrs, bin_offsets, resolution, num_records_limit){
    chr1_bin_indexes_df <- .create_mcool_chr1_indexes(
        chromosomes = chromosomes, ignore_chrs = ignore_chrs,
        matrix_chunk = matrix_chunk, chr_offsets = chr_offsets, 
        bin_offsets = bin_offsets, num_records_limit = num_records_limit)
    chromosome_chunks_list <- .convert_chr_to_chunks(chr_offsets = chr_offsets,
        chromosomes = chromosomes, ignore_chrs = ignore_chrs, 
        matrix_chunk = matrix_chunk)
    chromosomes <- chromosomes[!(seq_along(chromosomes) %in% ignore_chrs)]
    names(chromosome_chunks_list) <- chromosomes
    Files_df <- BrickContainer_list_files(Brick = Brick, chr1 = chromosomes,
        chr2 = chromosomes, resolution = resolution)
    chr1_chr2_pairs <- split(Files_df$chrom2, Files_df$chrom1)
    chr1_chr2_rows_list <- lapply(seq_len(nrow(chr1_bin_indexes_df)), 
        function(x){
        chr1_rows <- chr1_bin_indexes_df[x,]
        chr1 <- chr1_rows$chr1
        chr2s <- chr1_chr2_pairs[[chr1]]
        chr2_rows_list <- lapply(chr2s, function(chr2){
            chr2_chunks <- chromosome_chunks_list[[chr2]]
            chr2_rows <- data.frame(chr1 = chr1, chr2 = chr2, 
                chr1_start = chr1_rows$chr1_start, 
                chr1_end = chr1_rows$chr1_end, chr2_start = chr2_chunks$start,
                chr2_end = chr2_chunks$end, read_from = chr1_rows$read_from,
                read_to = chr1_rows$read_to, filepath = Files_df$filepaths[
                Files_df$chrom1 == chr1 & Files_df$chrom2 == chr2],
                stringsAsFactors = FALSE)
            chr2_rows
        })
        chr2_rows_df <- do.call(rbind, chr2_rows_list)
    })
    chr1_chr2_rows_df <- do.call(rbind, chr1_chr2_rows_list)
    chr1_chr2_rows_df_list <- split(chr1_chr2_rows_df, 
        paste(chr1_chr2_rows_df$chr1, 
            chr1_chr2_rows_df$chr1_start, 
            chr1_chr2_rows_df$chr1_end, sep = ":"))
    return(chr1_chr2_rows_df_list)
}

.add_to_hdf_array <- function(Brick = NULL, Path = NULL, name = NULL, 
    object = NULL, Start, Stride, Count){
    if(!(length(c(Brick,Path,name,object))>=4)){
        stop("All arguments are required!")
    }
    Brick.handler <- ._Brick_Get_Something_(Group.path = Path, Brick = Brick, 
        Name = name, return.what = "group_handle")
    h5writeDataset.array(h5loc = Brick.handler, obj = object, name = name, 
        start = Start, stride = Stride, count = Count)
    H5Gclose(Brick.handler)
}

.populate_with_empty_metrics <- function(Brick_filepath, chr1, chr2, 
    row_length, col_length){
    Reference_object <- GenomicMatrix$new()
    chr1_rows <- rep(0, times = row_length)
    chr2_rows <- rep(0, times = col_length)
    group_path <- Create_Path(c(Reference_object$hdf.matrices.root, 
        chr1, chr2))
    ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.rowSums, object = chr1_rows)
    ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.colSums, object = chr2_rows)
    ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.coverage, object = chr1_rows)
    ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.coverage.t, object = chr2_rows)
}

.final_metrics_processing <- function(Brick, metrics_df, resolution){
    Reference_object <- GenomicMatrix$new()
    Attributes <- Reference_object$matrices.chrom.attributes
    chr1_df_list <- split(metrics_df, metrics_df$chr1)
    chr1_extent_df_list <- lapply(names(chr1_df_list), function(chr1){
        chr1_df <- chr1_df_list[[chr1]]
        chr2_df_split <- split(chr1_df, chr1_df$chr2)
        chr2_extent_df_list <- lapply(names(chr2_df_split), function(chr2){
            chr2_df <- chr2_df_split[[chr2]]
            Brick_filepath <- BrickContainer_get_path_to_file(Brick = Brick, 
                chr1 = chr1, chr2 = chr2, resolution = resolution)
            Min <- min(chr2_df$min)
            Max <- max(chr2_df$max)
            mcool_name <- unique(chr2_df$mcool_name)
            Distance <- max(chr2_df$distance)
            group_path <- Create_Path(c(Reference_object$hdf.matrices.root, 
                chr1, chr2))
            chr1_bin_coverage <- .fetch_data(file = Brick_filepath, start = 1, 
                group = group_path, num_rows = max(chr2_df$nrow),
                name = Reference_object$hdf.matrix.coverage)
            chr1_bin_coverage <- chr1_bin_coverage/max(chr2_df$nrow)
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
            name = Reference_object$hdf.matrix.coverage, 
            object = chr1_bin_coverage)

            chr2_bin_coverage <- .fetch_data(file = Brick_filepath, start = 1,
                group = group_path, num_rows = max(chr2_df$ncol),
                name = Reference_object$hdf.matrix.coverage.t)
            chr2_bin_coverage <- chr2_bin_coverage/max(chr2_df$ncol)
            ._Brick_WriteArray_(Brick = Brick_filepath, Path = group_path,
            name = Reference_object$hdf.matrix.coverage.t, 
            object = chr2_bin_coverage)

            Attr.vals <- c(mcool_name, as.double(Min), as.double(Max),
                as.integer(FALSE), as.integer(Distance), as.integer(TRUE))
            WriteAttributes(Path = group_path, File = Brick_filepath,
            Attributes = Attributes, values = Attr.vals, on = "group")
            data.frame(chr1 = chr1, chr2 = chr2, min = Min, max = Max)
        })
    do.call(rbind, chr2_extent_df_list)
    })
    chr1_extent_df <- do.call(rbind, chr1_extent_df_list)
    return(chr1_extent_df)
}

.add_metrics <- function(Brick_filepath, group_path, a_matrix, 
    dimensions, chr1_start, chr2_start){
    Reference_object <- GenomicMatrix$new()
    matrix_mcols_keys <- Reference_object$hdf.matrix.meta.cols()
    a_matrix[is.na(a_matrix) | is.nan(a_matrix) | is.infinite(a_matrix)] <- 0
    Row_sums <- rowSums(a_matrix)
    Col_sums <- colSums(a_matrix)
    chr1_bin_coverage_submat <- rowSums(a_matrix > 0)
    chr2_bin_coverage_submat <- colSums(a_matrix > 0)

    chr1_positions <- seq(from = chr1_start, by = 1, 
        length.out = nrow(a_matrix))
    chr2_positions <- seq(from = chr2_start, by = 1, 
        length.out = ncol(a_matrix))

    chr1_row_sums <- .fetch_data(file = Brick_filepath, group = group_path, 
        name = Reference_object$hdf.matrix.rowSums, start = chr1_start, 
        num_rows = dimensions[1])
    chr1_row_sums <- chr1_row_sums + Row_sums
    .add_to_hdf_array(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.rowSums, object = chr1_row_sums,
        Start = chr1_start, Stride = 1, Count = dimensions[1])

    chr2_col_sums <- .fetch_data(file = Brick_filepath, group = group_path, 
        name = Reference_object$hdf.matrix.colSums, start = chr2_start, 
        num_rows = dimensions[2])
    chr2_col_sums <- chr2_col_sums + Col_sums
    .add_to_hdf_array(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.colSums, object = chr2_col_sums,
        Start = chr2_start, Stride = 1, Count = dimensions[2])

    chr1_bin_coverage <- .fetch_data(file = Brick_filepath, group = group_path, 
        name = Reference_object$hdf.matrix.coverage, start = chr1_start, 
        num_rows = dimensions[1])
    chr1_bin_coverage <- chr1_bin_coverage + chr1_bin_coverage_submat
    .add_to_hdf_array(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.coverage, 
        object = chr1_bin_coverage, Start = chr1_start, Stride = 1, 
        Count = dimensions[1])

    chr2_bin_coverage <- .fetch_data(file = Brick_filepath, group = group_path, 
        name = Reference_object$hdf.matrix.coverage.t, start = chr2_start, 
        num_rows = dimensions[2])
    chr2_bin_coverage <- chr2_bin_coverage + chr2_bin_coverage_submat
    .add_to_hdf_array(Brick = Brick_filepath, Path = group_path,
        name = Reference_object$hdf.matrix.coverage.t,
        object = chr2_bin_coverage, Start = chr2_start, Stride = 1, 
        Count = dimensions[2])
    return(range(a_matrix))
}

.fetch_data <- function(file,  group, name, start, num_rows){
    a_vector <- ._Brick_Get_Something_(Group.path = group, Brick = file, 
        Start = start, Stride = 1, Count = num_rows, Name = name, 
        return.what = "data")
    return(a_vector)
}

.make_mcool_iterations <- function(Start.pos = NULL, End.pos = NULL,
    step = NULL){
    Starts <- seq(from = Start.pos, to = End.pos, by = step)
    Starts <- Starts[Starts != End.pos]
    Ends <- c(Starts[-1] - 1, End.pos)
    return(list(start = Starts, end = Ends))
}

._mcool_remap_chromosomes = function(File = NULL, mcool.version = NULL,
    resolution = FALSE, binsize = NULL){
    Reference_object <- GenomicMatrix$new()
    Bintable.keys <- Reference_object$mcool.bintable.keys(
        version = mcool.version)
    Scaffold.keys <- Reference_object$mcool.scaffold.keys(
        version = mcool.version)
    Bintable.group <- Bintable.keys[1]
    Bintable.chr <- Bintable.keys[2]
    Scaffold.group <- Scaffold.keys[1]
    Scaffold.name <- Scaffold.keys[3]
    if(resolution){ 
        Bintable.group.path <- Create_Path(
            c(Reference_object$mcool.resolutions.name,
                binsize,
                Bintable.group))
        Scaffold.group.path <-  Create_Path(
            c(Reference_object$mcool.resolutions.name,
                binsize,
                Scaffold.group))
    }else{
        Bintable.group.path <- Create_Path(Bintable.group)
        Scaffold.group.path <-  Create_Path(Scaffold.group)
    }
    Chrom.names <- ._Brick_Get_Something_(Group.path = Scaffold.group.path,
        Brick = File, Name = Scaffold.name, return.what = "data")
    Chrom.lengths <- ._Brick_Get_Something_(Group.path = Scaffold.group.path,
        Brick = File, Name = Scaffold.keys[2], return.what = "data")
    Chrom.bintable.ends <- ._Brick_Get_Something_(
        Group.path = Bintable.group.path, Brick = File,
        Name = Bintable.keys[4], return.what = "data")
    if(!all(Chrom.lengths %in% Chrom.bintable.ends)){
        stop("Unable to remap chromosomes! Please contact the developer!\n")
    }
    Which.ends <- match(Chrom.lengths,Chrom.bintable.ends)
    Offset <- c(Which.ends[1],
        (Which.ends[-1] - Which.ends[-length(Which.ends)]))
    Order.of.chroms <- order(Which.ends)
    Cooler.chrm.df <- data.frame(remap.chrom = Order.of.chroms,
        chr.name = Chrom.names[Order.of.chroms],
        lengths = Chrom.lengths[Order.of.chroms],
        offset = Offset)
    return(Cooler.chrm.df)
}

._mcool_bintable_ranges = function(mcool.file = NULL, mcool.remap.chrom = NULL,
    mcool.version = NULL, resolution = FALSE, binsize = NULL){
    Reference_object <- GenomicMatrix$new()
    Bintable.keys <- Reference_object$mcool.bintable.keys(
        version = mcool.version)
    Bintable.group <- Bintable.keys[1]
    Bintable.chr <- Bintable.keys[2]
    Bintable.start <- Bintable.keys[3]
    Bintable.end <- Bintable.keys[4]
    if(resolution){
        Bintable.group.path <- Create_Path(
            c(Reference_object$mcool.resolutions.name,binsize,Bintable.group)) 
    }else{
        Bintable.group.path <- Create_Path(Bintable.group)
    }
    Chrom <- do.call(c,lapply(seq(1,nrow(mcool.remap.chrom)),function(x){
            chr <- mcool.remap.chrom[x,"chr.name"]
            Offset <- mcool.remap.chrom[x,"offset"]
            rep(chr, Offset)
        }))
    Start <- ._Brick_Get_Something_(Group.path = Bintable.group.path,
        Brick = mcool.file,
        Name = Bintable.start, return.what = "data")
    End <- ._Brick_Get_Something_(Group.path = Bintable.group.path,
        Brick = mcool.file,
        Name = Bintable.end, return.what = "data")       
    if(length(unique(length(Chrom),length(Start),length(End)))==1){
        message("All ok! Chrom, Start, End have matching lengths...\n")
    }else{
        stop("All is not ok! Chrom, Start, End don't have matching lengths...",
            c(length(Chrom),length(Start),length(End)),"\n")
    }
    if(any(End %in% Start)){
        Start <- Start + 1
    }
    Cooler.ranges.df <- data.frame(chr = Chrom, start = as.numeric(Start),
        end = as.numeric(End), stringsAsFactors = FALSE)
    return(Cooler.ranges.df)
}

mcool_list_resolutions <- function(mcool = NULL){
    Reference_object <- GenomicMatrix$new()
    Handler <- ReturnH5FileConnection(File = mcool)
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    CloseH5Con(Handle = Handler, type = "file")
    if(Reference_object$mcool.resolutions.name %in% GroupList){
        Handler <- ._Brick_Get_Something_(
            Group.path = Create_Path(Reference_object$mcool.resolutions.name),
            Brick = mcool, return.what = "group_handle")
        BinList <- h5ls(Handler,
            datasetinfo = FALSE, recursive = FALSE)[,"name"]
        CloseH5Con(Handle = Handler, type = "group")
    }else{
        BinList <- NULL
    }
    return(.format_resolution(BinList))
}

.process_mcool <- function(Brick = NULL, mcool_path = NULL, resolution = NULL, 
    has_resolution = FALSE, matrix_chunk = 2000, norm_factor = NULL,
    cooler_read_limit = 10000000, keep_chr1 = NA, keep_chr2 = NA) {
    is.sparse <- FALSE
    Reference_object <- GenomicMatrix$new()
    mcool_version <- GetAttributes(Path = NULL, File = mcool_path,
        Attributes="format-version", on = "file",
        ignore.fun.cast = TRUE)[,"format-version"]
    mcool_version <- as.numeric(as.character(mcool_version))

    chrom_names <- .return_chromosomes_list(mcool_path = mcool_path, 
        mcool_version = mcool_version, has_resolution = has_resolution, 
        resolution = resolution)
    offset_list <- .return_offsets(mcool_path = mcool_path, 
        mcool_version = mcool_version, has_resolution = has_resolution,
        resolution = resolution)
    bias.vector <- .return_bias_vector(mcool_path = mcool_path, 
        mcool_version = mcool_version, has_resolution = has_resolution,
        resolution = resolution, norm_factor = norm_factor)

    index_chrom_offset <- offset_list[["chrom_offset"]]
    names(index_chrom_offset) <- chrom_names
    chrom_index_differences <- diff(index_chrom_offset) - 1
    ignore_chrs <- NULL
    if(any(chrom_index_differences == 0)){
        ignore_chrs <- which(chrom_index_differences == 0)
    }
    index_bin_offset <- offset_list[["bin_offset"]]
    message("Making read indices...")
    chr1_chr2_rows_df_list <- .return_chr1_chr2_pairs(Brick = Brick, 
        chromosomes = chrom_names, ignore_chrs = ignore_chrs, 
        matrix_chunk = matrix_chunk, num_records_limit = cooler_read_limit, 
        chr_offsets = index_chrom_offset, bin_offsets = index_bin_offset, 
        resolution = resolution)
    chr1_chr2_df <- do.call(rbind, lapply(chr1_chr2_rows_df_list, 
        function(a_row){
            a_row[,c("chr1", "chr2", "filepath")] 
        }))
    chr1_chr2_df <- unique(chr1_chr2_df)
    lapply(seq_len(nrow(chr1_chr2_df)), function(x){
        Brick_filepath <- chr1_chr2_df[x,"filepath"]
        Dimensions <- Brick_matrix_dimensions(Brick = Brick, 
            chr1 = chr1_chr2_df[x,"chr1"], chr2 = chr1_chr2_df[x,"chr2"], 
            resolution = resolution)
        .populate_with_empty_metrics(Brick_filepath = Brick_filepath, 
            chr1 = chr1_chr2_df[x,"chr1"], chr2 = chr1_chr2_df[x,"chr2"],
            row_length = Dimensions[1], col_length = Dimensions[2])
    })
    if(!is.na(keep_chr1)){
        keep_chr2_split <- split(keep_chr2, keep_chr1)
    }

    chr1_chr2_extent_df_list <- lapply(chr1_chr2_rows_df_list, 
        function(a_row){
        chr1 <- unique(a_row$chr1)
        chr1_start <- unique(a_row$chr1_start)
        chr1_end <- unique(a_row$chr1_end)
        read_from <- unique(a_row$read_from)
        read_to <- unique(a_row$read_to)
        chr2_df_split <- split(a_row[,c("chr2", "chr2_start", "chr2_end")],
            a_row[,"chr2"])
        if(!is.na(keep_chr1)){
            if(!(chr1 %in% names(keep_chr2_split))){
                return(NULL)
            }
            chr2_df_split <- chr2_df_split[names(chr2_df_split) %in% 
                keep_chr2_split[[chr1]]]
        }
        values_list <- .return_values_list(mcool_path = mcool_path, 
            mcool_version = mcool_version, read_from = read_from,
            read_to = (read_to - read_from + 1), resolution = resolution,
            has_resolution = has_resolution)
        message("Read ", read_to - read_from + 1, " records...")
        Bin1_id <- values_list[["bin1_id"]] + 1
        Bin2_id <- values_list[["bin2_id"]] + 1
        Counts <- values_list[["counts"]]
        chr2_extent_df_list <- lapply(chr2_df_split, function(chr2_df) {
            chr2 <- unique(chr2_df$chr2)
            chr2_starts <- chr2_df$chr2_start
            chr2_ends <- chr2_df$chr2_end
            chr2_group_path <- Create_Path(
                c(Reference_object$hdf.matrices.root, chr1, chr2))
            chr2_start_df_list <- lapply(seq_along(chr2_starts), function(x){
                chr2_start <- chr2_starts[x]
                chr2_end <- chr2_ends[x]

                Brick_filepath <- unique(a_row$filepath[a_row$chr1 == chr1 & 
                                    a_row$chr2 == chr2])
                Matrix_dimensions <- Brick_matrix_dimensions(Brick = Brick,
                    chr1 = chr1, chr2 = chr2, resolution = resolution)
                Bin2_id_filter <- Bin2_id >= chr2_start & Bin2_id <= chr2_end
                if(!any(Bin2_id_filter)){
                    return(NULL)
                }
                Counts_submat <- Counts[Bin2_id_filter]

                if(!is.null(bias.vector)){
                    Counts_submat <- Counts_submat *
                    bias.vector[Bin1_id[Bin2_id_filter]] *
                    bias.vector[Bin2_id[Bin2_id_filter]]
                }

                Dim <- c(chr1_end - chr1_start + 1, chr2_end - chr2_start + 1)
                Matrix <- matrix(data = 0, nrow = Dim[1], ncol = Dim[2])

                chr1_mat_start <- (chr1_start - 1) - index_chrom_offset[chr1]
                chr2_mat_start <- (chr2_start - 1) - index_chrom_offset[chr2]
                Bin1_id_mat_id <- Bin1_id[Bin2_id_filter] - 
                    index_chrom_offset[chr1] - chr1_mat_start
                Bin2_id_mat_id <- Bin2_id[Bin2_id_filter] - 
                    index_chrom_offset[chr2] - chr2_mat_start
                message("Row segment: ", 
                    paste(chr1, chr1_start, chr1_end, sep = ":"), 
                    "; Col segment: ", 
                    paste(chr2, chr2_start, chr2_end, sep = ":"))
                Matrix[cbind(Bin1_id_mat_id, Bin2_id_mat_id)] <- Counts_submat
                chr1_file_start <- chr1_start - index_chrom_offset[chr1]
                chr2_file_start <- chr2_start - index_chrom_offset[chr2]
                Start <- c(chr1_file_start, chr2_file_start)
                Stride <- c(1, 1)
                Count <- Dim
                if(chr1 == chr2){
                    if(chr1_file_start == chr2_file_start){
                        Matrix[cbind(Bin2_id_mat_id, Bin1_id_mat_id)] <- 
                        Counts_submat
                    }
                    ._Brick_Put_Something_(Group.path = chr2_group_path,
                        data = t(Matrix), Brick = Brick_filepath,
                        Start = rev(Start), Stride = Stride,
                        Count = rev(Count),
                        Name = Reference_object$hdf.matrix.name)
                }
                ._Brick_Put_Something_(Group.path = chr2_group_path, 
                    data = Matrix, Brick = Brick_filepath, Start = Start, 
                    Name = Reference_object$hdf.matrix.name, Stride = Stride, 
                    Count = Count)
                matrix_extent <- .add_metrics(Brick_filepath = Brick_filepath,
                    group_path = chr2_group_path, a_matrix = Matrix, 
                    dimensions = dim(Matrix), 
                    chr1_start = chr1_start - index_chrom_offset[chr1], 
                    chr2_start = chr2_start - index_chrom_offset[chr2])
                data.frame(mcool_name = basename(mcool_path), chr1 = chr1, 
                    chr2 = chr2, min = matrix_extent[1], 
                    max = matrix_extent[2], nrow = Matrix_dimensions[1], 
                    ncol = Matrix_dimensions[2], 
                    distance = max(Matrix_dimensions),
                    stringsAsFactors = FALSE)
            })
            chr2_start_df_list <- chr2_start_df_list[
                !vapply(chr2_start_df_list, is.null, TRUE)]
            chr2_start_df <- do.call(rbind, chr2_start_df_list)
        })
        chr2_extent_df <- do.call(rbind, chr2_extent_df_list)
    })
    chr1_chr2_extent_df <- do.call(rbind, chr1_chr2_extent_df_list)
    return_df <- .final_metrics_processing(Brick = Brick, 
        metrics_df = chr1_chr2_extent_df, 
        resolution = resolution)
    return(return_df)
}