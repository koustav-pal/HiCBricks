._Process_mcool <- function(Brick = NULL, File = NULL,
    cooler.batch.size = 1000000, binsize = NULL, resolution = FALSE,
    matrix.chunk = 2000, group = NULL, chr1 = NULL,
    chr2 = NULL, dont.look.for.chr2 = FALSE, norm.factor = NULL){
    is.sparse <- FALSE
    Reference.object <- GenomicMatrix$new()
    Indexes.Keys <- Reference.object$mcool.index.keys()
    Indexes.group <- Indexes.Keys[1]
    Indexes.chrom.key <- Indexes.Keys[3]
    Indexes.bin.key <- Indexes.Keys[2]

    mcool.version <- GetAttributes(Path = NULL, File=File,
        Attributes="format-version", on = "file",
        ignore.fun.cast = TRUE)[,"format-version"]
    mcool.version <- as.numeric(as.character(mcool.version))
    Bintable.keys <- Reference.object$mcool.bintable.keys(
        version = mcool.version)
    cooler.remap.chrom <- ._mcool_remap_chromosomes(File = File,
        resolution = resolution, binsize = binsize,
        mcool.version = mcool.version)
    Chrom.names <- cooler.remap.chrom[,"chr.name"]
    Chrom.pos <- cooler.remap.chrom[,"remap.chrom"]
    names(Chrom.pos) <- Chrom.names
    chr1.pos <- Chrom.pos[chr1]
    chr2.pos <- Chrom.pos[chr2]
    if(chr2.pos < chr1.pos){
        Which.less <- which(chr2.pos < chr1.pos)
        Left.over <- chr1[Which.less]
        chr1[Which.less] <- chr2[Which.less]
        chr2[Which.less] <- Left.over
    }

    Bintable.group <- Bintable.keys[1]
    if(resolution){
        Indexes.group.path <- Create_Path(
            c(
                Reference.object$mcool.resolutions.name,
                binsize,Indexes.group
                )
            ) 
        Bintable.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Bintable.group))
        # Scaffold.group.path <-  Create_Path(
        # c(Reference.object$mcool.resolutions.name,binsize,Scaffold.group))
    }else{
        Indexes.group.path <- Create_Path(c(Indexes.group))
        Bintable.group.path <- Create_Path(Bintable.group)
        # Scaffold.group.path <-  Create_Path(Scaffold.group)
    }
    bias.vector <- NULL
    if(!is.null(norm.factor)){
        bias.vector <- ._Brick_Get_Something_(Group.path = Bintable.group.path,
            Brick = File, Name = norm.factor, return.what = "data")
    }
    Index.chrom.offset<-._Brick_Get_Something_(Group.path = Indexes.group.path,
        Brick = File, Name = Indexes.chrom.key, return.what = "data")
    names(Index.chrom.offset) <- Chrom.names
    Index.bin.offset <- ._Brick_Get_Something_(Group.path = Indexes.group.path,
        Brick = File, Name = Indexes.bin.key, return.what = "data")
    names(Index.bin.offset) <- Chrom.names

    chrom1.loc <- which(names(Index.chrom.offset) %in% chr1)
    Start.of.chrom1.reads <- ifelse(Index.chrom.offset[chrom1.loc]==0, 1,
        Index.bin.offset[Index.chrom.offset[chrom1.loc]]+1)

    Chrom1.iter.list <- make_mcool_iterations(
        Start.pos = Index.chrom.offset[chrom1.loc]+1,
        End.pos = Index.chrom.offset[chrom1.loc+1],
        step = matrix.chunk)
    End.of.chrom1.reads <- Index.bin.offset[Index.chrom.offset[chrom1.loc+1]]

    chrom2.loc <- which(names(Index.chrom.offset) %in% chr2)
    group.path <- Create_Path(c(Reference.object$hdf.matrices.root, 
        chr1, chr2))
    Starts <- Chrom1.iter.list[["start"]]
    Ends <- Chrom1.iter.list[["end"]]
    Chrom2.iter.list <- make_mcool_iterations(
        Start.pos = Index.chrom.offset[chrom2.loc]+1,
        End.pos = Index.chrom.offset[chrom2.loc+1],
        step = matrix.chunk)
    Chrom2.starts <- Chrom2.iter.list[["start"]]
    Chrom2.ends <-  Chrom2.iter.list[["end"]]

    ChromInfo <- Brick_get_chrominfo(Brick = Brick)
    Max.End.Offset <- (max(Ends) - Index.chrom.offset[chrom1.loc])
    if(!(Max.End.Offset %in% ChromInfo[ChromInfo[,"chr"]==chr1,"nrow"])){
        stop("Matrix dimensions do not ",
            "match to the provided Brick dimensions!\n")
    }
    Iterations <- make_mcool_iterations(Start.pos=Start.of.chrom1.reads,
        End.pos = End.of.chrom1.reads, step=cooler.batch.size)
    Iter.starts <- Iterations[["start"]]
    Iter.ends <- Iterations[["end"]]


    if(!dont.look.for.chr2){
        Start.of.chrom2.reads <- find_chr2_start_position(File = File,
            iterations = Iterations, start.1 = Starts[1],
            end.1 = Ends[length(Ends)], start.2 = Chrom2.starts[1],
        end.2 = Chrom2.ends[length(Chrom2.ends)],
        mcool.version = mcool.version, resolution = resolution,
        binsize = binsize)

        Iterations <- make_mcool_iterations(Start.pos=Start.of.chrom2.reads,
            End.pos = End.of.chrom1.reads, step=cooler.batch.size)
        Iter.starts <- Iterations[["start"]]
        Iter.ends <- Iterations[["end"]]
    }

    metrics.list <- prepare_empty_metrics_list(starts.1 = Starts,
        ends.1 = Ends, starts.2 = Chrom2.starts, ends.2 = Chrom2.ends,
        chrom1 = chr1, chrom2 = chr2)
    metrics.list <- populate_matrix_with_values(Brick = Brick, File = File,
        group.path = group.path, starts.1 = Starts, ends.1 = Ends,
        starts.2 = Chrom2.starts, ends.2 = Chrom2.ends, chrom1 = chr1,
        chrom1.offset = Index.chrom.offset[chrom1.loc],
        chrom2.offset = Index.chrom.offset[chrom2.loc], chrom2 = chr2,
        iter.start = Iter.starts, iter.end = Iter.ends,
        mcool.version = mcool.version, bias.vector = bias.vector,
        metrics.list = metrics.list, resolution = resolution, 
        binsize = binsize)

    distance <- max(Chrom2.ends) - Index.chrom.offset[chrom2.loc]
    Matrix.range <- metrics.list[["extent"]]
    Attributes <- Reference.object$matrices.chrom.attributes
    Attr.vals <- c(basename(File),as.double(Matrix.range),
        as.integer(is.sparse),as.integer(distance),
        as.integer(TRUE))

    if(is.list(metrics.list[["row.sums"]])){
        chr1.length <- length(metrics.list[["row.sums"]][[chr1]])
        chr2.length <- length(metrics.list[["row.sums"]][[chr2]])
        ._Brick_WriteArray_(Brick = Brick, Path = group.path,
            name = Reference.object$hdf.matrix.rowSums,
            object = metrics.list[["row.sums"]][[chr1]])
        ._Brick_WriteArray_(Brick = Brick, Path = group.path,
            name = Reference.object$hdf.matrix.rowSums,
            object = metrics.list[["bin.coverage"]][[chr1]]/chr1.length)
        WriteAttributes(Path = group.path, File = Brick,
            Attributes = Attributes, values = Attr.vals, on = "group")
        group.path <- Create_Path(
            c(Reference.object$hdf.matrices.root,
                chr2,
                chr1))
        ._Brick_WriteArray_(Brick = Brick, Path = group.path,
            name = Reference.object$hdf.matrix.rowSums,
            object = metrics.list[["row.sums"]][[chr2]])
        ._Brick_WriteArray_(Brick = Brick, Path = group.path,
            name = Reference.object$hdf.matrix.coverage,
            object = metrics.list[["bin.coverage"]][[chr2]]/chr2.length)
        WriteAttributes(Path = group.path, File = Brick,
            Attributes = Attributes, values = Attr.vals, on = "group")
    }else{
        chr1.length <- length(metrics.list[["row.sums"]])
        ._Brick_WriteArray_(Brick = Brick, Path = group.path,
            name = Reference.object$hdf.matrix.rowSums,
            object = metrics.list[["row.sums"]])
        ._Brick_WriteArray_(Brick = Brick, Path = group.path,
            name = Reference.object$hdf.matrix.coverage,
            object = metrics.list[["bin.coverage"]]/chr1.length)
        WriteAttributes(Path = group.path, File = Brick,
            Attributes = Attributes, values = Attr.vals, on = "group")
    }
    return(TRUE)
}

make_mcool_iterations <- function(Start.pos = NULL, End.pos = NULL,
    step = NULL){
    Starts <- seq(from = Start.pos, to = End.pos, by = step)
    Starts <- Starts[Starts != End.pos]
    Ends <- c(Starts[-1] - 1, End.pos)
    return(list(start = Starts, end = Ends))
}

find_chr2_start_position <- function(File = NULL, iterations = NULL,
    resolution = FALSE, binsize = NULL, start.1 = NULL, end.1 = NULL,
    start.2 = NULL, end.2 = NULL, mcool.version = NULL){
    message("Finding start of chr2 values...\n")
    Reference.object <- GenomicMatrix$new()
    Matrix.Keys <- Reference.object$mcool.matrix.keys(version = mcool.version)
    if(resolution){
        Matrix.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Matrix.Keys[1]))
    }else{
        Matrix.group.path <- Create_Path(c(Matrix.Keys[1]))
    }
    Start <- NULL
    k<- 1
    iter.start <- iterations[["start"]]
    iter.end <- iterations[["end"]]
    while(is.null(Start)){
        Bin2_id <- ._Brick_Get_Something_(Group.path = Matrix.group.path,
            Brick = File, Index = list(c(iter.start[k]:iter.end[k])),
            Name = Matrix.Keys[3], return.what = "data")
        Bin2_id <- Bin2_id + 1
        Bin1_id <- ._Brick_Get_Something_(Group.path = Matrix.group.path,
            Brick = File, Index = list(c(iter.start[k]:iter.end[k])),
            Name = Matrix.Keys[2], return.what = "data")
        Bin1_id <- Bin1_id + 1
        if(any(Bin2_id > start.2)){
            if(all(Bin1_id > end.1)){
                message(start.1,"\n")
                stop("bin1 ids exceeded chr1 reads!")
            }
            if(all(Bin2_id > end.2)){
                stop("bin2 ids exceeded chr2 reads!")
            }
            Start <- (iter.start[k] - 1) + min(which(Bin2_id > start.2))
        }
        k <- k + 1
    }
    return(Start)
}

populate_matrix_with_values <- function(Brick = NULL, File = NULL,
    group.path = NULL, chrom1 = NULL, starts.1 = NULL, chrom1.offset = NULL,
    ends.1 = NULL, starts.2 = NULL, ends.2 = NULL, chrom2 = NULL,
    chrom2.offset = NULL, iter.start = NULL, iter.end = NULL,
    mcool.version = NULL, metrics.list = NULL,
    resolution = FALSE, binsize = NULL, bias.vector = NULL){
    Reference.object <- GenomicMatrix$new()
    Matrix.Keys <- Reference.object$mcool.matrix.keys(version = mcool.version)
    x <- 1
    if(resolution){
        Matrix.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Matrix.Keys[1]))
    }else{
        Matrix.group.path <- Create_Path(c(Matrix.Keys[1]))
    }
    while(x <= length(starts.1)){
        y <- 1
        starts.2.sub <- starts.2[starts.2 >= starts.1[x]]
        ends.2.sub <- ends.2[starts.2 >= starts.1[x]]
        while(y <= length(starts.2.sub)){
            Message <- paste(paste("Pair: ",paste(chrom1,
                paste(chrom2,sep = "."),sep = ", "),".",sep = ""),
                "Reading row",paste(starts.1[x] - chrom1.offset,
                    "to", ends.1[x] - chrom1.offset), "of",
                max(ends.1) - chrom1.offset,"rows,", "and col",
                paste(starts.2.sub[y] - chrom2.offset, "to",
                    ends.2.sub[y] - chrom2.offset),
                "of", max(ends.2.sub) - chrom2.offset, "cols.")

            message(Message,"\n")
            rowspan <- (ends.1[x] - starts.1[x]) + 1
            colspan <- (ends.2.sub[y] - starts.2.sub[y]) + 1
            col.offset <- starts.2.sub[y] - chrom2.offset - 1
            row.offset <- starts.1[x] - chrom1.offset - 1
            Matrix <- matrix(0, nrow = rowspan, ncol = colspan)
            m <- 1
            next_segment <- FALSE
            Bin1_id.c <- NULL
            Bin2_id.c <- NULL
            while(!next_segment & m <= length(iter.start)){
                Index <- c(iter.start[m]:iter.end[m])
                Bin1_id <- ._Brick_Get_Something_(
                    Group.path = Matrix.group.path,
                    Brick = File, Index = list(Index), Name = Matrix.Keys[2],
                    return.what = "data")
                Bin1_id <- Bin1_id + 1
                Bin2_id <- ._Brick_Get_Something_(
                    Group.path = Matrix.group.path,
                    Brick = File, Index = list(Index), Name = Matrix.Keys[3],
                    return.what = "data")
                Bin2_id <- Bin2_id + 1
                message("Read ",length(c(iter.start[m]:iter.end[m])),
                    " records.\n")
                Filter.1 <- Bin1_id >= starts.1[x] & Bin1_id <= ends.1[x]
                Filter.2 <- Bin2_id >= starts.2.sub[y] &
                Bin2_id <= ends.2.sub[y]
                Filter <- Filter.1 & Filter.2
                Bin1_id_norm <- Bin1_id[Filter]
                Bin2_id_norm <- Bin2_id[Filter]
                Bin1_id <- Bin1_id[Filter] - chrom1.offset - row.offset
                Bin2_id <- Bin2_id[Filter] - chrom2.offset - col.offset
                if(any(length(Bin1_id) > 0)){
                    Counts <- ._Brick_Get_Something_(
                        Group.path = Matrix.group.path, Brick = File,
                    Index = list(Index[Filter]), Name = Matrix.Keys[4],
                    return.what = "data")
                    if(!is.null(bias.vector)){

                        Counts <- Counts *
                        bias.vector[Bin1_id_norm] *
                        bias.vector[Bin2_id_norm]
                    }
                    Matrix[cbind(Bin1_id,Bin2_id)] <- Counts
                    Bin1_id.c <- c(Bin1_id.c,Bin1_id)
                    Bin2_id.c <- c(Bin2_id.c,Bin2_id)
                }else{
                    next_segment <- TRUE
                }
                m <- m+1
            }
            Start <- c(starts.1[x] - chrom1.offset, starts.2.sub[y] -
                chrom2.offset)
            Stride <- c(1,1)
            mat.Count <- dim(Matrix)
            if(m > 1){
                metrics.list <- insert_data_and_computemetrics_both_matrices(
                    Brick = Brick, Matrix = Matrix, group.path = group.path,
                    chrom1 = chrom1, chrom2 = chrom2, row.offset = row.offset,
                    col.offset = col.offset, row.pos = Bin1_id.c,
                    col.pos = Bin2_id.c, metrics.list = metrics.list)
            }
            y <- y + 1
        }
        x <- x + 1
    }
    return(metrics.list)
}

add_to_data <- function(Vector = NULL, start = NULL, end = NULL,
    data = NULL){
    Vector[start:end] <- Vector[start:end] + data
    return(Vector)
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
        Start <- rev(Start)
        Stride <- c(1,1)
        Count <- rev(Count)
        ._Brick_Put_Something_(Group.path = Create_Path(
            c(Reference.object$hdf.matrices.root,chrom2,chrom1)),
            Brick = Brick, Name = dataset.name, data = t(Matrix), 
            Start = Start, Stride = Stride, Count = Count)
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

._mcool_remap_chromosomes = function(File=NULL, mcool.version = NULL,
    resolution = FALSE, binsize = NULL){
    Reference.object <- GenomicMatrix$new()
    Bintable.keys <- Reference.object$mcool.bintable.keys(
        version = mcool.version)
    Scaffold.keys <- Reference.object$mcool.scaffold.keys(
        version = mcool.version)
    Bintable.group <- Bintable.keys[1]
    Bintable.chr <- Bintable.keys[2]
    Scaffold.group <- Scaffold.keys[1]
    Scaffold.name <- Scaffold.keys[3]
    if(resolution){
        Bintable.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Bintable.group)) 
        Scaffold.group.path <-  Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Scaffold.group))
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
    Reference.object <- GenomicMatrix$new()
    Bintable.keys <- Reference.object$mcool.bintable.keys(
        version = mcool.version)
    Bintable.group <- Bintable.keys[1]
    Bintable.chr <- Bintable.keys[2]
    Bintable.start <- Bintable.keys[3]
    Bintable.end <- Bintable.keys[4]
    if(resolution){
        Bintable.group.path <- Create_Path(
            c(Reference.object$mcool.resolutions.name,binsize,Bintable.group)) 
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

mcool_list_resolutions <- function(mcool = NULL){
    Reference.object <- GenomicMatrix$new()
    Handler <- ReturnH5FileConnection(File = mcool)
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)[,"name"]
    CloseH5Con(Handle = Handler, type = "file")
    if(Reference.object$mcool.resolutions.name %in% GroupList){
        Handler <- ._Brick_Get_Something_(
            Group.path = Create_Path(Reference.object$mcool.resolutions.name),
            Brick = mcool, return.what = "group_handle")
        BinList <- h5ls(Handler,
            datasetinfo = FALSE, recursive = FALSE)[,"name"]
        CloseH5Con(Handle = Handler, type = "group")
    }else{
        BinList <- NULL
    }
    return(.format_resolution(BinList))
}
