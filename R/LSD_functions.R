._get_first_nonzero_bin <- function(Brick = NULL, chr = NULL, resolution = NA){
    RowSums <- Brick_get_matrix_mcols(Brick = Brick, chr1 = chr, chr2 = chr, 
        resolution = resolution, what = "chr1_row_sums")
    return(min(which(RowSums > 0)))
}
._get_sparsity_index <- function(Brick = NULL, chr = NULL, resolution = NA){
    Sparsity.index <- Brick_get_matrix_mcols(Brick = Brick, chr1 = chr, 
        chr2 = chr, resolution = resolution, what = "sparsity")
    return(Sparsity.index)
}
Backwards.Difference <- function(Vector=NULL,sparse=FALSE,sparsity.idx=NULL,
    sparsity_threshold=NULL){
        if(is.null(Vector)){
            stop("Vector variable cannot be empty.")
        }
        if(sparse){
            VectorDiff <- rep(NA,length(Vector))
            temp.diff <- rev(diff(
                rev(Vector[sparsity.idx > sparsity_threshold])))
            temp.diff[length(temp.diff)+1] <- 0
            VectorDiff[sparsity.idx > sparsity_threshold] <- temp.diff
        }else{
            VectorDiff<-rev(diff(rev(Vector)))
            VectorDiff[length(Vector)]=0
        }
        return(VectorDiff)
}
Forwards.Difference <- function(Vector=NULL,sparse=FALSE,sparsity.idx=NULL,
    sparsity_threshold=NULL){
        if(is.null(Vector)){
            stop("Vector variable cannot be empty.")
        }
        if(sparse){
            VectorDiff <- rep(NA,length(Vector))
            temp.diff <- diff(Vector[sparsity.idx > sparsity_threshold])
            temp.diff <- c(0,temp.diff)
            VectorDiff[sparsity.idx > sparsity_threshold] <- temp.diff
        }else{
        VectorDiff<-diff(Vector)
        VectorDiff=c(0,VectorDiff)
        }
        return(VectorDiff)
}
ComputeOutlierOverIQRInsideWindow <- function(lookup_window=NULL,
    diff.values=NULL,values=NULL, row.sums = NULL,
        min_sum = NULL, tukeys_constant=NULL,tail=NULL,
        sparse=FALSE,strict = FALSE,sparsity.idx=NULL,sparsity_threshold=NULL){
        seq.over.value <- seq_along(diff.values)
        Filter <- row.sums > min_sum
        if(sparse){
            Filter <- Filter & sparsity.idx > sparsity_threshold
        }
        seq.over.value <- seq.over.value[Filter]
        seq.over.seq <- seq_along(seq.over.value)
        outlier.list <- lapply(seq.over.seq,function(x.seq){
            lookup_window.range <- (
                (x.seq - lookup_window) : (x.seq + lookup_window))
            lookup_window.range <- seq.over.value[lookup_window.range[
            lookup_window.range>0 & lookup_window.range<=max(seq.over.seq)]]
            offset <- (min(lookup_window.range)-1)
            diff.value.window <- diff.values[lookup_window.range]
            value.window <- values[lookup_window.range]
            value.quartile<-quantile(diff.value.window,na.rm=TRUE)
            InterQuartile <- value.quartile[4]-value.quartile[2]
            if(tail=="lower.tail"){
                #Calculate Inner fences based on accepted formula
                fences <- value.quartile[2] - (InterQuartile*tukeys_constant)
                Outlier.Filter <- !is.na(value.window) &
                diff.value.window <= fences &
                diff.value.window < value.window
                if(strict){
                    Outlier.Filter <- Outlier.Filter & (value.window < 0)
                }
            }else if(tail=="upper.tail"){
                fences <- value.quartile[4] + (InterQuartile*tukeys_constant)
                Outlier.Filter <- !is.na(value.window) &
                diff.value.window >= fences &
                diff.value.window > value.window
                if(strict){
                    Outlier.Filter <- Outlier.Filter & (value.window > 0)
                }
            }
            outliers<-which(Outlier.Filter)
            if(length(outliers)>0){
                outliers <- lookup_window.range[outliers]
            }
            outliers
        })
        outlier.vector.dups <- do.call(c,outlier.list)
        outlier.vector.uniq.sorted <- sort(unique(outlier.vector.dups))
        return(outlier.vector.uniq.sorted)
}
CreateDomainlist <- function(start.vector=NULL,end.vector=NULL,fill_gaps=NULL){
    Domains.by.start.list <- lapply(start.vector,function(x){
        data.frame(startbin=x, endbin=min(end.vector[end.vector > x]))
    })
    Domains.by.start.df <- do.call(rbind,Domains.by.start.list)
    Domains.by.start.df$gap.fill <- as.numeric(FALSE)
    Domains.by.start.df$level <- 2
    Domains.by.end.df <- NULL
    Domains.by.assumption.df <- NULL
    if(fill_gaps){
        uncurated.ends <- end.vector[
        !(end.vector %in% Domains.by.start.df[,"endbin"])]
        if(length(uncurated.ends) > 0){
            Domains.by.end.list <- lapply(uncurated.ends,function(x){
                data.frame(startbin=(
                    max(Domains.by.start.df$endbin[
                        Domains.by.start.df$endbin<x])+1),endbin=x)
            })
            Domains.by.end.df <- do.call(rbind,Domains.by.end.list)
            Domains.by.end.df$gap.fill <- as.numeric(TRUE)
            Domains.by.end.df$level <- 1
        }
    }
    All.Domains <- rbind(Domains.by.start.df,
        Domains.by.end.df,
        Domains.by.assumption.df)
    All.Domains.sorted <- All.Domains[order(All.Domains$startbin),]
    return(All.Domains.sorted)
}
ComputeDirectionalityIndex <- function(Matrix = NULL, Window.size=NULL, 
    filter = NULL, start = NULL, end = NULL){
    Sequence <- seq_len(nrow(Matrix))
    Sequence <- Sequence[start:end]
    DI.Data <- rep(NA,length(Sequence))
    Bins.to.process <- Sequence[filter[start:end]]
    All.bins <- seq_len(nrow(Matrix))
    All.bins <- All.bins[filter]
    DI.list <- vapply(Bins.to.process,function(i){
            Upstream<-0
            Downstream<-0
            My.DI.Data <- NA
            Relative.mid <- which(All.bins == i)
            Window.range <- c(
                (Relative.mid - Window.size) : (Relative.mid + Window.size))
            Window.range <- All.bins[
            Window.range[Window.range >= 1 & Window.range <= length(All.bins)]]
            Upstream.range <- Window.range[Window.range < i]
            Downstream.range <- Window.range[Window.range > i]
            Row.vector <- Matrix[i,]
            Row.vector[is.na(Row.vector) | is.infinite(Row.vector)] <- 0
            if(length(Upstream.range) > 0){
                Upstream <- sum(Row.vector[Upstream.range])
            }
            if(length(Downstream.range) > 0){
                Downstream <- sum(Row.vector[Downstream.range])
            }
            Expected <- (Upstream + Downstream)/2
            if( Expected == 0 | Upstream == Downstream ){
                My.DI.Data <- 0
            }else{
                My.DI.Data <- (
                    (Downstream - Upstream)/abs(Downstream - Upstream)
                    ) * (
                    ((Upstream - Expected)^2)/Expected + (
                        (Downstream - Expected)^2)/Expected)
            }
        },1)
    DI.Data[Sequence %in% Bins.to.process] <- DI.list
    return(DI.Data)
}
get_directionality_index_by_chunks <- function(Brick = NULL, chr = NULL,
    resolution = NA, di_window = NULL, distance = NULL, chunk_size = 500, 
    sparse = FALSE, sparsity_threshold = 0.8, min_sum = -1, force = FALSE){
    Ranges <- Brick_get_bintable(Brick = Brick, chr = chr, 
        resolution = resolution)
    First.non.zero.bin <- ._get_first_nonzero_bin(Brick = Brick, chr = chr,
        resolution = resolution)
    chr.length <- length(Ranges)
    RowSums <- Brick_get_matrix_mcols(Brick = Brick, chr1 = chr, chr2 = chr, 
        resolution = resolution, what = "chr1_row_sums")
    if(sparse){
        SparsityIndex <- Brick_get_matrix_mcols(Brick = Brick, chr1 = chr, 
            chr2 = chr, resolution = resolution, what = "sparsity")
    }
    if((chunk_size - (di_window*2))/di_window < 10){
        stop("chunk_size is too small for this di_window\n")
    }
    if(any(di_window > distance)){
        stop("di_window cannot be larger than distance\n")
    }
    Span <- (chr.length - First.non.zero.bin)
    Iterations <- Span/chunk_size
    Starts <- seq(from = First.non.zero.bin, to = chr.length, by = chunk_size)
    Starts <- Starts[Starts != chr.length]
    Ends <- c(Starts[-1] -1, chr.length)
    DI.data.list <- lapply(seq_along(Starts), function(x){
        Start <- Starts[x]
        End <- Ends[x]
        Position.start <- Start
        Position.end <- End
        Start <- ifelse((Start - di_window) < First.non.zero.bin, 
            First.non.zero.bin, (Start - di_window))
        End <- ifelse((End + di_window) > chr.length, 
            chr.length, (End + di_window))
        RowSums.subset <- RowSums[Start:End]
        Filter <- RowSums.subset > min_sum
        if(sparse){
            Sparsity.index.subset <- SparsityIndex[Start:End]
            Filter <- Filter & (Sparsity.index.subset > sparsity_threshold)
        }
        Total.length <- length(Filter)
        Filter.extend.length <- length(which(!Filter))
        True.length <- length(which(Filter))
        extend <- Filter.extend.length
        while(True.length < Total.length){
            Start <- ifelse((Start - extend) < First.non.zero.bin, 
                First.non.zero.bin, Start - extend)
            End <- ifelse((End + extend) > chr.length, 
                chr.length, (End + extend))
            RowSums.subset <- RowSums[Start:End]
            Filter <- RowSums.subset > min_sum
            if(sparse){
                Sparsity.index.subset <- SparsityIndex[Start:End]
                Filter <- Filter & (Sparsity.index.subset > sparsity_threshold)
            }
            True.length <- length(which(Filter))
            extend <- extend + 1
        }
        Matrix <- Brick_get_vector_values(Brick = Brick, chr1 = chr, 
            resolution = resolution, chr2 = chr, xaxis=c(Start:End), 
            yaxis=c(Start:End), force = force)
        # cat((Start - 1),"\n")
        # message(Position.start," ",Position.end,"\n")
        # message(Position.start - (Start - 1),Position.end - (Start - 1),"\n")
        DI.data <- ComputeDirectionalityIndex(Matrix = Matrix, 
            Window.size = di_window, filter = Filter, 
            start = Position.start - (Start - 1), 
            end = Position.end - (Start - 1))
        return(DI.data)
    })
    DI.data <- do.call(c, DI.data.list)
    DI.data <- c(rep(NA,First.non.zero.bin - 1),DI.data)
    Ranges$DI.Data <- DI.data
    return(Ranges)
}
MakeBoundaries <- function(chr = NULL, Ranges = NULL, Binsize = NULL){
    Ends <- end(Ranges) 
    Ends <- Ends[seq_len(length(Ends)-1)]
    Starts <- start(Ranges)
    Starts <- Starts[seq_len(length(Starts))[-1]] - 1 
    Domain.boundaries <- unique(Starts,Ends)
    Boundary.ranges <- Brick_make_ranges(
        chrom=rep(chr,length(Domain.boundaries)),
        start=(Domain.boundaries-(Binsize/2))+1,
        end=Domain.boundaries+(Binsize/2))
}

#' Do TAD Calls with Local Score Differentiator on a Hi-C matrix
#' 
#' `Local_score_differentiator` calls topologically associated domains on Hi-C 
#' matrices. Local score differentiator at the most fundamental level is a 
#' change point detector, which detects change points in the directionality 
#' index using various thresholds defined on a local directionality index 
#' distributions.
#' The directionality index (DI) is calculated as defined by Dixon et al., 2012 
#' Nature. Next, the difference of DI is calculated between neighbouring bins to
#' get the change in DI distribution in each bin. When a DI value goes from a
#' highly negative value to a highly positive one as expected to occur at domain
#' boundaries, the ensuing DI difference distribution becomes a very flat 
#' distribution interjected by very large peaks signifying regions where such
#' a change may take place. We use two difference vectors, one is the difference
#' vector between a bin and its adjacent downstream bin and another is the 
#' difference between a bin and its adjacent upstream bin. Using these vectors,
#' and the original directionality index, we define domain borders as outliers.
#' 
#' To define an outlier, fences are first defined. The fences are defined using
#' tukeys_constant x inter-quartile range of the directionality index. The upper
#' fence used for detecting domain starts is the 75th quartile + 
#' (IQR x tukeys_constant), while the lower fence is the 
#' 25th quartile - (IQR x tukeys_constant). For domain starts the DI difference
#' must be greater than or equal to the upper fence, it must be greater than the
#' DI and the DI must be a finite real value. If strict is TRUE, DI will also
#' be required to be greater than 0. Similarly, for domain ends the 
#' DI difference must be lower than or equal to the lower fence, it must be 
#' lower than the DI and the DI must be a finite real value. If strict is TRUE,
#' DI will also be required to be lower than 0. 
#' 
#' After defining outliers, each domain start will be associated to its 
#' nearest downstream domain end. If \emph{fill_gaps} is defined as TRUE and
#' there are domain ends which remain unassociated to a domain start, These 
#' domain ends will be associated to the bin adjacent to their nearest upstream
#' domain end. This associations will be marked by metadata columns, gap.fill= 1
#' and level = 1.
#' 
#' This function provides the capability to call very accurante TAD definitions
#' in a very fast way. 
#' 
#' @inheritParams Brick_get_chrominfo
#' @inheritParams Brick_add_ranges
#' 
#' @param chrs \strong{Optional}. Default NULL
#' If present, only TAD calls for elements in \emph{chrs} will be done.
#' 
#' @param min_sum \strong{Optional}. Default -1
#' Process bins in the matrix with row.sums greater than \emph{min_sum}.
#' 
#' @param di_window \strong{Optional}. Default 200
#' Use \emph{di_window} to define the directionality index.
#' 
#' @param lookup_window \strong{Optional}. Default 200
#' Use \emph{lookup_window} local window to call borders. At smaller 
#' \emph{di_window} values we recommend setting this to 2*\emph{di_window}
#' 
#' @param tukeys_constant \strong{Optional}. Default 1.5
#' \emph{tukeys_constant}*IQR (inter-quartile range) defines the lower and upper
#' fence values.
#' 
#' @param strict \strong{Optional}. Default TRUE
#' If TRUE, \emph{strict} creates an additional filter on the directionality 
#' index requiring it to be either greater than or less than 0 on the right tail
#' or left tail respectively.  
#' 
#' @param fill_gaps \strong{Optional}. Default TRUE
#' If TRUE, this will affect the TAD stiching process. All Border starts are 
#' stiched to the next downstream border ends. Therefore, at times border ends 
#' remain unassociated to a border start. These border ends are stiched to the 
#' adjacent downstream bin from their upstream border end when \emph{fill_gaps} 
#' is true. 
#' 
#' TADs inferred in this way will be annotated with two metadata columns in the 
#' GRanges object. \emph{gap.fill} will hold a value of 1 and \emph{level} will 
#' hold a value 1. TADs which were not filled in will hold a gap.fill value of
#' 0 and a level value of 2.
#' 
#' @param ignore_sparse \strong{Optional}. Default TRUE
#' If TRUE, a matrix which has been defined as sparse during the matrix loading
#' process will be treated as a dense matrix. The \emph{sparsity_threshold} 
#' filter will not be applied. Please note, that if a matrix is defined as 
#' sparse and fill_gaps is TRUE, fill_gaps will be turned off.
#' 
#' @param sparsity_threshold \strong{Optional}. Default 0.8
#' Sparsity threshold relates to the sparsity index, which is computed as the 
#' number of non-zero bins at a certain distance from the diagonal. If a matrix
#' is sparse and ignore_sparse is FALSE, bins which have a sparsity index value
#' below this threshold will be discarded from DI computation.
#' 
#' @param remove_empty Not implemented.
#' After implementation, this will ensure that the presence of centromeric 
#' regions is accounted for.
#' 
#' @param chunk_size \strong{Optional}. Default 500
#' The size of the matrix chunk to process. This value should be larger than 2x
#' di_window.
#' 
#' @param force_retrieve \strong{Optional}. Default TRUE
#' If TRUE, this will force the retrieval of a matrix chunk even when the 
#' retrieval includes interaction points which were not loaded into a Brick 
#' store (larger chunks). Please note, that this does not mean that DI can be 
#' computed at distances larger than max distance. Rather, this is meant to aid
#' faster computation.
#' 
#' @return A ranges object containing domain definitions. The starts and ends
#' of the ranges coincide with the starts and ends of their contained bins from 
#' the bintable. 
#' 
#' @examples
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "lsd_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr3R.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr3R", 
#' chr2 = "chr3R", matrix_file = Matrix_file, delim = " ",
#' remove_prior = TRUE, resolution = 100000)
#' 
#' TAD_ranges <- Brick_local_score_differentiator(Brick = My_BrickContainer, 
#' chrs = "chr3R", resolution = 100000, di_window = 10, lookup_window = 30, 
#' strict = TRUE, fill_gaps = TRUE, chunk_size = 500)
Brick_local_score_differentiator <- function(Brick, chrs = NULL, 
    resolution = NA, all_resolutions = FALSE, min_sum = -1, di_window = 200L, 
    lookup_window = 200L, tukeys_constant=1.5, strict = TRUE, fill_gaps=TRUE, 
    ignore_sparse=TRUE, sparsity_threshold=0.8, remove_empty = NULL, 
    chunk_size = 500, force_retrieve = TRUE){
    BrickContainer_resolution_check(resolution, all_resolutions)
    ChromInfo <- Brick_get_chrominfo(Brick = Brick, 
        resolution = resolution)
    Chromosomes <- ChromInfo[,'chr']
    if(!is.null(chrs)){
        Chromosomes <- ChromInfo[ChromInfo[,'chr'] %in% chrs,'chr']
    }
    chr_done_filter <- vapply(Chromosomes, function(chr){
        Brick_matrix_isdone(Brick = Brick, chr1 = chr, chr2 = chr, 
            resolution = resolution)
    })
    if(!all(chr_done_filter)){
        message("Skipping intra-chromosomal maps containing no data...")
        message(paste(Chromosomes[!chr_done_filter], collapse = ", "), 
            " will be skipped")
    }
    Chrom.domains.ranges.list <- lapply(Chromosomes, function(chr){
        Ranges <- Brick_get_bintable(Brick = Brick, chr = chr, 
            resolution = resolution)
        sparse <- Brick_matrix_issparse(Brick = Brick, chr1 = chr, chr2 = chr, 
            resolution = resolution)
        max.distance <- Brick_matrix_maxdist(Brick = Brick, chr1 = chr, 
            chr2 = chr, resolution = resolution)
        if(ignore_sparse){
            sparse=FALSE
        }
        if(sparse & fill_gaps){
            fill_gaps=FALSE
        }
        message("[1] Computing DI for ",chr,"\n")
        Ranges <- get_directionality_index_by_chunks(Brick = Brick, 
            chr = chr, 
            resolution = resolution,
            di_window = di_window, 
            distance = max.distance, 
            chunk_size = chunk_size, 
            sparse=sparse, 
            sparsity_threshold=sparsity_threshold,
            min_sum = min_sum, force = force_retrieve)

        RowSums <- Brick_get_matrix_mcols(Brick = Brick, chr1 = chr, 
            chr2 = chr, resolution = resolution, what = "chr1_row_sums")
        Ranges$row.sums <- RowSums
        message("[2] Computing DI Differences for ",chr,"\n")
        if(sparse){
            SparsityIndex <- Brick_get_matrix_mcols(Brick = Brick, 
                chr1 = chr, 
                chr2 = chr, 
                resolution = resolution, 
                what = "sparsity")
            Backwards.DI.Difference <- Backwards.Difference(
                Vector = Ranges$DI.Data, 
                sparse = sparse,
                sparsity.idx = SparsityIndex, 
                sparsity_threshold = sparsity_threshold)
            Forwards.DI.Difference <- Forwards.Difference(
                Vector = Ranges$DI.Data, 
                sparse = sparse,
                sparsity.idx = SparsityIndex, 
                sparsity_threshold = sparsity_threshold)
        }else{
            Backwards.DI.Difference <- Backwards.Difference(
                Vector=Ranges$DI.Data)
            Forwards.DI.Difference <- Forwards.Difference(
                Vector=Ranges$DI.Data)
        }
        Ranges$backward.Differences <- Backwards.DI.Difference
        Ranges$forward.Differences <- Forwards.DI.Difference
        message("[2] Done\n")
        message("[3] Fetching Outliers ",chr,"\n")
        if(sparse){
            Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(
                lookup_window = lookup_window,
                diff.values = Backwards.DI.Difference, 
                values = Ranges$DI.Data, 
                sparse = sparse, 
                row.sums = Ranges$row.sums,
                min_sum = min_sum, 
                sparsity.idx = SparsityIndex, 
                sparsity_threshold = sparsity_threshold, 
                tukeys_constant = tukeys_constant, 
                tail = "lower.tail",
                strict = strict)
            Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(
                lookup_window = lookup_window,
                diff.values = Forwards.DI.Difference, 
                values = Ranges$DI.Data, 
                sparse = sparse, 
                row.sums = Ranges$row.sums,
                min_sum = min_sum, 
                sparsity.idx = SparsityIndex, 
                sparsity_threshold = sparsity_threshold,
                tukeys_constant = tukeys_constant, 
                tail = "upper.tail", 
                strict = strict)
        }else{
            Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(
                lookup_window=lookup_window,
                diff.values=Backwards.DI.Difference,
                values=Ranges$DI.Data, 
                row.sums = Ranges$row.sums,
                min_sum = min_sum, 
                tukeys_constant=tukeys_constant,
                tail="lower.tail",strict=strict)
            Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(
                lookup_window=lookup_window,
                diff.values=Forwards.DI.Difference,
                values=Ranges$DI.Data, 
                row.sums = Ranges$row.sums,
                min_sum = min_sum, 
                tukeys_constant=tukeys_constant,
                tail="upper.tail",
                strict=strict)
        }
        Domain.start.candidates <- Domain.start.candidates[
        Domain.start.candidates != length(Ranges)]
        Domain.end.candidates <- Domain.end.candidates[
        Domain.end.candidates != 1]
        message("[3] Done\n")
        message("[4] Creating Domain list for ",chr,"\n")

        if(!(1 %in% Domain.start.candidates)){
            Domain.start.candidates <- c(1,Domain.start.candidates)
        }
        if(!(length(Ranges) %in% Domain.end.candidates)){
            Domain.end.candidates <- c(Domain.end.candidates,length(Ranges))
        }
        Domain.list <- CreateDomainlist(start.vector=Domain.start.candidates,
            end.vector=Domain.end.candidates, fill_gaps=fill_gaps)
        Domain.Ranges <- Brick_make_ranges(chrom=rep(chr,nrow(Domain.list)),
            start=start(Ranges[Domain.list$startbin]),
            end=end(Ranges[Domain.list$endbin]))
        message("[4] Done\n")
        Domain.Ranges$gap.fill <- Domain.list$gap.fill
        Domain.Ranges$level <- Domain.list$level
        Domain.Ranges$window.size <- di_window
        Domain.Ranges$lookup_window <- lookup_window
        return(Domain.Ranges)
    })
    Chrom.domains.ranges <- do.call(c,unlist(Chrom.domains.ranges.list, 
        use.names = TRUE))
    return(Chrom.domains.ranges)
}
