._get_first_nonzero_bin <- function(Lego = NULL, chr = NULL){
	RowSums <- Lego_get_matrix_mcols(Lego = Lego, chr1 = chr, chr2 = chr, what = "row.sums")
    return(min(which(RowSums > 0)))
}
._get_sparsity_index <- function(Lego = NULL, chr = NULL){
    Sparsity.index <- Lego_get_matrix_mcols(Lego = Lego, chr1 = chr, chr2 = chr, what = "sparsity")
    return(Sparsity.Index)
}
Backwards.Difference <- function(Vector=NULL,sparse=FALSE,sparsity.idx=NULL,sparsity.threshold=NULL){
        if(is.null(Vector)){
            stop("Vector variable cannot be empty.")
        }
        if(sparse){
            VectorDiff <- rep(NA,length(Vector))
            temp.diff <- rev(diff(rev(Vector[sparsity.idx > sparsity.threshold])))
            temp.diff[length(temp.diff)+1] <- 0
            VectorDiff[sparsity.idx > sparsity.threshold] <- temp.diff
        }else{
            VectorDiff<-rev(diff(rev(Vector)))
            VectorDiff[length(Vector)]=0
        }
        return(VectorDiff)
}
Forwards.Difference <- function(Vector=NULL,sparse=FALSE,sparsity.idx=NULL,sparsity.threshold=NULL){
        if(is.null(Vector)){
            stop("Vector variable cannot be empty.")
        }
        if(sparse){
            VectorDiff <- rep(NA,length(Vector))
            temp.diff <- diff(Vector[sparsity.idx > sparsity.threshold])
            temp.diff <- c(0,temp.diff)
            VectorDiff[sparsity.idx > sparsity.threshold] <- temp.diff
        }else{
        VectorDiff<-diff(Vector)
        VectorDiff=c(0,VectorDiff)
        }
        return(VectorDiff)
}
ComputeOutlierOverIQRInsideWindow <- function(lookup.window=NULL,diff.values=NULL,values=NULL, row.sums = NULL,
        min.sum = NULL, tukeys.constant=NULL,tail=NULL,sparse=FALSE,strict = FALSE,sparsity.idx=NULL,sparsity.threshold=NULL){
        seq.over.value <- seq_along(diff.values)
        Filter <- row.sums > min.sum
        if(sparse){
            Filter <- Filter & sparsity.idx > sparsity.threshold
        }
        seq.over.value <- seq.over.value[Filter]
        seq.over.seq <- seq_along(seq.over.value)
        outlier.list <- lapply(seq.over.seq,function(x.seq){
            lookup.window.range <- ((x.seq - lookup.window) : (x.seq + lookup.window))
            lookup.window.range <- seq.over.value[lookup.window.range[lookup.window.range>0 & lookup.window.range<=max(seq.over.seq)]]
            offset <- (min(lookup.window.range)-1)
            diff.value.window <- diff.values[lookup.window.range]
            value.window <- values[lookup.window.range]
            value.quartile<-quantile(diff.value.window,na.rm=TRUE)
            InterQuartile <- value.quartile[4]-value.quartile[2]
            if(tail=="lower.tail"){
                #Calculate Inner fences based on accepted formula
                fences <- value.quartile[2] - (InterQuartile*tukeys.constant)
                Outlier.Filter <- !is.na(value.window) &
                diff.value.window <= fences &
                diff.value.window < value.window
                if(strict){
                    Outlier.Filter <- Outlier.Filter & (value.window < 0)
                }
            }else if(tail=="upper.tail"){
                fences <- value.quartile[4] + (InterQuartile*tukeys.constant)
                Outlier.Filter <- !is.na(value.window) &
                diff.value.window >= fences &
                diff.value.window > value.window
                if(strict){
                    Outlier.Filter <- Outlier.Filter & (value.window > 0)
                }
            }
            outliers<-which(Outlier.Filter)
            if(length(outliers)>0){
                outliers <- lookup.window.range[outliers]
            }
            outliers
         })
        outlier.vector.dups <- do.call(c,outlier.list)
        outlier.vector.uniq.sorted <- sort(unique(outlier.vector.dups))
        return(outlier.vector.uniq.sorted)
}
CreateDomainlist <- function(start.vector=NULL,end.vector=NULL,fill.gaps=NULL){
    Domains.by.start.list <- lapply(start.vector,function(x){
        data.frame(startbin=x, endbin=min(end.vector[end.vector > x]))
    })
    Domains.by.start.df <- do.call(rbind,Domains.by.start.list)
    Domains.by.start.df$gap.fill <- as.numeric(FALSE)
    Domains.by.start.df$level <- 2
    Domains.by.end.df <- NULL
    Domains.by.assumption.df <- NULL
    if(fill.gaps){
        uncurated.ends <- end.vector[!(end.vector %in% Domains.by.start.df[,"endbin"])]
        if(length(uncurated.ends) > 0){
            Domains.by.end.list <- lapply(uncurated.ends,function(x){
                data.frame(startbin=(max(Domains.by.start.df$endbin[Domains.by.start.df$endbin<x])+1),endbin=x)
            })
            Domains.by.end.df <- do.call(rbind,Domains.by.end.list)
            Domains.by.end.df$gap.fill <- as.numeric(TRUE)
            Domains.by.end.df$level <- 1
        }
    }
    All.Domains <- rbind(Domains.by.start.df,Domains.by.end.df,Domains.by.assumption.df)
    All.Domains.sorted <- All.Domains[order(All.Domains$startbin),]
    return(All.Domains.sorted)
}
ComputeDirectionalityIndex <- function(Matrix = NULL, Window.size=NULL, filter = NULL, start = NULL, end = NULL){
    Sequence <- 1:nrow(Matrix)
    Sequence <- Sequence[start:end]
    DI.Data <- rep(NA,length(Sequence))
    Bins.to.process <- Sequence[filter[start:end]]
    All.bins <- 1:nrow(Matrix)
    All.bins <- All.bins[filter]
    DI.list <- vapply(Bins.to.process,function(i){
            Upstream<-0
            Downstream<-0
            My.DI.Data <- NA
            Relative.mid <- which(All.bins == i)
            Window.range <- c((Relative.mid - Window.size) : (Relative.mid + Window.size))
            Window.range <- All.bins[Window.range[Window.range >= 1 & Window.range <= length(All.bins)]]
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
                # $DI = ( ($B - $A)/abs($B - $A) ) *( (($A - $E)**2)/$E + (($B - $E)**2)/$E);
                My.DI.Data <- ((Downstream - Upstream)/abs(Downstream - Upstream)) * (((Upstream - Expected)^2)/Expected + ((Downstream - Expected)^2)/Expected)
            }
        },1)
    DI.Data[Sequence %in% Bins.to.process] <- DI.list
    return(DI.Data)
}
get_directionality_index_by_chunks <- function(Lego = NULL, chr = NULL, di.window = NULL, distance = NULL,
    chunk.size = 500, sparse = FALSE, sparsity.threshold = 0.8, min.sum = -1, force = FALSE){
    Ranges <- Lego_get_bintable(Lego = Lego, chr = chr)
    First.non.zero.bin <- ._get_first_nonzero_bin(Lego = Lego, chr = chr)
    chr.length <- length(Ranges)
    RowSums <- Lego_get_matrix_mcols(Lego = Lego, chr1 = chr, chr2 = chr, what = "row.sums")
    if(sparse){
        SparsityIndex <- Lego_get_matrix_mcols(Lego = Lego, chr1 = chr, chr2 = chr, what = "sparsity")
    }
    if((chunk.size - (di.window*2))/di.window < 10){
        stop("chunk.size is too small for this di.window\n")
    }
    if(any(di.window > distance)){
        stop("di.window cannot be larger than distance\n")
    }
    Iterations.number <- (chr.length - (First.non.zero.bin - 1))/chunk.size
    Iterations <- rep(chunk.size,floor(Iterations.number))
    if(floor(Iterations.number)!=ceiling(Iterations.number)){
        cumulative <- sum(Iterations)
        Iterations <- c(Iterations,((chr.length - (First.non.zero.bin - 1))-cumulative))
    }
    Starts <- First.non.zero.bin - 1

    if(length(Iterations)>1){
        Starts.cumsum <- cumsum(Iterations)
        Starts <- c(First.non.zero.bin - 1, Starts.cumsum[-length(Starts.cumsum)])
    }
    DI.data.list <- lapply(seq_along(Starts), function(x){
        End <- Starts[x] + Iterations[x]
        Start <- Starts[x] + 1
        Position.start <- Start
        Position.end <- End
        Start <- ifelse((Start - di.window) < First.non.zero.bin, First.non.zero.bin, (Start - di.window))
        End <- ifelse((End + di.window) > chr.length, chr.length, (End + di.window))
        RowSums.subset <- RowSums[Start:End]
        Filter <- RowSums.subset > min.sum
        if(sparse){
            Sparsity.index.subset <- SparsityIndex[Start:End]
            Filter <- Filter & (Sparsity.index.subset > sparsity.threshold)
        }
        Total.length <- length(Filter)
        Filter.extend.length <- length(which(!Filter))
        True.length <- length(which(Filter))
        while(True.length < Total.length){
            if(i == 1){
                extend <- Filter.extend.length
            }
            Start <- ifelse((Start - extend) < First.non.zero.bin, First.non.zero.bin, Start - extend)
            End <- ifelse((End + extend) > chr.length, chr.length, (End + extend))
            RowSums.subset <- RowSums[Start:End]
            Filter <- RowSums.subset > min.sum
            if(sparse){
                Sparsity.index.subset <- SparsityIndex[Start:End]
                Filter <- Filter & (Sparsity.index.subset > sparsity.threshold)
            }
            True.length <- length(which(Filter))
            extend <- extend + 1
        }
        Matrix <- Lego_get_vector_values(Lego = Lego, chr1 = chr, chr2 = chr, xaxis=c(Start:End), yaxis=c(Start:End), force = force)
        DI.data <- ComputeDirectionalityIndex(Matrix = Matrix, Window.size = di.window, filter = Filter, 
            start = Position.start - (Start - 1), end = Position.end - (Start - 1))
        return(DI.data)
    })
    DI.data <- do.call(c, DI.data.list)
    DI.data <- c(rep(NA,First.non.zero.bin - 1),DI.data)
    Ranges$DI.Data <- DI.data
    return(Ranges)
}
MakeBoundaries <- function(chr = NULL, Ranges = NULL, Binsize = NULL){
    Ends <- end(Ranges) 
    Ends <- Ends[1:(length(Ends)-1)]
    Starts <- start(Ranges)
    Starts <- Starts[2:length(Starts)] - 1 
    Domain.boundaries <- unique(Starts,Ends)
    Boundary.ranges <- MakeGRangesObject(chr=rep(chr,length(Domain.boundaries)),
        Start=(Domain.boundaries-(Binsize/2))+1,End=Domain.boundaries+(Binsize/2))
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
#' tukeys.constant x inter-quartile range of the directionality index. The upper
#' fence used for detecting domain starts is the 75th quartile + 
#' (IQR x tukeys.constant), while the lower fence is the 
#' 25th quartile - (IQR x tukeys.constant). For domain starts the DI difference
#' must be greater than or equal to the upper fence, it must be greater than the
#' DI and the DI must be a finite real value. If strict is TRUE, DI will also
#' be required to be greater than 0. Similarly, for domain ends the 
#' DI difference must be lower than or equal to the lower fence, it must be 
#' lower than the DI and the DI must be a finite real value. If strict is TRUE,
#' DI will also be required to be lower than 0. 
#' 
#' After defining outliers, each domain start will be associated to its 
#' nearest downstream domain end. If \emph{fill.gaps} is defined as TRUE and
#' there are domain ends which remain unassociated to a domain start, These 
#' domain ends will be associated to the bin adjacent to their nearest upstream
#' domain end. This associations will be marked by metadata columns, gap.fill= 1
#' and level = 1.
#' 
#' This function provides the capability to call very accurante TAD definitions
#' in a very fast way. 
#' 
#' @inheritParams Lego_get_chrominfo
#' 
#' @param chrs \strong{Optional}. Default NULL
#' If present, only TAD calls for elements in \emph{chrs} will be done.
#' 
#' @param min.sum \strong{Optional}. Default -1
#' Process bins in the matrix with row.sums greater than \emph{min.sum}.
#' 
#' @param di.window \strong{Optional}. Default 200
#' Use \emph{di.window} to define the directionality index.
#' 
#' @param lookup.window \strong{Optional}. Default 200
#' Use \emph{lookup.window} local window to call borders. At smaller 
#' \emph{di.window} values we recommend setting this to 2*\emph{di.window}
#' 
#' @param tukeys.constant \strong{Optional}. Default 1.5
#' \emph{tukeys.constant}*IQR (inter-quartile range) defines the lower and upper
#' fence values.
#' 
#' @param strict \strong{Optional}. Default TRUE
#' If TRUE, \emph{strict} creates an additional filter on the directionality 
#' index requiring it to be either greater than or less than 0 on the right tail
#' or left tail respectively.  
#' 
#' @param fill.gaps \strong{Optional}. Default TRUE
#' If TRUE, this will affect the TAD stiching process. All Border starts are 
#' stiched to the next downstream border ends. Therefore, at times border ends 
#' remain unassociated to a border start. These border ends are stiched to the 
#' adjacent downstream bin from their upstream border end when \emph{fill.gaps} 
#' is true. 
#' 
#' TADs inferred in this way will be annotated with two metadata columns in the 
#' GRanges object. \emph{gap.fill} will hold a value of 1 and \emph{level} will 
#' hold a value 1. TADs which were not filled in will hold a gap.fill value of 0 
#' and a level value of 2.
#' 
#' @param ignore.sparse \strong{Optional}. Default TRUE
#' If TRUE, a matrix which has been defined as sparse during the matrix loading
#' process will be treated as a dense matrix. The \emph{sparsity.threshold} 
#' filter will not be applied. Please note, that if a matrix is defined as 
#' sparse and fill.gaps is TRUE, fill.gaps will be turned off.
#' 
#' @param sparsity.threshold \strong{Optional}. Default 0.8
#' Sparsity threshold relates to the sparsity index, which is computed as the 
#' number of non-zero bins at a certain distance from the diagonal. If a matrix
#' is sparse and ignore.sparse is FALSE, bins which have a sparsity index value
#' below this threshold will be discarded from DI computation.
#' 
#' @param remove.empty Not implemented.
#' After implementation, this will ensure that the presence of centromeric 
#' regions is accounted for.
#' 
#' @param chunk.size \strong{Optional}. Default 500
#' The size of the matrix chunk to process. This value should be larger than 2x
#' di.window.
#' 
#' @param force.retrieve \strong{Optional}. Default TRUE
#' If TRUE, this will force the retrieval of a matrix chunk even when the 
#' retrieval includes interaction points which were not loaded into a Lego 
#' store (larger chunks). Please note, that this does not mean that DI can be 
#' computed at distances larger than max distance. Rather, this is meant to aid
#' faster computation.
#' 
#' @return A ranges object containing domain definitions. The starts and ends
#' of the ranges coincide with the starts and ends of their contained bins from 
#' the bintable. 
#' 
Lego_local_score_differentiator <- function(Lego = NULL, chrs = NULL, min.sum = -1, 
    di.window = 200L, lookup.window = 200L, tukeys.constant=1.5, strict = TRUE, 
    fill.gaps=TRUE, ignore.sparse=TRUE, sparsity.threshold=0.8,
    remove.empty = NULL, chunk.size = 500, force.retrieve = TRUE){

    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    Chromosomes <- ChromInfo[,'chr']
    if(!is.null(chrs)){
        Chromosomes <- ChromInfo[ChromInfo[,'chr'] %in% chrs,'chr']
    }

    Chrom.domains.ranges.list <- lapply(Chromosomes, function(chr){
        Ranges <- Lego_get_bintable(Lego = Lego, chr = chr)
        sparse <- Lego_matrix_issparse(Lego = Lego, chr1 = chr, chr2 = chr)
        max.distance <- Lego_matrix_maxdist(Lego = Lego, chr1 = chr, chr2 = chr)
        if(ignore.sparse){
            sparse=FALSE
        }
        if(sparse & fill.gaps){
            fill.gaps=FALSE
        }
        cat("[1] Computing DI for",chr,"\n")
        Ranges <- get_directionality_index_by_chunks(Lego = Lego, chr = chr, di.window = di.window, 
            distance = max.distance, chunk.size = chunk.size, sparse=sparse, sparsity.threshold=sparsity.threshold,
            min.sum = min.sum, force = force.retrieve)

        RowSums <- Lego_get_matrix_mcols(Lego = Lego, chr1 = chr, chr2 = chr, what = "row.sums")
        Ranges$row.sums <- RowSums
        cat("[2] Computing DI Differences for",chr,"\n")
        if(sparse){
            SparsityIndex <- Lego_get_matrix_mcols(Lego = Lego, chr1 = chr, chr2 = chr, what = "sparsity")
            Backwards.DI.Difference <- Backwards.Difference(Vector = Ranges$DI.Data, sparse = sparse,
                sparsity.idx = SparsityIndex, sparsity.threshold = sparsity.threshold)
            Forwards.DI.Difference <- Forwards.Difference(Vector = Ranges$DI.Data, sparse = sparse,
                sparsity.idx = SparsityIndex, sparsity.threshold = sparsity.threshold)
        }else{
            Backwards.DI.Difference <- Backwards.Difference(Vector=Ranges$DI.Data)
            Forwards.DI.Difference <- Forwards.Difference(Vector=Ranges$DI.Data)
        }
        Ranges$backward.Differences <- Backwards.DI.Difference
        Ranges$forward.Differences <- Forwards.DI.Difference
        cat("[2] Done\n")
        cat("[3] Fetching Outliers ",chr,"\n")
        if(sparse){
            Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(lookup.window=lookup.window,
                diff.values=Backwards.DI.Difference, values=Ranges$DI.Data, sparse=sparse, row.sums = Ranges$row.sums,
                min.sum = min.sum, sparsity.idx=SparsityIndex, sparsity.threshold=sparsity.threshold, 
                tukeys.constant=tukeys.constant, tail="lower.tail",strict=strict)
            Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(lookup.window=lookup.window,
                diff.values=Forwards.DI.Difference, values=Ranges$DI.Data, sparse=sparse, row.sums = Ranges$row.sums,
                min.sum = min.sum, sparsity.idx=SparsityIndex, sparsity.threshold=sparsity.threshold,
                tukeys.constant=tukeys.constant, tail="upper.tail", strict=strict)
        }else{
            Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(lookup.window=lookup.window,
                diff.values=Backwards.DI.Difference,values=Ranges$DI.Data, row.sums = Ranges$row.sums,
                min.sum = min.sum, tukeys.constant=tukeys.constant,tail="lower.tail",strict=strict)
            Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(lookup.window=lookup.window,
                diff.values=Forwards.DI.Difference,values=Ranges$DI.Data, row.sums = Ranges$row.sums,
                min.sum = min.sum, tukeys.constant=tukeys.constant,tail="upper.tail",strict=strict)
        }
        Domain.start.candidates <- Domain.start.candidates[Domain.start.candidates != length(Ranges)]
        Domain.end.candidates <- Domain.end.candidates[Domain.end.candidates != 1]
        cat("[3] Done\n")
        cat("[4] Creating Domain list for",chr,"\n")

        if(!(1 %in% Domain.start.candidates)){
            Domain.start.candidates <- c(1,Domain.start.candidates)
        }
        if(!(length(Ranges) %in% Domain.end.candidates)){
            Domain.end.candidates <- c(Domain.end.candidates,length(Ranges))
        }
        Domain.list <- CreateDomainlist(start.vector=Domain.start.candidates,
            end.vector=Domain.end.candidates,fill.gaps=fill.gaps)
        Domain.Ranges <- Lego_make_ranges(Chrom=rep(chr,nrow(Domain.list)),
            Start=start(Ranges[Domain.list$startbin]),
            End=end(Ranges[Domain.list$endbin]))
        cat("[4] Done\n")
        Domain.Ranges$gap.fill <- Domain.list$gap.fill
        Domain.Ranges$level <- Domain.list$level
        Domain.Ranges$window.size <- di.window
        Domain.Ranges$lookup.window <- lookup.window
        return(Domain.Ranges)
    })
    Chrom.domains.ranges <- do.call(c,unlist(Chrom.domains.ranges.list, use.names = TRUE))
    return(Chrom.domains.ranges)
}