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
ComputeOutlierOverIQRInsideWindow <- function(Lookup.window=NULL,diff.values=NULL,values=NULL, row.sums = NULL,
        min.sum = NULL, Tukeys.constant=NULL,tail=NULL,sparse=FALSE,strict = FALSE,sparsity.idx=NULL,sparsity.threshold=NULL){
        seq.over.value <- seq_along(diff.values)
        Filter <- row.sums > min.sum
        if(sparse){
            Filter <- Filter & sparsity.idx > sparsity.threshold
        }
        seq.over.value <- seq.over.value[Filter]
        seq.over.seq <- seq_along(seq.over.value)
        outlier.list <- lapply(seq.over.seq,function(x.seq){
            lookup.window.range <- ((x.seq - Lookup.window) : (x.seq + Lookup.window))
            lookup.window.range <- seq.over.value[lookup.window.range[lookup.window.range>0 & lookup.window.range<=max(seq.over.seq)]]
            offset <- (min(lookup.window.range)-1)
            diff.value.window <- diff.values[lookup.window.range]
            value.window <- values[lookup.window.range]
            value.quartile<-quantile(diff.value.window,na.rm=TRUE)
            InterQuartile <- value.quartile[4]-value.quartile[2]
            if(tail=="lower.tail"){
                #Calculate Inner fences based on accepted formula
                fences <- value.quartile[2] - (InterQuartile*Tukeys.constant)
                Outlier.Filter <- !is.na(value.window) &
                diff.value.window <= fences &
                diff.value.window < value.window
                if(strict){
                    Outlier.Filter <- Outlier.Filter & (value.window < 0)
                }
            }else if(tail=="upper.tail"){
                fences <- value.quartile[4] + (InterQuartile*Tukeys.constant)
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
CreateDomainlist <- function(start.vector=NULL,end.vector=NULL,Fill.gaps=NULL){
    Domains.by.start.list <- lapply(start.vector,function(x){
        data.frame(startbin=x, endbin=min(end.vector[end.vector > x]))
    })
    Domains.by.start.df <- do.call(rbind,Domains.by.start.list)
    Domains.by.start.df$gap.fill <- as.numeric(FALSE)
    Domains.by.start.df$level <- 2
    Domains.by.end.df <- NULL
    Domains.by.assumption.df <- NULL
    if(Fill.gaps){
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
    Sequence <- start:end
    Sequence <- Sequence[filter[start:end]]
    Seq.Over.Sequence <- seq_along(Sequence)
    DI.Data <- rep(NA,length(Sequence))
    DI.list <- vapply(Seq.Over.Sequence,function(seq.i){
            Upstream<-0
            Downstream<-0
            My.DI.Data <- NA
            i <- Sequence[seq.i]
            Window.range <- c((seq.i-Window.size):(seq.i+Window.size))
            Window.range <- Sequence[Window.range]
            Relative.mid <- which(Window.range==i)
            Upstream.range <- Window.range[Window.range < i]
            Downstream.range <- Window.range[Window.range > i]
            if(i >= FirstNonZeroBin)
            {
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
            }else{
                NA
            }
        })
    DI.Data[Sequence] <- do.call(c,DI.list)
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
    if(any(di.window > distance){
        stop("di.window cannot be larger than distance\n")
    }
    Iterations.number <- (chr.length - (Position.start - 1))/chunk.size
    Iterations <- rep(chunk.size,floor(Iterations.number))
    if(floor(Iterations.number)!=ceiling(Iterations.number)){
        cumulative <- sum(Iterations)
        Iterations <- c(Iterations,((chr.length - (Position.start - 1))-cumulative))
    }
    Starts<-Position.start
    if(length(Iterations)>1){
        Starts.cumsum <- cumsum(Iterations)
        Starts <- c(Position.start,Starts.cumsum[-length(Starts.cumsum)])
    }
    DI.data.list <- lapply(seq_along(Starts), function(x){
        Start <- Starts[x] + 1
        End <- Start + Iterations[x]
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
        Matrix <- Lego_get_vector_values(Lego = Lego, chr1 = chr, chr2 = chr, xaxis=c(Start,End), yaxis=c(Start,End), force = force)
        DI.data <- ComputeDirectionalityIndex(Matrix = Matrix, Window.size = di.window, filter = Filter, 
            start = Position.start - (Start - 1), end = Position.end - (Start - 1))
        return(DI.data)
    })
    DI.data <- do.call(c, DI.data.list)
    DI.data <- c(rep(NA,First.non.zero.bin - 1),DI.data)
    Ranges$DI.Data <- DI.data
    return(Ranges)
}
LocalScoreDifferentiator <- function(Lego = NULL, chrs = NULL, min.sum = -1, di.window = 200L, lookup.window = 200L,
    tukeys.constant=1.5, strict = TRUE, fill.gaps=TRUE, ignore.sparse=TRUE, sparsity.threshold=0.8, 
    remove.empty = NULL, chunk.size = 500, force.retrieve = FALSE){

    ChromInfo <- Lego_get_chrominfo(Lego = Lego)
    Chromosomes <- ChromInfo[,'chr']
    if(!is.null(chrs))}{
        Chromosomes <- ChromInfo[ChromInfo[,'chr'] %in% chrs,'chr']
    }

    Chrom.domains.ranges.list <- lapply(Chromosomes, function(chr){
        Ranges <- Lego_get_bintable(Lego = Lego, chr = chr)
        sparse <- Lego_matrix_issparse(Lego = Lego, chr1 = chr, chr2 = chr)
        max.distance <- Lego_matrix_maxdist(Lego = Lego, chr1 = chr, chr2 = chr)
        if(ignore.sparse){
            sparse=FALSE
        }
        if(sparse & Fill.gaps){
            Fill.gaps=FALSE
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
            Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
                diff.values=Backwards.DI.Difference, values=Ranges$DI.Data, sparse=sparse, row.sums = Ranges$row.sums,
                min.sum = min.sum, sparsity.idx=SparsityIndex, sparsity.threshold=sparsity.threshold, 
                Tukeys.constant=Tukeys.constant, tail="lower.tail",strict=strict)
            Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
                diff.values=Forwards.DI.Difference, values=Ranges$DI.Data, sparse=sparse, row.sums = Ranges$row.sums,
                min.sum = min.sum, sparsity.idx=SparsityIndex, sparsity.threshold=sparsity.threshold,
                Tukeys.constant=Tukeys.constant, tail="upper.tail", strict=strict)
        }else{
            Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
                diff.values=Backwards.DI.Difference,values=Ranges$DI.Data, row.sums = Ranges$row.sums,
                min.sum = min.sum, Tukeys.constant=Tukeys.constant,tail="lower.tail",strict=strict)
            Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
                diff.values=Forwards.DI.Difference,values=Ranges$DI.Data, row.sums = Ranges$row.sums,
                min.sum = min.sum, Tukeys.constant=Tukeys.constant,tail="upper.tail",strict=strict)
        }
        Domain.start.candidates <- Domain.start.candidates[Domain.start.candidates != length(Ranges)]
        Domain.end.candidates <- Domain.end.candidates[Domain.end.candidates != 1]
        cat("[3] Done\n")
        cat("[4] Creating Domain list for",chr.name,"\n")

        if(!(1 %in% Domain.start.candidates)){
            Domain.start.candidates <- c(1,Domain.start.candidates)
        }
        if(!(length(Ranges) %in% Domain.end.candidates)){
            Domain.end.candidates <- c(Domain.end.candidates,length(Ranges))
        }
        Domain.list <- CreateDomainlist(start.vector=Domain.start.candidates,
            end.vector=Domain.end.candidates,Fill.gaps=Fill.gaps)
        Domain.Ranges <- MakeGRangesObject(chr=rep(chr.name,nrow(Domain.list)),
            Start=start(Ranges[Domain.list$startbin]),
            End=end(Ranges[Domain.list$endbin]))
        cat("[4] Done\n")
        Domain.Ranges$gap.fill <- Domain.list$gap.fill
        Domain.Ranges$level <- Domain.list$level
        Domain.Ranges$window.size <- DI.Window.size
        Domain.Ranges$lookup.window <- Lookup.window
        return(Domain.Ranges)
    })
    Chrom.domains.ranges <- do.call(c,unlist(Chrom.domains.ranges.list, use.names = TRUE))
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