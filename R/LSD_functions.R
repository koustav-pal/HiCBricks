GetFirstNonZeroBin <- function(Matrix=NULL){
	min(which(rowSums(Matrix) > 0))
}
Return.sparsity.index <- function(Matrix = NULL,sparsity.bins = 100){
    SparsityIndex = function(x=NULL){
        x[is.na(x) | is.infinite(x)] <- 0
        return(length(x[x!=0])/length(x))
    }
    Sparsity.Index <- sapply(1:nrow(Matrix),function(x){
        Range <- (x - sparsity.bins) : (x + sparsity.bins)
        Range <- Range[Range>0 & Range<nrow(Matrix)]
        Row.vector <- Matrix[x,Range]
        SparsityIndex(Row.vector)
    })
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
ComputeOutlierOverIQRInsideWindow <- function(Lookup.window=NULL,diff.values=NULL,values=NULL,
        Tukeys.constant=NULL,tail=NULL,sparse=FALSE,strict = FALSE,sparsity.idx=NULL,sparsity.threshold=NULL){
        seq.over.value <- seq_along(diff.values)
        if(sparse){
            seq.over.value <- seq.over.value[sparsity.idx > sparsity.threshold]
        }
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
ComputeDirectionalityIndex <- function(Matrix=NULL, Window.size=NULL, Lookup.window = NULL, 
	sparse=FALSE,sparsity.threshold=0.8,filter.chromosomes=NULL){
    FirstNonZeroBin <- GetFirstNonZeroBin(Matrix=Matrix)
    Sequence <- 1:nrow(Matrix)
    if(sparse){
    	sparsity.idx <- Return.sparsity.index(Matrix=Matrix)
        Sequence <- Sequence[sparsity.idx > sparsity.threshold]
    }
    Seq.Over.Sequence <- seq_along(Sequence)
    DI.Data <- rep(NA,nrow(Matrix))
    DI.list <- mclapply(Seq.Over.Sequence,function(seq.i){
            Upstream<-0
            Downstream<-0
            My.DI.Data <- NA
            i <- Sequence[seq.i]
            Window.range <- c((seq.i-Window.size):(seq.i+Window.size))
            Window.range <- Sequence[Window.range[Window.range>=1 & Window.range<=max(Seq.Over.Sequence)]]
            Window.range <- Window.range[Window.range >= FirstNonZeroBin & Window.range <= nrow(Matrix)]
            Relative.mid <- which(Window.range==i)
            Upstream.range <- Window.range[Window.range < i]
            Downstream.range <- Window.range[Window.range > i]
            if(i >= FirstNonZeroBin)
            {
                # Row.vector <- unlist(private$HiCData.Obj$FetchRowOrColVector(Chrom1=Chrom,Chrom2=Chrom,by="position",vector=i),use.names=FALSE)
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
        },mc.cores = 32)
    DI.Data[Sequence] <- do.call(c,DI.list)
    return(DI.Data)
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
LocalScoreDifferentiator <- function(Chrom = NULL, Ranges = NULL, Matrix=NULL, DI.Window.size=200L,Lookup.window=200L,
    Tukeys.constant=1.5,strict = TRUE,Fill.gaps=TRUE,ignore.sparse=TRUE,sparsity.threshold=0.8){
    ChromName <- Chrom
	sparse <- TRUE
    if(ignore.sparse){
        sparse=FALSE
    }
    if(sparse & Fill.gaps){
        Fill.gaps=FALSE
    	Ranges$sparsity.idx <- Return.sparsity.index(Matrix=Matrix)
    }
    cat("[1] Computing DI for",ChromName,"\n")
    DI.Data <- ComputeDirectionalityIndex(Matrix = Matrix, Window.size=DI.Window.size, 
    	sparse=sparse, sparsity.threshold=sparsity.threshold)
    Ranges$DI.Data <- DI.Data
    # return(Ranges)
    cat("[2] Computing DI Differences for",ChromName,"\n")

    if(sparse){
        Backwards.DI.Difference <- Backwards.Difference(Vector=Ranges$DI.Data, sparse=sparse,
            sparsity.idx=Ranges$sparsity.idx, sparsity.threshold=sparsity.threshold)
        Forwards.DI.Difference <- Forwards.Difference(Vector=Ranges$DI.Data, sparse=sparse,
            sparsity.idx=Ranges$sparsity.idx, sparsity.threshold=sparsity.threshold)
    }else{
        Backwards.DI.Difference <- Backwards.Difference(Vector=Ranges$DI.Data)
        Forwards.DI.Difference <- Forwards.Difference(Vector=Ranges$DI.Data)
    }
    Ranges$backward.Differences <- Backwards.DI.Difference
    Ranges$forward.Differences <- Forwards.DI.Difference
    cat("[2] Done\n")
    cat("[3] Fetching Outliers ",ChromName,"\n")
    if(sparse){
        Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
            diff.values=Backwards.DI.Difference,values=Ranges$DI.Data,sparse=sparse,
            sparsity.idx=Ranges$sparsity.idx,sparsity.threshold=sparsity.threshold,
            Tukeys.constant=Tukeys.constant,tail="lower.tail",strict=strict)
        Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
            diff.values=Forwards.DI.Difference,values=Ranges$DI.Data,sparse=sparse,
            sparsity.idx=Ranges$sparsity.idx,sparsity.threshold=sparsity.threshold,
            Tukeys.constant=Tukeys.constant,tail="upper.tail",strict=strict)
    }else{
        Domain.end.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
            diff.values=Backwards.DI.Difference,values=Ranges$DI.Data,
            Tukeys.constant=Tukeys.constant,tail="lower.tail",strict=strict)
        Domain.start.candidates <- ComputeOutlierOverIQRInsideWindow(Lookup.window=Lookup.window,
            diff.values=Forwards.DI.Difference,values=Ranges$DI.Data,
            Tukeys.constant=Tukeys.constant,tail="upper.tail",strict=strict)
    }
    Domain.start.candidates <- Domain.start.candidates[Domain.start.candidates!=length(Ranges)]

    Domain.end.candidates <- Domain.end.candidates[Domain.end.candidates!=1]
    cat("[3] Done\n")
    cat("[4] Creating Domain list for",ChromName,"\n")
    if(!(1 %in% Domain.start.candidates)){
        Domain.start.candidates <- c(1,Domain.start.candidates)
    }
    if(!(length(Ranges) %in% Domain.end.candidates)){
        Domain.end.candidates <- c(Domain.end.candidates,length(Ranges))
    }
    Domain.list <- CreateDomainlist(start.vector=Domain.start.candidates,
        end.vector=Domain.end.candidates,Fill.gaps=Fill.gaps)
    Domain.Ranges <- MakeGRangesObject(Chrom=rep(ChromName,nrow(Domain.list)),
        Start=start(Ranges[Domain.list$startbin]),
        End=end(Ranges[Domain.list$endbin]))
    cat("[4] Done\n")
    Domain.Ranges$gap.fill <- Domain.list$gap.fill
    Domain.Ranges$level <- Domain.list$level
    Domain.Ranges$window.size <- DI.Window.size
    Domain.Ranges$lookup.window <- Lookup.window
    return(Domain.Ranges)
}
MakeBoundaries <- function(Chrom = NULL, Ranges = NULL, Binsize = NULL){
	Ends <- end(Ranges) 
	Ends <- Ends[1:(length(Ends)-1)]
	Starts <- start(Ranges)
	Starts <- Starts[2:length(Starts)] - 1 
	Domain.boundaries <- unique(Starts,Ends)
	Boundary.ranges <- MakeGRangesObject(Chrom=rep(Chrom,length(Domain.boundaries)),
		Start=(Domain.boundaries-(Binsize/2))+1,End=Domain.boundaries+(Binsize/2))
}