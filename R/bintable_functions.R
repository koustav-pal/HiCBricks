._ColIndexError_ <- function(x){
    Error.col <- c("chr","start",NA,NA,NA,"more")
    ColClasses <- list(is.character,is.numeric,is.numeric,
        is.character,is.character)
    ColNames <- c("chr","start","end","strand","names")
    Index.type.error <- c("index is missing.")
    Norm.x <- (x - min(x))+1
    Len.x <- length(Norm.x)

    if(Len.x < 3 | Len.x > 5){
        if(Len.x > 5){
            Index.type.error <- "indices were provided."
        }
        stop(paste("col.index expects as bare minimum chr,start,end.",
            Error.col[x],Index.type.error,"\n"))
    }
    Alist <- list("Names" = ColNames[Norm.x], "Classes" = ColClasses[Norm.x])
    return(Alist)
}
Read_bintable = function(Filename = NULL, read.delim = " ", 
    exec = "cat", 
    col.index = c(1,2,3), 
    chromosomes=NULL, impose.discontinuity=TRUE){
    if(is.null(exec)) {
        stop("exec is not allowed to be null")
    }
    ColMetrics <- ._ColIndexError_(col.index)
    Colnames<-ColMetrics[["Names"]]
    ColClasses<- ColMetrics[["Classes"]]
    
    Table <- Filename
    if(is.character(Filename)){
        Command <- paste(exec,Filename,sep=" ")
        Table <- fread(input=Command, sep=read.delim,
            stringsAsFactors=FALSE, verbose=FALSE, showProgress=FALSE, 
            data.table=FALSE)
    }
    colnames(Table) <- Colnames
    Ranges.table <- Table[,col.index]
    is.stranded <- col.index[4] != NA
    has.names <- col.index[5] != NA
    Validate_table(Table=Table, colClasses=ColClasses, colnames = Colnames, 
        col.index=col.index, chrom=chromosomes)
    if(impose.discontinuity){
        CheckContinuousRanges(Table=Table,StartCol=c("start"),EndCol=c("end"))
    }
    Ranges.table <- Table[order(Table[,'chr'],Table[,'start']),]
    Table.list <- list('main.tab' = Ranges.table, 'stranded' = is.stranded, 
        'named' = has.names)
    return(Table.list)
}
Validate_table = function(Table = NULL, colnames = NULL, colClasses = NULL, 
    col.index = NULL, chrom = NULL) {
    for (i in seq_len(length(colnames))) {
        if(!colClasses[[i]](Table[,i])){
            stop(paste("Values expected for",colnames[i],"at col",col.index[i],"
                found values of class",class(Table[,i])))
        }
    }
    UniqueChromNames<-unique(Table[,'chr'])
    if(any( !(UniqueChromNames %in% chrom) )){
        stop("Some chromosome names are not defined in the chromosome table")
    }
    if(any( !(chrom %in% UniqueChromNames) )){
        stop("Some chromosome names are not defined in the chromosome table")
    }
    if(any(!(Table[,'start'] %% 1 == 0)) | any(!(Table[,'end'] %% 1 == 0))) {
        stop("Genomic coordinates at col",col.index[2],"and",col.index[3],
            "cannot have float values")
    }
    if( any( Table[,'start'] > Table[,'end'] ) ){
        stop("start coordinates cannot be greater than end coordinates")
    }
    # if( is.unsorted(Table[,'chr']) ){
    #     stop("Table must be sorted by chromosome!")
    # }
}
CheckContinuousRanges = function(Table = NULL, StartCol = NULL, EndCol = NULL){
    Starts<-Table[,StartCol]
    Starts<-Starts[seq_len(length(Starts))[-1]]
    End<-Table[,EndCol]
    End<-End[seq_len(length(End)-1)]
    if( any(Starts==End) ){
        stop("Found continuous ranges in file! Cannot proceed further!
            Use impose.discontinuity = FALSE to load continuous ranges.")
    }
}
get_chrom_info <- function(bin.table = NULL, chrom = NULL, FUN = NULL, 
    col.name = NULL){
    if(is.null(chrom)){
        chrom <- unique(bin.table[,"chr"])
    }
    Info <- sapply(chrom,function(x){
        FUN(bin.table[bin.table[,'chr']==x,col.name])
    })
    names(Info) <- chrom
    return(Info)
}
Split_genomic_coordinates = function(Coordinate = NULL){
    Reference.object <- GenomicMatrix$new()
    Sep <- Reference.object$Ranges.separator
    Coord.Split<-stringr::str_split(pattern=Sep,string=Coordinate)
    if(length(Coord.Split[[1]])!=3 | length(Coord.Split[[1]])!=3){
        stop("Coordinate must be separated by :")
    }
    Chrom<-Coord.Split[[1]][1]
    start<-as.numeric(Coord.Split[[1]][2])
    stop<-as.numeric(Coord.Split[[1]][3])
    if(any(class(Chrom)!="character" | 
        class(start)!="numeric" | 
        class(stop)!="numeric")){
        stop("Provided chromosome,start,end do not match expected ",
            "class definitions of character, numeric, numeric")
    }
    return(Coord.Split)
}