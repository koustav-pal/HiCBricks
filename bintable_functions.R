_ColIndexError_ <- function(x){
    Error.col <- c("chr","start",NA,NA,NA,"more")
    ColClasses <- c("character","numeric","numeric","character","character")
    ColNames <- c("chr","start","end","strand","names")
    Index.type.error <- c("index is missing.")
    if(x < 3 | x > 5){
        if(x > 5){
            Index.type.error <- "indices were provided."
        }
        stop(paste("col.index expects as bare minimum chr,start,end.",Error.col[x],Second.part,"\n"))
    }
    Alist <- list("Names" = ColNames[1:x], "Classes" = ColClasses[1:x])
    return(Alist)
}


Read_bintable = function(Filename=NULL,read.delim=" ",exec="cat", col.index=c(1,2,3), 
    chromosomes=NULL, impose.discontinuity=TRUE){
    require(data.table)
    if(is.null(RangeKey)) {
        stop("Variable name not provided while reading table.")
    }
    if(is.null(exec)) {
        stop("exec is not allowed to be null")
    }
    ColMetrics <- _ColIndexError_(length(col.index))
    Colnames<-ColMetrics[["Names"]]
    ColClasses<- ColMetrics[["Classes"]]
    
    Table <- Filename
    if(class(Filename)=="character"){
        Command <- paste(exec,Filename,sep=" ")
        Table <- fread(input=Command, sep=read.delim,
            stringsAsFactors=FALSE, verbose=FALSE, showProgress=FALSE, data.table=FALSE)
    }
    Ranges.table <- Table[,col.index]
    is.stranded <- col.index[4] != NA
    has.names <- col.index[5] != NA
    Validate_table(Table=Table,stranded=is.stranded,colClasses=ColClasses,named=has.names,chrom=chromosomes)
    if(impose.discontinuity){
        CheckContinuousRanges(Table=Ranges.table,StartCol=c("start"),EndCol=c("end"))
    }
    Ranges.table <- Ranges.table[order(Ranges.table[,'chr'],Ranges.table[,'start']),]
    Table.list <- list('main.tab' = ValidatedTable, 'stranded' = is.stranded, 'named' = has.names)
    return(Table.list)
}
Validate_table = function(Table=NULL,colnames=NULL,colClasses=NULL,col.index=NULL,chrom=NULL) {
    for (i in 1:length(colnames)) {
        if(class(Table[,i])!=colClasses[i]){
            stop(paste(colClasses[i],"values expected for",colnames[i],"at col",col.index[i]))
        }
    }
    UniqueChromNames<-unique(Table[,'chr'])
    if(any( !(UniqueChromNames %in% chrom) )){
        stop("Some chromosome names are not defined in the chromosome table")
    }
    if(any(!(Table[,'start'] %% 1 == 0)) | any(!(Table[,'end'] %% 1 == 0))) {
        stop("Genomic coordinates at col,"col.index[2],"and",col.index[3],"cannot have float values")
    }
    if( any( Table[,'start'] > Table[,'end'] ) ){
        stop("start coordinates cannot be greater than end coordinates")
    }
}

CheckContinuousRanges = function(Table=NULL, StartCol=NULL, EndCol=NULL){
    Starts<-Table[,StartCol]
    Starts<-Starts[2:length(Starts)]
    End<-Table[,EndCol]
    End<-End[1:(length(End)-1)]
    if( any(Starts==End) ){
        stop("Found continuous ranges in file! Cannot proceed further!
            Use impose.discontinuity = FALSE to load continuous ranges.")
    }
}

get_chrom_info <- function(bin.table = NULL, chrom = NULL, FUN = NULL, col.name = NULL){
    Info <- sapple(chrom,function(x){
        FUN(bin.table[,col.name][bin.table[,'chr']==x])
    }) 
    names(Info) <- chrom
    return(Info)
}