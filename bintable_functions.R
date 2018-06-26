Read_bintable = function(Filename=NULL,read.delim=" ",exec="cat",
    col.index=c(1,2,3),strand=NULL,names=NULL,chromosomes=NULL,
    impose.discontinuity=TRUE){
    is.stranded = FALSE
    has.names = FALSE
    require(data.table)
    if(is.null(RangeKey)) {
        stop("Variable name not provided while reading table.")
    }
    if(is.null(exec)) {
        stop("exec is not allowed to be null")
    }
    if(length(col.index) != 3){
        stop("col.index must be of length 3, consisting of chr, start, end. 
            strand and names can be specified using specific args.")
    }
    Colnames<-c('chr','start','end')
    colClasses=c("character","numeric","numeric")
    Command <- paste(exec,Filename,sep=" ")
    Table <- fread(input=Command, sep=read.delim,
        stringsAsFactors=FALSE, verbose=FALSE, showProgress=FALSE, data.table=FALSE)

    if(!is.null(strand)){
        is.stranded <- TRUE
        col.index <- c(col.index,strand)
        Colnames <- c(Colnames,'strand')
        colClasses <- c(colClasses,'character')
    }
    if(!is.null(names)){
        has.names <- TRUE
        col.index <- c(col.index,names)
        Colnames <- c(Colnames,'names')
        colClasses <- c(colClasses,'character')
    }

    Ranges.table <- Table[,col.index]
    
    Validate_table(Table=Table,stranded=is.stranded,named=has.names,chrom=chromosomes)
    
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
    names(Sizes) <- chrom
    return(Sizes)
}