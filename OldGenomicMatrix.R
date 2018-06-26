require(R6)
GenomicMatrix <- R6Class("GenomicMatrix",
    public = list (
        ChromosomeList=NA,
        initialize = function(){
            cat("Initialized all variable names")
        },
        LoadRangesObject = function(RangesObject=NULL,RangeKey=NULL){
            require(GenomicRanges)
            if(any(is.null(c(RangesObject,RangeKey)))) {
                stop("RangesObject or RangeKey cannot be null..")
            }
            if(class(RangesObject) != "GRanges") {
                stop("RangesObject must be of class GenomicRanges")
            }
            if(any(!(GenomeInfoDb::seqnames(RangesObject) %in% self$ChromosomeList))) {
                stop("RangesObject chromosome names do not match")
            }
            if(length(RangesObject)==0 | length(RangeKey) > 1) {
                stop("RangesObject cannot be of length 0. RangeKey length cannot be greater than 1")
            }
            if((RangeKey %in% private$Protected.Ranges.Keys)) {
                stop("RangeKey cannot have a protected key value.",
                    paste(private$Protected.Ranges.Keys,collapse=" "))
            }
            private$RangesObjects[[RangeKey]] <- S4Vectors::split(RangesObject,GenomeInfoDb::seqnames(RangesObject))
            private$Ranges.Col.Keys[[RangeKey]] <- list()
            private$Ranges.Keys <- c(private$Ranges.Keys,RangeKey)
        },
        LoadRangesByVectors = function(RangeKey=NULL,Chrom=NULL,Start=NULL,End=NULL,Strand=NULL,Names=NULL){
            if(any(is.null(c(length(Chrom),length(Start),length(End))))){
                stop("Chrom, Start, End cannot be NULL")
            }
            if(length(unique(length(Chrom),length(Start),length(End)))>1){
                stop("Chrom, Start, End seem to be of unequal length")
            }
            if((RangeKey %in% private$Protected.Ranges.Keys)) {
                stop("RangeKey cannot have a protected key value.",
                    paste(private$Protected.Ranges.Keys,collapse=" "))
            }
            Temp.Ranges <- private$MakeGRangesObject(Chrom=Chrom,Start=Start,End=End,Strand=Strand,Names=Names)
            private$RangesObjects[[RangeKey]] <- split(Temp.Ranges,seqnames(Temp.Ranges))
            private$Ranges.Col.Keys[[RangeKey]] <- list()
            private$Ranges.Keys <- c(private$Ranges.Keys,RangeKey)
        },
        UnloadRangesObject = function(RangeKey=NULL){
            require(GenomicRanges)
            if(is.null(RangeKey)) {
                stop(" RangeKey cannot be null..")
            }
            if(length(RangesObject)==0 | length(RangeKey) > 1) {
                stop("RangeKey length cannot be 0 or greater than 1")
            }
            if(RangeKey %in% private$Protected.Ranges.Keys) {
                stop("RangeKey cannot have a protected key value.",
                    paste(private$Protected.Ranges.Keys,collapse=" "))
            }
            if(!(RangeKey %in% private$Ranges.Keys)) {
                stop("RangeKey not found")
            }
            private$RangesObjects[[RangeKey]]=NULL
            private$Ranges.Keys[private$Ranges.Keys==RangeKey]=NULL
            private$Ranges.Col.Keys[[RangeKey]]=NULL
        },
        # Figure out what is the use of this particular function.
        FetchOverlapBetweenRanges = function(QueryKey=NULL,SubjectKey=NULL,Chrom=NULL,type="any"){
            AllTypes<-c("any","within")
            if( any(!(type %in% AllTypes)) ){
                stop("type takes one of two arguments: c(\"any\",\"within\")")
            }
            if(any(is.null(c(QueryKey,SubjectKey,Chrom)))){
                stop("QueryKey and SubjectKey cannot be empty.")
            }
            if(any(!(c(QueryKey,SubjectKey) %in% private$Ranges.Keys)) | 
                any(c(QueryKey,SubjectKey) %in% private$Protected.Ranges.Keys)) {
                stop("One of either QueryKey or SubjectKey not found.",
                    " Or one of the protected keys were provided.",
                    paste(private$Protected.Ranges.Keys,collapse=" "))
            }
            ChromosomeNames<-self$ChromosomeList
            if( !(Chrom %in% ChromosomeNames) ){
                stop(Chrom," does not exist in object..")
            }
            if(length(Chrom) > 1 | length(QueryKey) > 1 | length(SubjectKey) > 1){
                stop("QueryKey SubjectKey and Chrom cannot have length greater than 1")
            }
            if( any( !(Chrom %in% names(private$RangesObjects[[QueryKey]])) | 
                !(Chrom %in% names(private$RangesObjects[[SubjectKey]])) ) ){
                stop(Chrom,"doesnt exist for one of Query or Subject")
            }
            
            Query<-private$RangesObjects[[QueryKey]][[Chrom]]
            Subject<-private$RangesObjects[[SubjectKey]][[Chrom]]
                HitsObject<-findOverlaps(Query,Subject,type=type)
                UniqueQueries<-unique(queryHits(HitsObject))
                MatchingIndexes<-lapply(1:length(UniqueQueries),function(x){
                    A.Query<-UniqueQueries[x]
                    MatchingRanges<-subjectHits(HitsObject)[queryHits(HitsObject)==A.Query]
                    ListObj<-list()
                    ListObj[["SubjectInfo"]]<-MatchingRanges
                    ListObj[["QueryInfo"]]<-Query[A.Query]
                    ListObj
                    })
            names(OverlapByChromosome)<-Chrom
            return(OverlapByChromosome)
        },
        GetAttributes = function(){
            Attribute<-private$Attribute.List
            return(Attribute)
        },
        GetFilelist = function(){
            FileList<-private$File.List
            return(FileList)
        },
        GetRangeskeys = function(){
            Rangekeys<-private$Ranges.Keys
            return(Rangekeys)
        },
        GetDoneList = function(){
            DoneList<-private$Matrice.done
            return(DoneList)
        },
        GetRangesByChromosome = function(RangeKey=NULL,Chrom=NULL){
            require(GenomicRanges)
            if(any(is.null(c(RangeKey,Chrom)))){
                stop("RangeKey and Chrom are mandatory variables")
            }
            ChromosomeNames<-self$ChromosomeList
            RangeKeys<-private$Ranges.Keys
            if( !(RangeKey %in% RangeKeys) ){
                stop(RangeKey," does not exist in object..")
            }
            if( !(Chrom %in% ChromosomeNames) ){
                stop(Chrom," does not exist in object..")
            }
            ReturnObj<-private$RangesObjects[[RangeKey]][[Chrom]]
            return(ReturnObj)
        },
        Dimensions = function(Chrom1=NULL,Chrom2=NULL){
            if(is.null(Chrom1) | is.null(Chrom2)){
                stop("Chrom2 and Chrom2 key cannot be empty!")
            }
            ChromosomeList<-self$ChromosomeList
            if(length(Chrom1)>1 | length(Chrom2)>1){
                stop("Chrom1 length, Chrom2 length exceeded allowed length of 1.",
                    " This module is not iterative!")
            }
            if(!all(c(Chrom1,Chrom2) %in% ChromosomeList)){
                stop("Provided chromosomes were not found in chromosome list.")
            }
            if(!private$Matrice.done[Chrom1,Chrom2]){
                stop(Chrom1,Chrom2," matrix is yet to be loaded into the class.")
            }
            Chrom1.ranges <- self$GetRangesByChromosome(RangeKey="Bintable",Chrom=Chrom1)
            Chrom2.ranges <- self$GetRangesByChromosome(RangeKey="Bintable",Chrom=Chrom2)
            return(c(length(Chrom1.ranges),length(Chrom2.ranges)))
        },
        FetchRangeIndex = function(Chrom=NULL,Start=NULL,End=NULL,Names=NULL,type="any"){
            AllTypes<-c("any","within")
            if( any(!(type %in% AllTypes)) ){
                stop("type takes one of two arguments: c(\"any\",\"within\")")
            }
            if(is.null(Chrom) | is.null(Start) | is.null(End)){
                stop("Chrom, Start, End cannot be empty")
            }
            ChromosomeList<-self$ChromosomeList
            if( any(!(Chrom %in% ChromosomeList)) ){
                stop("Provided Chromosome does not exist in chromosome list.")
            }
            if(any(class(Chrom)!="character" | class(Start)!="numeric" | class(End)!="numeric")){
                stop("Provided chromosome,start,end do not match expected class definitions of character, numeric, numeric")
            }
            UniqueChromosomes<-unique(Chrom)
            OverlapByChromosome<-lapply(UniqueChromosomes,function(CurChrom){
                Cur.Chrom<-Chrom[Chrom==CurChrom]
                Cur.Start<-Start[Chrom==CurChrom]
                Cur.End<-End[Chrom==CurChrom]
                Cur.Names<-Names[Chrom==CurChrom]

                SubjectRanges<-self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=CurChrom)
                if( any(!(Cur.End <= max(end(SubjectRanges)) & Cur.Start >= min(start(SubjectRanges)))) ){
                    stop("Start or End is out of ranges for Bintable")
                }
                QueryRanges<-private$MakeGRangesObject(Chrom=Cur.Chrom,Start=Cur.Start,End=Cur.End,Names=Cur.Names)
                require(GenomicRanges)
                HitsObject<-GenomicRanges::findOverlaps(SubjectRanges,QueryRanges,type=type)
                UniqueQueries<-unique(S4Vectors::subjectHits(HitsObject))
                MatchingIndexes<-lapply(1:length(UniqueQueries),function(x){
                    A.Query<-UniqueQueries[x]
                    MatchingQueries<-S4Vectors::queryHits(HitsObject)[S4Vectors::subjectHits(HitsObject)==A.Query]
                    ListObj<-list()
                    ListObj[["Indexes"]]<-MatchingQueries
                    ListObj[["SubjectInfo"]]<-QueryRanges[A.Query]
                    ListObj
                    })
            })
            names(OverlapByChromosome)<-UniqueChromosomes
            return(OverlapByChromosome)
        },
        FetchMatrixMinMax = function(Chrom1=NULL,Chrom2=NULL){
            if(is.null(Chrom1) | is.null(Chrom2)) {
                stop("Either of Chrom1, Chrom2 were provided as NULL values")
            }
            private$CheckForChromosomes(Chrom1=Chrom1,Chrom2=Chrom2)
            if(!private$Matrice.done[Chrom1,Chrom2]){
                stop(Chrom1,Chrom2," matrix is yet to be loaded into the class.")
            }
            Ret.vector <- (private$Matrix.range[[Chrom1]][[Chrom2]])
            names(Ret.vector) <- c('min','max')
            return(Ret.vector)
        },
        FetchMatrixWithinCoords = function(XCoords=NULL,YCoords=NULL,type="within",FUN=NULL){
            if( (is.null(XCoords)) | (is.null(YCoords)) ){
                stop("XCoords and YCoords cannot be NULL")
            }
            if( !(length(XCoords)==1) | !(length(YCoords)==1) ){
                stop("This function processes single process calls at a time. ",
                    "Setup an Iterator for more functionality")
            }
            if( class(XCoords)!="character" | class(YCoords)!="character" ){
                stop("Two string variables were expected for XCoords & YCoords,\nfound XCoords class ",
                    class(XCoords),
                    " and YCoords class ",
                    class(YCoords))
            }            
            XCoord.split <- private$SplitGenomicCoordinates(Coordinate=XCoords)
            XCoord.Chrom <- XCoord.split[[1]][1]
            XCoord.start <- as.numeric(XCoord.split[[1]][2])
            XCoord.stop <- as.numeric(XCoord.split[[1]][3])
            YCoord.split <- private$SplitGenomicCoordinates(Coordinate=YCoords)
            YCoord.Chrom <- YCoord.split[[1]][1]
            YCoord.start <- as.numeric(YCoord.split[[1]][2])
            YCoord.stop <- as.numeric(YCoord.split[[1]][3])
            ChromosomeList<-self$ChromosomeList
            if( any(!(c(XCoord.Chrom,YCoord.Chrom) %in% ChromosomeList)) ){
                stop("Provided chromosomes were not found in chromosome list.")
            }
            Chrom1 = XCoord.Chrom
            Chrom2 = YCoord.Chrom
            if(!private$Matrice.done[Chrom1,Chrom2]){
                stop(Chrom1,Chrom2," matrix is yet to be loaded into the class.")
            }
            XCoords.Ranges<-self$FetchRangeIndex(Chrom=Chrom1,Start=XCoord.start,End=XCoord.stop,Names=NULL,type=type)
            YCoords.Ranges<-self$FetchRangeIndex(Chrom=Chrom2,Start=YCoord.start,End=YCoord.stop,Names=NULL,type=type)
            if( is.null(XCoords.Ranges[[Chrom1]][[1]][["Indexes"]]) | 
                is.null(YCoords.Ranges[[Chrom2]][[1]][["Indexes"]]) ){
                stop("Overlap operation was unsuccessful! Please check coordinates ",XCoords," & ",YCoords)
            }
            XVector<-XCoords.Ranges[[Chrom1]][[1]][["Indexes"]]
            YVector<-YCoords.Ranges[[Chrom2]][[1]][["Indexes"]]
            Matrix<-self$ExtractMatrix(Chrom1=Chrom1,Chrom2=Chrom2,XVector=as.numeric(XVector),YVector=as.numeric(YVector),FUN=FUN)
            return(Matrix)
        },
        FetchValuesByDistance = function(Chrom=NULL,Distance=NULL,Constrain.region=NULL,batch.size=500,
            sparsity.threshold=0.8,FUN=NULL){
            if(is.null(Distance) | is.null(Chrom)){
                stop("Distance and Chrom key cannot be empty!")
            }
            ChromosomeList<-self$ChromosomeList
            if(length(Chrom)>1 | length(Distance)>1 | length(Constrain.region)>1){
                stop("Chrom length, distance length. Constrain.region length exceeded allowed length of 1.",
                    " This module is not iterative!")
            }
            if(!(Chrom %in% ChromosomeList)){
                stop("Provided chromosome was not found in chromosome list.")
            }
            if(!private$Matrice.done[Chrom,Chrom]){
                stop(Chrom,Chrom," matrix is yet to be loaded into the class.")
            }           
            if(class(Distance)!="numeric"){
                stop("Provided Distance does not appear to be numeric")
            }
            Bin.Ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom)
            if(Distance < 0 | Distance >= length(Bin.Ranges)){
                stop("Provided Distance is out of bounds for Chromosome",Chrom,
                    "\nAs this is a symmetric matrices, please target the upper triangle.\n",
                    "P.S. Distance 0 corresponds to 1 as this module works on indices.")
            }
            Vector.start <- 1
            Vector.stop <- length(Bin.Ranges)
            if(!is.null(Constrain.region)){
                Vector.coordinates <- private$ReturnRegionPosition(Constrain.region=Constrain.region,Chrom=Chrom)

                if(is.null(Vector.coordinates)){
                    stop("Overlap operation was unsuccessful! Please check coordinates ",Constrain.region)
                }
                Vector.start <- min(Vector.coordinates)
                Vector.stop <- max(Vector.coordinates)
            }
            Starting.col <- Vector.start + Distance
            Count <- ((Vector.stop - Vector.start) + 1) - Distance 
            tot.mat.extracted <- Count * Count
            # cat("Boo ",Count,"\n")
            Start <- c(Vector.start,Starting.col)
            Stride <- c(1,1)
            Counts <- Count
            CumSums <- 0
            Groups <- c(private$hdf.matrices,Chrom)
            if(Count > batch.size){
                repeated <- floor(Count/batch.size)
                Counts <- rep(batch.size,repeated)
                if(repeated != ceiling(Count/batch.size)){
                    cumulative <- sum(Counts)
                    Counts <- c(Counts,(Count-cumulative))
                }
                CumSums <- cumsum(c(0,Counts[1:(length(Counts)-1)]))
            }
            DistancesVector.list <- lapply(1:length(Counts),function(x){
                Count <- Counts[x]
                Offset <- CumSums[x]
                cur.start <- Start+Offset
                # cat(Count,Offset,cur.start,"\n")
                diag(private$FetchFromDataset(Groups=Groups,Chrom=Chrom,Start=cur.start,Stride=c(1,1),Count=c(Count,Count)))
            })

            DistancesVector <- do.call(c,DistancesVector.list)

            # Indices <- list(c(Distance+1))
            # Groups <- c(private$hdf.distances.offset)
            # Starting.positions <- private$FetchFromDataset(Groups=Groups,Chrom=Chrom,Index=Indices)
            # Start <- Vector.start + Starting.positions
            # End <- (Vector.stop - Distance) + Starting.positions

            # Indices <- list(c(Start:End))

            # Groups=c(private$hdf.distances)
            # DistancesVector <- private$FetchFromDataset(Groups=Groups,Chrom=Chrom,Index=Indices)
            # Indices <- list(Vector.start:(Vector.stop-Distance))
            # Groups <- c(private$hdf.offset,Chrom)
            # Indices.offset <- private$FetchFromDataset(Groups=Groups,Chrom=Chrom,Index=Indices)
            # Indices <- list((Indices.offset + c(Vector.start:(Vector.stop-Distance))) + Distance)
            # Groups=c(private$hdf.matrices,Chrom)
            # DistancesVector <- private$FetchFromDataset(Groups=Groups,Chrom=Chrom,Index=Indices)

            if(is.null(FUN)){
                return(DistancesVector)
            }else{
                return(FUN(DistancesVector))
            }
        },
        ExtractMatrix = function(Chrom1=NULL,Chrom2=NULL,XVector=NULL,YVector=NULL,FUN=NULL){
            # cat(" Rows: ",XVector," Cols: ",YVector,"\n")
            if(any(!(class(XVector) %in% c("numeric","integer")) | !(class(YVector) %in% c("numeric","integer")))){
                stop("XVector and YVector must be numeric")
            }
            if(is.null(Chrom1) | is.null(Chrom2) | is.null(XVector) | is.null(YVector)){
                stop("Either of Chrom1, Chrom2, XVector or YVector were provided as NULL values")
            }
            private$CheckForChromosomes(Chrom1=Chrom1,Chrom2=Chrom2)
            if(!private$Matrice.done[Chrom1,Chrom2]){
                stop(Chrom1,Chrom2," matrix is yet to be loaded into the class.")
            }
            Group <- c(private$hdf.matrices,Chrom1)
            Chrom1.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom1)
            Chrom2.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom2)         
            if(any(XVector > length(Chrom1.ranges)) | any(YVector > length(Chrom2.ranges)) | min(XVector,YVector) < 1 ) {
                stop("XVector or YVector falls outside the bounds of loaded Bintables") 
            }
            Matrix<-self$FetchVectorValues(Chrom1=Chrom1,Chrom2=Chrom2,xaxis=XVector,yaxis=YVector)
            if(is.null(FUN)){
                return(Matrix)              
            }else{
                return(FUN(Matrix)) 
            }
        },
        FetchRowOrColVector = function(Chrom1=NULL,Chrom2=NULL,by=NULL,vector=NULL,Constrain.regions=NULL,FUN=NULL){
            Chrom.all <- c(Chrom1,Chrom2)
            if(is.null(Chrom1) | is.null(Chrom2) | is.null(by) | is.null(vector)) {
                stop("Either of Chrom1, Chrom2, by or vector were provided as NULL values")
            }
            if(!(by %in% c("position","ranges")) | length(by) != 1){
                stop("by expects a vector of type character, length 1 and takes either one of position or ranges as values")
            }
            private$CheckForChromosomes(Chrom1=Chrom1,Chrom2=Chrom2)
            if(!private$Matrice.done[Chrom1,Chrom2]){
                stop(Chrom1,Chrom2," matrix is yet to be loaded into the class.")
            }
            if(by=="position"){
                Class.type <- c("numeric","integer")
            }
            if(by=="ranges"){
                Class.type <- "character"
            }
            if(!(class(vector) %in% Class.type)){
                stop("vector must be of class",
                    ifelse(length(Class.type)>1,paste(Class.type,collapse=" or "),paste(Class.type)),"when by has value",by)
            }
            if(length(Chrom.all)>2){
                stop("This module is not iterable")
            }
            Chrom1.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom1)
            Chrom2.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom2)
            if(by=="ranges"){
                if(any(!(vector %in% names(Chrom1.ranges)))){
                    stop("All ranges not found in ",Chrom1," Bintable")
                }
                Positions<-which(names(Chrom1.ranges) %in% vector)
                names(Positions) <- vector
            }else{
                if(any(max(vector)>max(length(Chrom1.ranges)) | min(vector)<1)) {
                    stop("Position vector falls outside the bounds of ",Chrom1," Bintable") 
                }               
                Positions <- vector
                names(Positions) <- names(Chrom1.ranges[vector])
            }
            PikaPika<-FALSE
            if(!is.null(Constrain.regions) & length(Constrain.regions)==length(vector)){
                PikaPika<-TRUE
            }
            Seq.indices <- seq_along(Positions)
            names(Seq.indices) <- names(Positions)
            Vector.values<-lapply(Seq.indices,function(ind){
                x <- Positions[ind]
                if(PikaPika){
                    Constrain.region <- Constrain.regions[ind]
                    y <- private$ReturnRegionPosition(Constrain.region=Constrain.region,Chrom=Chrom2)
                }else{
                    y <- 1:length(Chrom2.ranges)
                }
                Chrom2.names <- names(Chrom2.ranges[y])

                Values <- self$FetchVectorValues(Chrom1=Chrom1,Chrom2=Chrom2,xaxis=x,yaxis=y,FUN=FUN)
                return(Values)
            }) 
            return(Vector.values)
        },
        FetchVectorValues = function(Chrom1=NULL,Chrom2=NULL,xaxis=NULL,yaxis=NULL,FUN=NULL){
            if(is.null(Chrom1) | is.null(Chrom2)){
                stop("Chrom1 and Chrom2 keys cannot be empty!")
            }
            ChromosomeList<-self$ChromosomeList
            if(class(Chrom1)!="character" | class(Chrom2)!="character"){
                stop("Provided Chromosomes does not appear to be of class character")
            }
            if(is.null(xaxis) & is.null(yaxis)){
                stop("Both xaxis and yaxis cannot be null")
            }
            Start <- c(min(xaxis),min(yaxis))
            Stride <- c(1,1)
            Count <- c(length(xaxis),length(yaxis))
            Vector <- private$FetchFromDataset(Groups=c(private$hdf.matrices,Chrom1),Chrom=Chrom2,Start=Start,
                Stride=Stride,Count=Count)
            if(is.null(FUN)){
                return(Vector)
            }else{
                return(FUN(Vector))
            }
        },
        FetchMatrixHandler = function(Chrom1=NULL,Chrom2=NULL){
            if(is.null(Chrom1) | is.null(Chrom2)){
                stop("Chrom1 and Chrom2 keys cannot be empty!")
            }
            ChromosomeList<-self$ChromosomeList
            if(max(unique(length(Chrom1),length(Chrom2)))>1){
                stop("Chrom1 or Chrom2 length exceeded allowed length of 1. This module is not iterative!")
            }
            private$CheckForChromosomes(Chrom1=Chrom1,Chrom2=Chrom2)
            if(class(Chrom1)!="character" | class(Chrom2)!="character"){
                stop("Provided Chrom1 or Chrom2 does not appear to be of class character")
            }

            Chrom2.Dataset <-private$ReturnH5DatasetHandler(Group=private$hdf.matrices,
                Chrom1=Chrom1,Chrom2=Chrom2)

            return(Chrom2.Dataset)
        },
        ComputeRowOrColCoverage = function(){
            warning("DEPRECATION NOTICE: Parts of this function are now done on the fly.\nBasically it's useless. For now.")
            ChromosomeList<-self$ChromosomeList
            for (i in 1:length(ChromosomeList)) {
                Chrom<-ChromosomeList[i]
                if(!private$Matrice.done[Chrom,Chrom]){
                    stop("File is yet to be loaded for Chromosome",Chrom)
                }
            }
            PercentGTZero = function(x){
                LengthOfRow<-length(x)
                LengthOfGTZero<-length(x[x>0])
                Fraction<-LengthOfGTZero/LengthOfRow
                return(Fraction)
            }
            for (i in 1:length(ChromosomeList)) {
                Chrom<-ChromosomeList[i]
                CurrentRanges<-self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom)
                RangesNames<-names(CurrentRanges)
                VectorList<-self$FetchRowOrColVector(Chrom1=Chrom,Chrom2=Chrom,by="ranges",vector=RangesNames,FUN=PercentGTZero)
                VectorList[is.na(VectorList) | is.infinite(VectorList)] <- 0
                private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chrom,Column=FractionNonZero,
                    Col.Name="Bin.coverage",Replace=FALSE,na.function=as.numeric)
                VectorSums <- self$FetchRowOrColVector(Chrom1=Chrom,Chrom2=Chrom,by="ranges",vector=RangesNames,FUN=sum)
                VectorSums[is.na(VectorList) | is.infinite(VectorList)] <- 0
                private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chrom,Column=VectorSums,
                    Col.Name="Row.sums",Replace=FALSE,na.function=as.numeric)
            }
        },
        RecomputeSparsityIndex = function(sparsity.compute.bins=100){
            if(self$is.sparse()){
                if(length(sparsity.compute.bins)>1 | !is.numeric(sparsity.compute.bins)){
                    stop("sparsity.compute.bins expects numeric value of length 1\n")
                }
                private$Sparsity.compute.bins=sparsity.compute.bins
            }
            if(self$is.sparse()){
                private$ComputeSparsityIndex()
            }
        },
        is.sparse = function(){
            return(private$ComputeSparsity)
        },
        sparsity.index.parameter = function(){
            return(private$Sparsity.compute.bins)  
        }
    ),
    private = list(
        Attribute.List=NA,
        File.List=NA,
        RangesObjects=NA,
        Ranges.Keys=NA,
        Ranges.Col.Keys=NA,
        WhichStarter=NA,
        Working.Dir=NA, 
        Working.File=NA,
        HDF.Connection=NA,
        ComputeSparsity=FALSE,
        Sparsity.compute.bins=NA,
        hdf.matrices = "matrices",
        hdf.ranges = "base.ranges.tables",
        Max.vector.size=104857600,
        Matrix.range=NA,
        Num.lines=1,
        hdf.root.folders=c("matrices","base.ranges.tables"),
        Protected.Ranges.Keys=c("Bintable"),
        Bintable.Key = "Bintable",
        NonStrandedColNames=c("chr","start","end"),
        Matrice.done = NA,
        Ranges.separator=":",
        CreateFileConnection = function(Filename=NULL,mode="r"){
            if(grepl('$/.gz',Filename)){
                connection=gzfile(Filename,open=mode)
            }else if(grepl('$/.bz2',Filename)){
                connection=bzfile(Filename,open=mode)
                }else {
                    connection=file(Filename,mode)
                }
                return(connection)
        },
        ReturnH5FileConnection = function(){
            require(rhdf5)
                HDF.File <- file.path(private$Working.Dir,private$Working.File)
                private$HDF.Connection = H5Fopen(name=HDF.File)
                return(private$HDF.Connection)
        },
        CreateTempH5FileConnection = function(File=NULL){
            require(rhdf5)
            HDF.File <- File
            private$HDF.Connection = H5Fopen(name=HDF.File)
            return(private$HDF.Connection)
        },
        CloseH5FileConnection = function(){
            H5Fclose(private$HDF.Connection)
        },
        ReturnH5GroupHandler = function(Groups=NULL){
                Temp.connection = private$ReturnH5FileConnection()
                Path.to.file <- ""
                for (Group in Groups) {
                    Path.to.file <- file.path(Path.to.file,Group)
                }
                Group.connection=Temp.connection&Path.to.file
                return(Group.connection)
        },
        ReturnH5DatasetHandler = function(Group=NULL,Chrom1=NULL,Chrom2=NULL){
                Temp.connection = private$ReturnH5FileConnection()
                Path.to.file <- file.path("",Group,Chrom1,Chrom2)
                Chrom2.data.connection=Temp.connection&Path.to.file
                return(Chrom2.data.connection)
        },
        InsertIntoDataset = function(Groups=NULL,Connection=NULL,Chrom=NULL,Data=NULL,Index=NULL,
            Start=NULL,Stride=NULL,Count=NULL){
            if(!is.null(Connection) & is.null(Groups)){
                Temp.connection <- Connection
            }else if(is.null(Connection) & !is.null(Groups)){
                Temp.connection = private$ReturnH5GroupHandler(Groups=Groups)
            }else {
                stop("InsertIntoDataset expects one of Groups or Connection to be NULL and !NULL")
            }
            if(!is.null(Index)){
                h5writeDataset(obj=Data, h5loc=Temp.connection, name=Chrom, index=Index)
            }else if(!is.null(Count) & !is.null(Stride) & !is.null(Start)){
                h5writeDataset(obj=Data, h5loc=Temp.connection, name=Chrom, start=Start, stride=Stride, count=Count)
            }else{
                stop("FetchFromDataset expects one of Start, Stride, Count if Index is NULL")
            }
            private$Flush.HDF()
        },
        FetchFromDataset = function(Groups=NULL,Connection=NULL,Chrom=NULL,Index=NULL,
            Start=NULL,Stride=NULL,Count=NULL,Block=NULL){

            if(!is.null(Connection) & is.null(Groups)){
                Temp.connection <- Connection
            }else if(is.null(Connection) & !is.null(Groups)){
                Temp.connection = private$ReturnH5GroupHandler(Groups=Groups)
            }else {
                stop("FetchFromDataset expects one of Groups or Connection to be NULL and !NULL")
            }
            if(!is.null(Index)){
                Data<- h5read(file=Temp.connection,name=Chrom,index=Index)              
            }else if(is.null(Count) | is.null(Stride) | is.null(Start)){
                stop("FetchFromDataset expects one of Start, Stride, Count if Index is NULL")   
            }
            Data <- h5read(file=Temp.connection,name=Chrom,start=Start,stride=Stride,count=Count,block=Block)
            if(is.null(Connection) & !is.null(Groups)){
                H5Gclose(Temp.connection)
                private$CloseH5FileConnection()
            }
            return(Data)
        },
        ReturnRegionPosition = function(Constrain.region=NULL,Chrom=NULL){
            if(!is.character(Constrain.region) | length(Constrain.region) > 1){
                stop("Constrain.region must be a character vector of length 1")
            }
            Coord.Split<-private$SplitGenomicCoordinates(Coordinate=Constrain.region)
            Region.Chrom<-Coord.Split[[1]][1]
            Region.start<-as.numeric(Coord.Split[[1]][2])
            Region.stop<-as.numeric(Coord.Split[[1]][3])
            if(Chrom!=Region.Chrom){
                stop("Constrain.region Chrom and provided Chrom must be same.") 
            }
            Region.Ranges<-self$FetchRangeIndex(Chrom=Chrom,Start=Region.start,End=Region.stop,
                Names=NULL,type="within")
            Vector.coordinates <- Region.Ranges[[Chrom]][[1]][["Indexes"]]
            return(Vector.coordinates)
        },
        Flush.HDF = function(){
            require(rhdf5)
            H5Fflush(private$HDF.Connection)            
        },
        AddFileToList = function(Key1=NULL,Key2=NULL,Value=NULL){
            private$File.List[Key1,Key2]<-Value
        },
        AddDoneForChromosome = function(Key1=NULL,Key2=NULL){
            private$Matrice.done[Key1,Key2] <- as.logical(1)
        },
        AddAttribute = function(Key=NULL,Value=NULL){
            private$Attribute.List[[Key]]<-Value
        },
        MakeGRangesObject = function(Chrom=NULL, Start=NULL, End=NULL, Strand=NULL, Names=NULL){
            require(GenomicRanges)
            if(is.null(Names)){
                Names<-paste(Chrom,as.integer(Start),as.integer(End),sep=private$Ranges.separator)
            }
            if(is.null(Strand)){
                Strand<-rep("*",length(Chrom))
            }
                Object<-GenomicRanges::GRanges(
                    seqnames=Rle(Chrom),
                    ranges=IRanges(Start,end=End,names=Names),
                    strand=Rle(strand( Strand )))
                return(Object)
        },
        CheckNonStrandedTableForErrors = function(Table=NULL, ObsColNum=NULL, ExpColNum=NULL){
            if(ObsColNum!=ExpColNum){
                stop("Expected Data format is chr, start, stop. Default delimiters are //s and //n")
            }else{
                colnames(Table)<-private$NonStrandedColNames
            }
            if(class(Table[,1])=="character"){
                UniqueChromNames<-unique(Table[,1])
                if(any( !(UniqueChromNames %in% self$ChromosomeList) )){
                    stop("Some chromosome names are not defined in the chromosome table")
                }
            }
            if(class(Table[,2])!="numeric" | class(Table[,3])!="numeric"){
                stop("Numeric values expected at column 2 and 3")
            }
            if(any(!(Table[,2] %% 1 == 0)) | any(!(Table[,3] %% 1 == 0))) {
                stop("Genomic coordinates in column 2 and 3 cannot have float values")
            }
            if( any( (Table[,2] > Table[,3]) ) ){
                stop("start coordinates cannot be greater than stop coordinates")
            }
            return(Table)
        },
        CheckContinuousRanges = function(Table=NULL, StartCol=NULL, EndCol=NULL){
            Starts<-Table[,StartCol]
            Starts<-Starts[2:length(Starts)]
            End<-Table[,EndCol]
            End<-End[1:(length(End)-1)]
            if( any(Starts==End) ){
                stop("Found continuous ranges in file! Cannot proceed further!")
            }
        },
        ReadNonStrandCoordsToGRanges = function(Filename=NULL,read.delim=" ",exec="cat",RangeKey=NULL){
            require(data.table)
            if(is.null(RangeKey)){
                stop("Variable name not provided while reading table.")
            }
            if(is.null(exec)){
                stop("exec is not allowed to be null")
            }
            Colnames<-private$NonStrandedColNames
            Command <- paste(exec,Filename,sep=" ")
            Table <- fread(input=Command, sep=read.delim, colClasses=c("character","numeric","numeric"),
                stringsAsFactors=FALSE, verbose=FALSE, showProgress=FALSE, col.names=Colnames, data.table=FALSE)
            NumberOfColumns<-ncol(Table)
            
            ValidatedTable<-private$CheckNonStrandedTableForErrors(Table=Table,ObsColNum=NumberOfColumns,ExpColNum=3)
            
            private$CheckContinuousRanges(Table=ValidatedTable,StartCol=c("start"),EndCol=c("end"))
            ValidatedTable<-ValidatedTable[order(ValidatedTable$chr,ValidatedTable$start),]
            
            Chrom<-ValidatedTable$chr
            Start<-ValidatedTable$start
            End<-ValidatedTable$end

            RangesObject<-private$MakeGRangesObject(Chrom=Chrom,Start=Start,End=End) 

            private$RangesObjects[[RangeKey]] = split(RangesObject,seqnames(RangesObject))
            private$Ranges.Col.Keys[[RangeKey]]<-list()
            private$Ranges.Keys<-c(private$Ranges.Keys,RangeKey)
        },
        CreateFileList = function(ChromNames){
            private$File.List = matrix(data="",nrow=length(ChromNames),ncol=length(ChromNames))
            rownames(private$File.List)<-ChromNames
            colnames(private$File.List)<-ChromNames
        },
        CreateDoneList = function(ChromNames){
            private$Matrice.done = matrix(nrow=length(ChromNames),ncol=length(ChromNames),data=as.logical(0))
            rownames(private$Matrice.done)<-ChromNames
            colnames(private$Matrice.done)<-ChromNames
        },
        InitiateAttribute = function(){
            private$Attribute.List=list()
        },
        AddData = function(Filename=NULL,delim=" ",exec=NULL,Chromosome1=NULL,Chromosome2=NULL,fix.num.rows.at=2000){
            require(data.table)
            PercentGTZero = function(x){
                x[is.na(x) | is.infinite(x)] <- 0
                LengthOfRow <- length(x)
                LengthOfGTZero <- length(x[x!=0])
                Fraction <- LengthOfGTZero/LengthOfRow
                return(Fraction)
            }
            ComputeRowSums = function(x){
                x[is.na(x) | is.infinite(x)] <- 0
                return(sum(x))
            }
            ComputeMinMax = function(x){
                x[is.na(x) | is.infinite(x)] <- 0
                return(c(min(x),max(x)))
            }
            if(self$is.sparse()){
                Sparsity.bins = self$sparsity.index.parameter()
                SparsityIndex = function(x=NULL,index=NULL,length=NULL){
                    x[is.na(x) | is.infinite(x)] <- 0
                    Range <- (index-private$Sparsity.compute.bins):(index+private$Sparsity.compute.bins)
                    Range <- Range[Range>0 & Range<length]
                    Rows <- x[Range]
                    return(length(Rows[Rows!=0])/length(Rows))
                }
            }
            Chromosomes.all <- c(Chromosome1,Chromosome2)
            if(is.null(exec)){
                stop("exec takes as input the shell command to be used by data.table for reading")
            }
            if((!is.null(Chromosome1) & !is.null(Chromosome2))){
                private$CheckForChromosomes(Chrom1=Chromosome1,Chrom2=Chromosome2)
            }
            if(!is.null(Filename)){
                Command <- paste(exec,Filename,sep=" ")
                Start.row <- 0
                Path.to.file <- file.path(private$Output.Directory,private$Output.Filename)
                Chrom1.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chromosome1)
                Chrom2.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chromosome2)
                Chrom1.len <- length(Chrom1.ranges)
                Chrom2.len <- length(Chrom2.ranges)
                Chrom1.Group <- private$ReturnH5GroupHandler(Groups=c(private$hdf.matrices,Chromosome1))
                Cumulative.data <- NULL
                Cumulative.distances.data <- NULL
                Cumulative.indices <- NULL
                NumLines <- private$FindLineNumbers(Row.len=Chrom1.len,Col.len=Chrom2.len)
                if(NumLines <= fix.num.rows.at){
                    NumLines <- fix.num.rows.at
                }
                if(Chrom1.len <= NumLines){
                    NumLines <- Chrom1.len
                }
                Bin.coverage <- NULL
                Row.sums <- NULL
                Sparsity.Index <- NULL
                Iterations.number <- Chrom1.len / NumLines
                Iterations <- rep(NumLines,floor(Iterations.number))
                if(floor(Iterations.number)!=ceiling(Iterations.number)){
                    cumulative <- sum(Iterations)
                    Iterations <- c(Iterations,(Chrom1.len-cumulative))
                }
                Skippity<-0
                if(length(Iterations)>1){
                    Skippity.cumsum <- cumsum(Iterations)
                    Skippity <- c(0,Skippity.cumsum[1:(length(Skippity.cumsum)-1)])
                }
                i<-1
                while(i<=length(Iterations)) {
                    Iter <- Iterations[i]
                    Skip <- Skippity[i]
                    Matrix <- as.matrix(fread(input=Command, sep=delim, nrows=Iter, na.strings="NA", 
                        stringsAsFactors=FALSE, skip=Skip,verbose=FALSE, dec=".",
                        showProgress=FALSE))
                    cat("Read",Iter,"lines after Skipping",Skip,"lines\n")
                    Bin.coverage <- c(Bin.coverage,sapply(1:nrow(Matrix),function(x){
                        Vec.sub <- Matrix[x,]
                        PercentGTZero(Vec.sub)
                    }))
                    Row.sums <- c(Row.sums,sapply(1:nrow(Matrix),function(x){
                        Vec.sub <- Matrix[x,]
                        ComputeRowSums(Vec.sub)
                    }))
                    if(self$is.sparse() & Chromosome1==Chromosome2){
                        Sparsity.Index <- c(Sparsity.Index,sapply(1:nrow(Matrix),function(x){
                            Vec.sub <- Matrix[x,]
                            sparsity.bin.idexes <- x + Skip
                            SparsityIndex(x=Vec.sub,index=sparsity.bin.idexes,length=Chrom2.len)
                        }))
                    }
                    Row.extent <- ComputeMinMax(Matrix)
                    if(private$Matrix.range[[Chromosome1]][[Chromosome2]][1] > Row.extent[1]) {
                        private$Matrix.range[[Chromosome1]][[Chromosome2]][1] <- Row.extent[1]
                    }
                    if(private$Matrix.range[[Chromosome1]][[Chromosome2]][2] < Row.extent[2]){
                        private$Matrix.range[[Chromosome1]][[Chromosome2]][2] <- Row.extent[2]  
                    }
                    Cumulative.data <- rbind(Cumulative.data,Matrix)
                    Obj.size <- object.size(Cumulative.data)
                    if(Obj.size>=private$Max.vector.size | i==length(Iterations)){
                        Start <- c(Start.row+1,1) 
                        Stride <- c(1,1)
                        Count <- c(nrow(Cumulative.data),ncol(Cumulative.data))
                        cat("Inserting Data at location:",Start[1],"\n")
                        cat("Data length:",Count[1],"\n")
                        private$InsertIntoDataset(Connection=Chrom1.Group,
                            Chrom=Chromosome2,Data=Cumulative.data,Start=Start,Stride=Stride,Count=Count)
                        Start.row <- Start.row+Count[1]
                        Cumulative.data <- NULL
                        cat("Loaded ",Obj.size," bytes of data...\n")
                        
                    }
                    cat("Read ",(Skip+Iter),"records...\n")
                    i<-i+1
                }
                H5Gclose(Chrom1.Group)
                private$CloseH5FileConnection()
                # cat(length(Bin.coverage),"\n")
                private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Bin.coverage,
                    Col.Name=paste(Chromosome2,"cov",sep="."),Replace=TRUE,na.function=as.numeric)
                # cat(length(Row.sums),"\n")
                private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Row.sums,
                    Col.Name=paste(Chromosome2,"sum",sep="."),Replace=TRUE,na.function=as.numeric)
                if((Chromosome1 == Chromosome2) & self$is.sparse()){
                    private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chromosome1,Column=Sparsity.Index,
                        Col.Name=paste("sparsity","idx",sep="."),Replace=TRUE,na.function=as.numeric)                    
                }
                private$AddFileToList(Key1=Chromosome1,Key2=Chromosome2,Value=normalizePath(Filename))
                private$AddDoneForChromosome(Key1=Chromosome1,Key2=Chromosome2)
            }
        },
        AddMetaDataToRanges = function(RangeKey=NULL,Chrom=NULL,Column=NULL,
            Col.Name=NULL,Replace=FALSE,na.function=as.numeric){
            require(GenomicRanges)
            ChromosomeNames<-self$ChromosomeList
            RangeKeys<-private$Ranges.Keys
            if( !(RangeKey %in% RangeKeys) ){
                stop(RangeKey," does not exist in object..")
            }
            if( !(Chrom %in% ChromosomeNames) ){
                stop(Chrom," does not exist in object..")
            }
            if(any(is.null(Column) | is.null(Col.Name))){
                stop("Column and Col.Name cannot be empty") 
            }
            RangeKeys<-private$Ranges.Keys
            if( !(RangeKey %in% RangeKeys) ){
                stop(RangeKey," does not exist in object..")
            }
            if(!(Col.Name %in% private$Ranges.Col.Keys[[RangeKey]][[Chrom]]) | Replace){
                private$Ranges.Col.Keys[[RangeKey]][[Chrom]]<-unique(c(private$Ranges.Col.Keys[[RangeKey]][[Chrom]],
                    Col.Name))
            }else if((Col.Name %in% private$Ranges.Col.Keys[[RangeKey]][[Chrom]])){
                stop("Col.Name already exists for ranges",RangeKey)
            }
            if(!(Col.Name %in% colnames(elementMetadata(private$RangesObjects[[RangeKey]][[Chrom]])))){
                TempObject <- unlist(private$RangesObjects[[RangeKey]],use.names=FALSE)
                S4Vectors::elementMetadata(TempObject)[Col.Name] <- na.function(NA)
                private$RangesObjects[[RangeKey]] <- split(TempObject, seqnames(TempObject))
            }
            elementMetadata(private$RangesObjects[[RangeKey]][[Chrom]])[Col.Name] <- Column
        },
        CheckForChromosomes = function(Chrom1=NULL,Chrom2=NULL){
            Chrom.all <- c(Chrom1,Chrom2)
            if(any(!(Chrom.all  %in% self$ChromosomeList))){
                stop(Chrom.all[which(!(Chrom.all  %in% self$ChromosomeList))]," was not found in ChromosomeList\n")
            }
        },
        SplitGenomicCoordinates = function(Coordinate=NULL){
            require(stringr)
            Sep <- private$Ranges.separator
            Coord.Split<-stringr::str_split(pattern=Sep,string=Coordinate)
            if(length(Coord.Split[[1]])!=3 | length(Coord.Split[[1]])!=3){
                stop("Coordinate must be separated by :")
            }

            Chrom<-Coord.Split[[1]][1]
            start<-as.numeric(Coord.Split[[1]][2])
            stop<-as.numeric(Coord.Split[[1]][3])
            if(any(class(Chrom)!="character" | class(start)!="numeric" | class(stop)!="numeric")){
                stop("Provided chromosome,start,end do not match expected class definitions of character, numeric, numeric")
            }
            return(Coord.Split)
        },
        FindLineNumbers = function(Row.len=NULL,Col.len=NULL){
            pixel.mem <- 48
            row.pixel.mem <- pixel.mem * Col.len
            Batch.size <-  floor(private$Max.vector.size/row.pixel.mem) 
            Batch.size[Batch.size<1] <- 1
            return(Batch.size)
        },
        ComputeSparsityIndex = function(){
            private$ComputeSparsityIndex()
            Sparsity.bins = self$sparsity.index.parameter()
            SparsityIndex = function(x=NULL){
                x[is.na(x) | is.infinite(x)] <- 0
                return(length(x[x!=0])/length(x))
            }
            for(Chrom1 in self$ChromosomeList){
                Chrom1.table <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom1)           
                if(private$Matrice.done[Chrom1,Chrom1]){
                    Sparsity.Index <- sapply(1:length(Chrom1.table),function(x){
                        Range <- (x - Sparsity.bins) : (x + Sparsity.bins)
                        Range <- Range[Range>0 & Range<length(Chrom1.table)]
                        Constrain.region <- paste(Chrom1,min(start(Chrom1.table[Range])),max(end(Chrom1.table[Range])),sep=private$Ranges.separator)
                        Row.vector <- unlist(self$FetchRowOrColVector(Chrom1=Chrom1,Chrom2=Chrom1,
                            by="position",vector=x,Constrain.region=Constrain.region),use.names=FALSE)
                        SparsityIndex(Row.vector)
                    })
                    private$AddMetaDataToRanges(RangeKey=private$Bintable.Key,Chrom=Chrom1,Column=Sparsity.Index,
                        Col.Name=paste("sparsity","idx",sep="."),Replace=TRUE,na.function=as.numeric)
                }
            }
        }
        # ,
        # CheckCoolerBintableCompatibility <- function(File=NULL,Group=NULL,Chrom=Bin.chrom.dataset,
        #     Start=NULL,End=NULL,from.HDF=FALSE){
        #     # Bin.group
        #     # Bin.start.dataset
        #     # Bin.end.dataset
        #     File.Connection <- private$CreateTempH5FileConnection(File=File)
        #     if(from.HDF){
        #         Chrom <- FetchFromDataset(Connection=File.Connection,Chrom=Group,)
        #     }
        # }
    )
)