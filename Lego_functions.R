# HDF structure
# Base.matrices <group>
#   chromsome <group>
#       chromosome <group> Attributes FileName, Done
#           matrix <dataset>
#           extent <dataset>
#           bla bla datasets
# Base.ranges <group>
#   Names <dataset> containing the names of all the ranges tables listed under Base.ranges
#   Bintable <group> containing binary attributes Strand, Names
#       Names <dataset>
#       ranges <dataset> based on value of attributes 3,4,5 column
#       BlaBla1 <dataset>
#       BlaBla2 <dataset>
# Base.metadata <group>
# chromosomes <dataset>

CreateLego <- function(ChromNames=NULL, BinTable=NULL, bin.delim="\t",
        col.index=c(1,2,3), strand=NULL, names=NULL, impose.discontinuity=TRUE,
        ChunkSize=NULL, Output.Filename=NULL, exec="cat", remove.existing=FALSE,
        sparse=FALSE, sparsity.compute.bins=100){
            require(rhdf5)
            H5close()

            Working.File <- normalizePath(Output.Filename)
			Reference.object <- GenomicMatrix$new()
            Root.folders <- Reference.object$GetRootFolders()
            if(is.null(ChromNames) | length(ChromNames) == 0){
                stop("Variable ChromNames cannot be empty")   
            }

            HDF.File <- file.path(Working.File)
            if(file.exists(HDF.File)){
                if(remove.existing){
                    file.remove(HDF.File)
                }else{
                    stop("Provided HDF file already exists. Please provide remove.existing = TRUE to overwrite it\n")
                }
            }

            ChromosomeList<-ChromNames

            if(is.null(BinTable)){
                stop("Variable Bintable cannot be empty. Binning information must be provided at startup")
            }
            
            cat("Reading Bintable:",BinTable,"\n")
            Bintable <- Read_bintable(Filename=BinTable,read.delim=bin.delim,exec=exec,
                        col.index=col.index,strand=strand,names=names,chromosomes=ChromNames,
                        impose.discontinuity=impose.discontinuity)

            h5createFile(HDF.File)
            for (Folder in Root.folders) {
                CreateGroups(Groups = Folder, File = HDF.File)
            }

            Chrom.lengths <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = length, col.name = 'chr')
            Chrom.sizes <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = max, col.name = 'end')
            Chrom.info.df <- data.frame(chr = names(Chrom.lengths),
                nrow = as.vector(Chrom.lengths),
                size = as.vector(Chrom.sizes))
            # Create metadata chromosome groups
            Lego_WriteDataFrame(Lego = HDF.File, Path = c(Root.folders['metadata'],
                Reference.object$metadata.chrom.dataset()), object = Chrom.info.df)

            # Create matrices groups
            for (chrom1 in ChromosomeList) {
                CreateGroups(Groups = c(Root.folders['matrices'],chrom1), File = HDF.File)
                for (chrom2 in ChromosomeList) {
                    CreateGroups(Groups = c(Root.folders['matrices'],chrom1,chrom2), File = HDF.File)
                    # h5createAttribute(obj = file.path(Root.folders['matrices'],chrom1,chrom2),file = HDF.File,
                    #     attr = "Filename", dims = 1)
                    # h5createAttribute(obj = file.path(Root.folders['matrices'],chrom1,chrom2),file = HDF.File,
                    #     attr = "Done", dims = 1)
                    Dims <-c(length(Chrom1.ranges),length(Chrom2.ranges))
                    if(is.null(ChunkSize)){
                        ChunkSize <- ceiling(Dims/100)  
                    }
                    h5createDataset(file=Chrom1.Group, dataset=Chrom2, dims=Dims, 
                        maxdims = Dims, H5type="H5T_NATIVE_DOUBLE", chunk = ChunkSize,
                        level = 6, fillValue=0, showWarnings = TRUE)
                }
            }
            Matrices.HDF <- H5Gopen(h5loc=HDF.Handler, name=Root.folders['matrices'])
            for (Chrom1 in ChromosomeList) {
                H5Gcreate(h5loc=Matrices.HDF, name=Chrom1)
                private$Matrix.range[[Chrom1]] <- list()
                Chrom1.Group <- H5Gopen(h5loc=Matrices.HDF, name=Chrom1)
                Chrom1.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom1)
                for (Chrom2 in self$ChromosomeList) {
                    private$Matrix.range[[Chrom1]][[Chrom2]] <- c(0,0)
                    Chrom2.ranges <- self$GetRangesByChromosome(RangeKey=private$Bintable.Key,Chrom=Chrom2)
                    # ChunkSize = floor(chunk.size)
                    Dims <-c(length(Chrom1.ranges),length(Chrom2.ranges))
                    if(is.null(ChunkSize)){
                        ChunkSize <- ceiling(Dims/100)  
                    }
                    h5createDataset(file=Chrom1.Group, dataset=Chrom2, dims=Dims, 
                        maxdims = Dims, H5type="H5T_NATIVE_DOUBLE", chunk = ChunkSize,
                        level = 6, fillValue=0, showWarnings = TRUE)
                }
                H5Gclose(Chrom1.Group)
            }
            H5Gclose(Matrices.HDF)
            H5Fclose(HDF.Handler)
}

Lego_WriteDataFrame <- function(Lego = NULL, Path = NULL, object = NULL){
    library(stringr)
    if(!(length(c(Lego,Path,object))>=4)){
        stop("All arguments are required!")
    }
    Dirs <- Path[1:(length(Path)-1)]
    Lego.handler <- ReturnH5GroupHandler(Groups = Dirs, File = Lego)
    h5writeDataset.data.frame(h5loc = Lego.handler, 
        obj = Chrom.len.df, 
        name = Path[length(Path)])
    H5Fclose(Lego.handler)
}





AddAttribute(Key=private$Bintable.Key,Value=normalizePath(BinTable))


if(sparse){
    if(length(sparsity.compute.bins)>1 | !is.numeric(sparsity.compute.bins)){
        stop("sparsity.compute.bins expects numeric value of length 1\n")
    }
    private$ComputeSparsity=sparse
    private$Sparsity.compute.bins=sparsity.compute.bins
}