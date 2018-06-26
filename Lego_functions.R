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
            # Read in the binning table
            cat("Reading Bintable:",BinTable,"\n")
            Bintable <- Read_bintable(Filename=BinTable,read.delim=bin.delim,exec=exec,
                        col.index=col.index,strand=strand,names=names,chromosomes=ChromNames,
                        impose.discontinuity=impose.discontinuity)
            # Create the 0 level directories in the HDF file
            h5createFile(HDF.File)
            for (Folder in Root.folders) {
                CreateGroups(Groups = Create_Path(Folder), File = HDF.File)
            }
            # Add the chromosome information into the metadata column
            Chrom.lengths <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = length, col.name = 'chr')
            Chrom.sizes <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = max, col.name = 'end')
            Chrom.info.df <- data.frame(chr = names(Chrom.lengths),
                nrow = as.vector(Chrom.lengths),
                size = as.vector(Chrom.sizes))
            # Create metadata chromosome groups
            Lego_WriteDataFrame(Lego = HDF.File, Path = c(Root.folders['metadata'],
                Reference.object$metadata.chrom.dataset), object = Chrom.info.df)
            # Add the chromosome information into the metadata column
            # Create matrices groups
            for (chrom1 in ChromosomeList) {
                CreateGroups(Groups = Create_Path(c(Root.folders['matrices'],chrom1)), File = HDF.File)
                for (chrom2 in ChromosomeList) {
                    Chrom2.path <- Create_Path(c(Root.folders['matrices'],chrom1,chrom2))
                    CreateGroups(Groups = Chrom2.path, File = HDF.File)
                    CreateAttributes(Path = Chrom2.path, File = HDF.File, 
                        Attributes = Reference.object$matrices.chrom.attributes)
                    Dims <-c(Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"],Chrom.info.df[Chrom.info.df$chr == chrom2,"nrow"])
                    if(is.null(ChunkSize)){
                        ChunkSize <- ceiling(Dims/100)  
                    }
                    CreateDataset(Path = c(Root.folders['matrices'],chrom1,chrom2), File = HDF.File, 
                        name = Reference.object$hdf.matrix.name, dims = Dims, maxdims = Dims)
                }
            }




            H5Fclose(HDF.Handler)
}

Lego_WriteDataFrame <- function(Lego = NULL, Path = NULL, object = NULL){
    library(stringr)
    if(!(length(c(Lego,Path,object))>=4)){
        stop("All arguments are required!")
    }
    Dirs <- Path[1:(length(Path)-1)]
    Path.to.group <- Create_Path(Dirs)
    Lego.handler <- ReturnH5Handler(Path = Path.to.group, File = Lego)
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