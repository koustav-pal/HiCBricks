# HDF structure
# Base.matrices <group>
#   chromsome <group>
#       chromosome <group> Attributes FileName, Done
#           matrix <dataset>
#           extent <dataset>
#           bla bla datasets
# Base.ranges <group>
#   Bintable <group> containing binary attributes Strand, Names
#       Names <dataset>
#       ranges <dataset> based on value of attributes 3,4,5 column
#       BlaBla1 <dataset>
#       BlaBla2 <dataset>
#   Other <group>
# Base.metadata <group>
# chromosomes <dataset>

CreateLego <- function(ChromNames=NULL, BinTable=NULL, bin.delim="\t",
    col.index=c(1,2,3), impose.discontinuity=TRUE, ChunkSize=NULL, 
    Output.Filename=NULL, exec="cat", remove.existing=FALSE,
    sparse=FALSE, sparsity.compute.bins=100){
    require(rhdf5)
    H5close()
    Dir.path <- dirname(Output.Filename)
    Filename <- basename(Output.Filename)
    Working.File <- file.path(normalizePath(Dir.path),Filename)
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
    Bintable.list <- Read_bintable(Filename = BinTable, read.delim = bin.delim, exec = exec,
                col.index = col.index, chromosomes = ChromNames,
                impose.discontinuity = impose.discontinuity)
    Bintable <- Bintable.list[['main.tab']]
    # Create the 0 level directories in the HDF file
    h5createFile(HDF.File)
    for (Folder in Root.folders) {
        CreateGroups(Group.path = Create_Path(Folder), File = HDF.File)
    }
    # Add the chromosome information into the metadata column
    Chrom.lengths <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = length, col.name = 'chr')
    Chrom.sizes <- get_chrom_info(bin.table = Bintable, chrom = ChromosomeList, FUN = max, col.name = 'end')
    Chrom.info.df <- data.frame(chr = names(Chrom.lengths),
        nrow = as.vector(Chrom.lengths),
        size = as.vector(Chrom.sizes),stringsAsFactors = FALSE)
    # Create metadata chromosome groups
    Lego_WriteDataFrame(Lego = HDF.File, Path = c(Root.folders['metadata'],
            Reference.object$metadata.chrom.dataset), object = Chrom.info.df)
    Base.ranges.folders <- Reference.object$GetBaseRangesFolders()
    for (Folder in Base.ranges.folders) {
        CreateGroups(Group.path = Create_Path(c(Root.folders['ranges'],Folder)), File = HDF.File)
    }
    Lego_WriteDataFrame(Lego = HDF.File, Path = c(Root.folders['ranges'],Base.ranges.folders['Bintable'],
        Reference.object$hdf.ranges.dataset.name), object = Bintable)
    # Add the chromosome information into the metadata column
    # Create matrices groups
    for (chrom1 in ChromosomeList) {
        CreateGroups(Group.path = Create_Path(c(Root.folders['matrices'],chrom1)), File = HDF.File)
        for (chrom2 in ChromosomeList) {
            Chrom2.path <- Create_Path(c(Root.folders['matrices'],chrom1,chrom2))
            # cat(Chrom.info.df$chr,chrom1,chrom2,"\n")
            CreateGroups(Group.path = Chrom2.path, File = HDF.File)
            CreateAttributes(Path = Chrom2.path, File = HDF.File, 
                Attributes = Reference.object$matrices.chrom.attributes, on = "group")
            Dims <- c(Chrom.info.df[Chrom.info.df$chr == chrom1,"nrow"], Chrom.info.df[Chrom.info.df$chr == chrom2,"nrow"])
            if(is.null(ChunkSize)){
                ChunkSize <- ceiling(Dims/100)
            }
            CreateDataset(Path = c(Root.folders['matrices'],chrom1,chrom2), File = HDF.File, 
                name = Reference.object$hdf.matrix.name, dims = Dims, maxdims = Dims)
        }
    }
}

# _ProcessLegoInfoTable_ <- function(LOCDIR = NULL, create.dir = NULL, create.recursively = NULL, remove.prior = NULL){
#     Reference.object <- GenomicMatrix$new()
#     InfoTable <- file.path(LOCDIR,paste(Basename,Reference.object$LegosInfoTable,sep = Reference.object$GeneralFileSeparator))
#     FileName <- _GenerateRandomName_()
#     BackingFile <- paste(FileName,"backingfile.bigmem",sep = Reference.object$GeneralFileSeparator)
#     if(!dir.exists(LOCDIR) & create.dir){
#         dir.create(LOCDIR, recursive = create.recursively)
#     }else if(!dir.exists(LOCDIR)){
#         stop(LOCDIR,"not found!\n")
#     }
#     if(file.exists(InfoTable)){
#         InfoTable.read <- read.table(file = InfoTable, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#         rownames(InfoTable.read) <- paste(InfoTable.read[,Reference.object$InfoTableHeader.chr1], 
#             InfoTable.read[,Reference.object$InfoTableHeader.chr2], sep = Reference.object$InfoTableIDSeparator)
#         Chrom.id <- paste(chr1, chr2, sep = Reference.object$InfoTableIDSeparator)
#         New.df <- data.frame(chr1, chr2, FileName, BackingFile)
#         colnames(New.df) <- Reference.object$InfoTableHeaders
#         rownames(New.df) <- Chrom.id
#         if(!(Chrom.id %in% rownames(InfoTable.read))) {
#             InfoTable.read <- rbind(InfoTable.read,New.df)
#         }else{
#             if(remove.prior){
#                 Prior.file <- file.path(LOCDIR, InfoTable.read[Chrom.id,Reference.object$InfoTableHeader.name])
#                 Prior.backing.file <- file.path(LOCDIR, InfoTable.read[Chrom.id,Reference.object$InfoTableHeader.backingfile])
#                 if(file.exists(Prior.file)){
#                     file.remove(c(Prior.file,Prior.backing.file))
#                 }
#             }
#             InfoTable.read[Chrom.id,Reference.object$InfoTableHeader.name] <- FileName
#             InfoTable.read[Chrom.id,Reference.object$InfoTableHeader.backingfile] <- BackingFile
#         }
#     }
#     return(InfoTable.read)
# }

#### It should be done after the fist write 
Lego_Get_ChromInfo <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Dataset <- ._Lego_Get_Something_(Group.path = Reference.object$hdf.metadata.root, 
        Lego = Lego, Name = Reference.object$metadata.chrom.dataset, handler = FALSE)
    return(Dataset)
}

Lego_Get_Bintable <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Dataset <- ._Lego_Get_Something_(Group.path = Create_Path(c(Reference.object$hdf.ranges.root, Reference.object$hdf.bintable.ranges.group)),
        Lego = Lego, Name = Reference.object$hdf.ranges.dataset.name, handler = FALSE)
    return(Dataset)
}
Lego_List_Matrices <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Group.Handle <- ReturnH5Handler(Path = Create_Path(c(Reference.object$hdf.matrices.root)),File = Lego)
    Dataset <- h5read(name = Reference.object$metadata.chrom.dataset,
        file = Group.Handle)
    H5Gclose(Group.Handle)
    return(Dataset)
}




Lego_LoadMatrix <- function(Lego = NULL, LOCDIR = NULL, chr1 = NULL, chr2 = NULL, Basename = NULL, 
    create.dir = FALSE, create.recursively = FALSE,  exec = NULL, Dataset = NULL, 
    delim = " ", remove.prior = FALSE, num.rows = 2000, is.sparse = FALSE, sparsity.bins = 100){

    ListVars <- list(Lego = Lego, LOCDIR = LOCDIR, create.dir = create.dir, create.recursively = create.recursively, 
        Basename = Basename, chr1 = chr1, chr2 = chr2, is.sparse = is.sparse, sparsity.bins = sparsity.bins, 
        exec = exec, delim = delim, Dataset = Dataset, remove.prior = remove.prior)
    sapply(1:length(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    sapply(1:length(ListVars[c("Lego","LOCDIR","chr1","chr2","Basename","Dataset")]),function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })
    Chrom.info.df <- Lego_GetChromInfo(Lego = Lego)
    if(!(all(c(chr1, chr2) %in% Chrom.info.df[,"chr"]))){
        stop("Provided chromosomes do not exist in the chrom table\n")
    }
    Chrom1.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr1,"nrow"]
    Chrom2.len <- Chrom.info.df[Chrom.info.df[,"chr"]==chr2,"nrow"]
    _ProcessMatrix_(Data = Dataset, delim = NULL, Matrix.file = Matrix.file.path, 
        exec = exec, chr1.len = Chrom1.len, chr2.len = Chrom2.len, 
        fix.num.rows.at = num.rows, is.sparse = is.sparse, sparsity.bins = sparsity.bins)
}

Lego_ListRangeKeys <- function(Lego = NULL){
    Reference.object <- GenomicMatrix$new()
    Handler <- ReturnH5Handler(Path = Create_Path(Reference.object$hdf.ranges.root),File = Lego)
    GroupList <- h5ls(Handler, datasetinfo = FALSE, recursive = FALSE)
    return(GroupList)
}



Lego_GetRanges <- function(Lego = NULL, chr = NULL, rangekey = NULL){
    Reference.object <- GenomicMatrix$new()
    if(is.null(chr)){
        GetAll <- TRUE
    }
    if(is.null(rangekey) | is.null(Lego)){
        stop("rangekey and Lego cannot remain empty!\n")
    }

}


Lego_attach_mcol <- function()






Lego_LoadMatrix <- function(Lego = NULL, chr1 = NULL, chr2 = NULL, FormatAs = "mxnMatrix", 
    delim = " ", Dataset = NULL, exec = NULL){
    
    ListVars <- list(Lego = Lego, FormatAs = FormatAs, chr1 = chr1, chr2 = chr2, exec = NULL, delim = delim, Dataset = Dataset)
    sapply(1:length(ListVars),function(x){
        if(length(ListVars[[x]]) > 1){
            stop(names(ListVars[x]),"had length greater than 1.\n")
        }
    })
    sapply(1:length(ListVars[c("Lego","chr1","chr2","Dataset")]),function(x){
        if(is.null(ListVars[[x]])){
            stop(names(ListVars[x]),"has no value.\n")
        }
    })

    # Process and then return 0 
    
}


Lego_WriteDataFrame <- function(Lego = NULL, Path = NULL, object = NULL){
    library(stringr)
    if(!(length(c(Lego,Path,object))>=4)){
        stop("All arguments are required!")
    }
    Dirs <- Path[1:(length(Path)-1)]
    Path.to.group <- Create_Path(Dirs)
    cat(Path.to.group,Path[length(Path)],"\n")
    Lego.handler <- ReturnH5Handler(Path = Path.to.group, File = Lego)
    h5writeDataset.data.frame(h5loc = Lego.handler, 
        obj = object, 
        name = Path[length(Path)])
    H5Gclose(Lego.handler)
}


# AddAttribute(Key=private$Bintable.Key,Value=normalizePath(BinTable))


# if(sparse){
#     if(length(sparsity.compute.bins)>1 | !is.numeric(sparsity.compute.bins)){
#         stop("sparsity.compute.bins expects numeric value of length 1\n")
#     }
#     private$ComputeSparsity=sparse
#     private$Sparsity.compute.bins=sparsity.compute.bins
# }