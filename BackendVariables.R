require(R6)
GenomicMatrix <- R6Class("GenomicMatrix",
    public = list (
        initialize = function(){},
        hdf.matrices.root = "Base.matrices",
        hdf.ranges.root = "Base.ranges",
        hdf.metadata.root = "Base.metadata",
        metadata.chrom.dataset = "chrominfo",
        GetRootFolders = function(){
            Folders <- c(self$hdf.matrices.root, self$hdf.ranges.root, self$hdf.metadata.root)
            names(Folders) <- c('matrices','ranges','metadata')
            return(Folders)
        }
        Get
        },

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
        Max.vector.size=104857600,
        Matrix.range=NA,
        Num.lines=1,
        hdf.root.folders=c("matrices","base.ranges.tables"),
        Protected.Ranges.Keys=c("Bintable"),
        Bintable.Key = "Bintable",
        NonStrandedColNames=c("chr","start","end"),
        Matrice.done = NA,
        Ranges.separator=":"
)