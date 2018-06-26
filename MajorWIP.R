# HDF structure
# Base.matrices <group>
# 	chromsome <group>
# 		chromosome <group> Attributes FileName, Done
#			matrix <dataset>		
# 			extent <dataset>
# 			bla bla datasets
# Base.ranges <group>
# 	Names <dataset> containing the names of all the ranges tables listed under Base.ranges 
#	Bintable <group> containing binary attributes Strand, Names
#		Names <dataset>
#		ranges <dataset> based on value of attributes 3,4,5 column 
#		BlaBla1 <dataset>
# 		BlaBla2	<dataset>
# Base.metadata <group> 
# chromosomes <dataset> 





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

