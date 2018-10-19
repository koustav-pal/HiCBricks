\name{HiCBlocks}
\alias{HiCBlocks}
\title{A package for storing, accessing and plotting Hi-C data}
\description{
HiCBlocks is a package allowing users to flexibly import and work with Hi-C data
}
\details{
Using HiCBlocks users are able to import Hi-C matrices stored in various 
formats into an HDF structure. This is the Block file. You can then access the 
Hi-C data using accessor functions. Since the data is stored in an HDF file, 
if you have the Block (HDF) file, you can keep on accessing the same file an 
infinite number of times.

Users can also associate different ranges objects with the HDF file.

The HDF file must have the same structure as followed by HiCBlocks

Users can then move forward and create analysis pipelines and statistical 
methods based on HiCBlocks HDF files without worrying about the underlying 
data structure. To showcase this, Local score differentiator (LSD) our novel 
TAD calling procedure comes packaged with HiCBlocks. 

You are also able to plot Hi-C data using HiCBlocks functions. There are a few
types. You can create, 
\itemize{
\item a square heatmap
\item a rotated heatmap
\item two group square/rotated heatmaps
\item both heatmaps until a certain distance
\item plot TADs on both heatmaps
}
}
\section{Block creation}{
\itemize{
\item \code{\link{CreateBlock}} - Create the complete HDF data structure. We 
refer to the HDF files as Block
\item \code{\link{CreateBlock_from_mcool}} - Create the complete Block data 
structure from an mcool file.
}
}
\section{Matrix loaders}{
\itemize{
\item \code{\link{Block_load_matrix}} - Load a complete nxm dimensional matrix.
\item \code{\link{Block_load_cis_matrix_till_distance}} - Load a sam chromosome
nxn dimensional matrix until a certain distance. 
\item \code{\link{Block_load_data_from_mcool}} - Load parts of the data from 
the 4DN consortium generated mcool files. 
}
}
\section{Matrix Accessors}{
\itemize{
\item \code{\link{Block_get_matrix_within_coords}} - Fetches a matrix within the
provided genomic coordinates.
\item \code{\link{Block_get_matrix}} - Fetches a matrix within the provided x
and y coordinates.
\item \code{\link{Block_get_values_by_distance}} - Fetch all values 
corresponding to interactions between genomic loci separated by the
corresponding value.
\item \code{\link{Block_fetch_row_vector}} - Fetch all values at a given row
or column.
}
All of the functions above can be subsetted and contain further value
transformations.
}
\section{Ranges operators}{
\itemize{
\item \code{\link{Block_get_bintable}} - All HiCBlocks Block files contain a 
binning table containing the coordinate information of the matrix. This fetches
the associated binning table.
\item \code{\link{Block_add_ranges}} - Add a ranges object to the Block file.
\item \code{\link{Block_get_ranges}} - Get a ranges object associated to a Block
file.
\item \code{\link{Block_fetch_range_index}} - Provided a set of coordinate 
vectors, get the corresponding rows/cols overlapping with those coordinates.
\item \code{\link{Block_make_ranges}} - Create a granges object from provided 
vectors.
\item \code{\link{Block_return_region_position}} - Get the row/col number
corresponding to coordinates spelled out in human readable format.
}
}
\section{Other functions}{
\itemize{
\item \code{\link{Block_local_score_differentiator}} - Use the LSD TAD calling
procedure to do some TAD calls.
\item \code{\link{Block_vizart_plot_heatmap}} - Plot pretty heatmaps.
}
}
\section{Utility functions}{
\itemize{
\item \code{\link{Block_get_chrominfo}} - Get the basic information regarding
the Block file. Which chromosomes are present, dimension of the matrix and the
total length of the chromosome.
\item \code{\link{Block_get_matrix_mcols}} - Get the matrix metadata 
information. Such as, row sums, coverage information and how sparse regions
near the diagonal are.
\item \code{\link{Block_list_matrices}} - List all the matrices present in
the Block file. Alongside, also provide information such as if the matrix has
been loaded or not, min max values, e.t.c
\item \code{\link{Block_list_rangekeys}} - List the names of the ranges present
in the Block file.
\item \code{\link{Block_rangekey_exists}} - Answers the question, is this 
rangekey present in the Block file?
\item \code{\link{Block_list_ranges_mcols}} - List the names of metadata columns
associated to a ranges object in the Block file.
\item \code{\link{Block_matrix_dimensions}} - Get the dimensions of a given 
matrix.
\item \code{\link{Block_matrix_exists}} - Answers the question, has a matrix 
been created for this Block store?
\item \code{\link{Block_matrix_filename}} - Answers the question, what is the 
name of the file used to load this particular matrix?
\item \code{\link{Block_matrix_isdone}} - Answers the question, has this matrix
been loaded already?
\item \code{\link{Block_matrix_issparse}} - Answers the question, was this 
matrix defined as a sparse matrix while loading?
\item \code{\link{Block_matrix_maxdist}} - If 
\code{\link{Block_load_cis_matrix_till_distance}} was used for loading data, 
then this function will tell you until what distance data was loaded.
\item \code{\link{Block_matrix_minmax}} - Outputs the value range of the matrix.
}
}
\section{mcool utility functions}{
\itemize{
\item \code{\link{Block_list_mcool_normalisations}} - List the names of
normalisation vectors that can be present in a mcool file.
\item \code{\link{Block_mcool_normalisation_exists}} - Check if a specific
normalisation vector exists in an mcool file.
\item \code{\link{Block_list_mcool_resolutions}} - List the resolutions present 
in an mcool file.
}
}
\section{Bioc Cache functions}{
\itemize{
\item \code{\link{Block_is_tracked}} - Check if a Block is being tracked
\item \code{\link{Block_untrack_block}} - Untrack the Block file at path, if
it is being tracked
\item \code{\link{Block_track_blocks}} - Track the Block file at path
\item \code{\link{Block_list_tracked_blocks}} - List all Blocks that are being
tracked by the HiCBlocks package.
\item \code{\link{Block_list_tracked_blocks}} - List all Blocks that are being
tracked by the HiCBlocks package.
\item \code{\link{Block_path_to_file}} - Returns the path to a file that is 
being tracked.
}
}




