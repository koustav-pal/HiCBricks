% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Brick_functions.R
\name{Brick_matrix_issparse}
\alias{Brick_matrix_issparse}
\title{Check if a matrix for a chromosome pair is sparse.}
\usage{
Brick_matrix_issparse(Brick, chr1, chr2, resolution = NA)
}
\arguments{
\item{Brick}{\strong{Required}.
A string specifying the path to the Brick store created with
Create_many_Brick.}

\item{chr1}{\strong{Required}.
A character vector of length 1 specifying the chromosome corresponding to
the rows of the matrix}

\item{chr2}{\strong{Required}.
A character vector of length 1 specifying the chromosome corresponding to
the columns of the matrix}

\item{resolution}{\strong{Optional}. Default NA
When an object of class BrickContainer is provided, resolution defines the
resolution on which the function is executed}
}
\value{
Returns a logical vector of length 1, specifying if a matrix was
loaded as a sparse matrix.
}
\description{
Check if a matrix for a chromosome pair is sparse.
}
\examples{
Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
package = "HiCBricks")

out_dir <- file.path(tempdir(), "matrix_issparse_test")
dir.create(out_dir)

My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
    bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
    experiment_name = "Vignette Test", resolution = 100000,
    remove_existing = TRUE)

Matrix_file <- system.file(file.path("extdata", 
"Sexton2012_yaffetanay_CisTrans_100000_corrected_chr2L.txt.gz"), 
package = "HiCBricks")

Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr2L", 
chr2 = "chr2L", matrix_file = Matrix_file, delim = " ", 
remove_prior = TRUE, resolution = 100000)

Brick_matrix_issparse(Brick = My_BrickContainer, chr1 = "chr2L", 
chr2 = "chr2L", resolution = 100000)

}
