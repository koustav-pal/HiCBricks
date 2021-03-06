% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BrickContainer_functions.R
\name{BrickContainer_change_experiment_name}
\alias{BrickContainer_change_experiment_name}
\title{Change the location of HDF files in the BrickContainer object}
\usage{
BrickContainer_change_experiment_name(Brick = NULL, experiment_name = NULL)
}
\arguments{
\item{Brick}{\strong{Required}.
A string specifying the path to the BrickContainer created using
\code{Create_many_Bricks} or \code{Load_BrickContainer}}

\item{experiment_name}{\strong{Required}. Default NULL
A string specifying the new experiment name}
}
\value{
An object of class BrickContainer where the experiment_name
has been changed
}
\description{
\code{BrickContainer_change_experiment_name} changes the location of
name of the experiment
}
\examples{
Bintable.path <- system.file("extdata",
"Bintable_100kb.bins", package = "HiCBricks")

out_dir <- file.path(tempdir(), "BrickContainer_expt_name_test")
dir.create(out_dir)

My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
    bin_delim = " ", output_directory = out_dir, file_prefix = "Test", 
    experiment_name = "Vignette Test", resolution = 100000, 
    remove_existing = TRUE)

BrickContainer_change_experiment_name(Brick = My_BrickContainer, 
experiment_name = "I change my mind")

}
