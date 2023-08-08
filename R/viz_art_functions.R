#' Create the entire HDF5 structure and load the bintable
#' 
#' `Brick_vizart_plot_heatmap` creates various heatmaps and plots TADs. 
#' 
#' This function provides the capability to plot various types of heatmaps from
#' Hi-C data. 
#' \itemize{
#'  \item One sample heatmap.
#'  \item Two sample heatmap (One sample on upper and other on lower).
#'  \item All of the above with 90 degree rotation.
#'  \item All of the above but with signal capped at a certain value.
#'  \item All of the above but filtered by distance.
#'  \item All of the above with TADs/TAD borders plotted on top.
#' }
#' 
#' @inheritParams Brick_add_ranges
#' 
#' @param File \strong{Required}
#' A character vector containing the output filename to write.
#' 
#' @param Bricks \strong{Required}
#' A list of length 1 (in case of one sample heatmaps) or 2 (in case
#' of two sample heatmaps) specifying the BrickContainers from where
#' to fetch the data.
#' 
#' @param x_coords \strong{Required}
#' A character vector of length 1 specifying the coordinates from where to fetch
#' the data.
#' 
#' @param y_coords \strong{Required}
#' A character vector of length 1 specifying the coordinates from where to fetch
#' the data.
#' 
#' @param FUN \strong{Optional}. Default NULL
#' If any sort of transformations should be applied to the data before plotting.
#' Such as, log10 or log2 transformations. 
#' 
#' @param value_cap \strong{Optional}. Default NULL
#' If present, values beyond a certain quantile will be capped to that quantile.
#' In Hi-C this helps to emphasize structural information. Please note, if this
#' parameter is present the greatest value will have a greater than sign append-
#' -ed to them. 
#' 
#' @param distance \strong{Optional}. Default NULL
#' If present, values beyond this distance will be filtered out. Please note,
#' that if a Brick store matrix was loaded until a certain distance, this 
#' parameter will result in an error if it is greater than the loaded distance.
#' 
#' @param rotate \strong{Optional}. Default FALSE
#' If TRUE, will rotate the heatmap by 90 degrees.
#' 
#' @param x_axis \strong{Optional}. Default TRUE
#' If FALSE, the x-axis will be removed (ticks, x-axis labels and title).
#'  
#' @param x_axis_title \strong{Optional}. Default NULL
#' If present, will be the \emph{x-axis} title. Else defaults to the provided
#' x_coords
#' 
#' @param y_axis \strong{Optional}. Default TRUE
#' If FALSE, the y-axis will be removed (ticks, y-axis labels and title).
#' 
#' @param y_axis_title \strong{Optional}. Default NULL
#' If present, will be the \emph{y-axis} title. Else defaults to the provided
#' y_coords
#' 
#' @param title \strong{Optional}. Default NULL
#' If present, will be the \emph{plot} title. Else defaults to the provided
#' x_coords vs y_coords
#' 
#' @param legend_title \strong{Optional}. Default NULL
#' If present will be the title of the legend. Else defaults to "Signal".
#' 
#' @param return_object \strong{Optional}. Default FALSE
#' If present the ggplot object will be returned
#' 
#' @param x_axis_num_breaks \strong{Optional}. Default 5
#' Number of ticks on the x axis
#' 
#' @param y_axis_num_breaks \strong{Optional}. Default 5
#' Number of ticks on the y axis
#' 
#' @param x_axis_text_size \strong{Optional}. Default 10
#' x-axis text size
#' 
#' @param y_axis_text_size \strong{Optional}. Default 10
#' y-axis text size
#' 
#' @param text_size \strong{Optional}. Default 10
#' text size of text elements in the plot.
#' 
#' @param legend_title_text_size \strong{Optional}. Default 8
#' text size of the legend title
#' 
#' @param legend_text_size \strong{Optional}. Default 8
#' text size of the legend text
#' 
#' @param title_size \strong{Optional}. Default 10
#' text size of the title
#' 
#' @param tad_ranges \strong{Optional}. Default NULL
#' A GenomicRanges object specifying the start and end coordinates of TADs to be
#' plotted on the heatmap.
#' 
#' @param group_col \strong{Optional}. Default NULL
#' Name of the column which will be used to categorize TADs as belonging to 
#' either the first or the second Brick stores. This must be a numeric value
#' ranging from 1 to 2. If NULL, TADs will be plotted on both Hi-C maps.
#' 
#' @param tad_colour_col \strong{Optional}. Default NULL
#' tad_colour_col takes as value the column name in the tad_ranges object 
#' corresponding to the column which should be used to define different TAD
#' categories.
#' 
#' @param line_width \strong{Optional}. Default 0.5
#' When plotting TADs set the width of the plotted lines
#' 
#' @param cut_corners \strong{Optional}. Default FALSE
#' if cut_corners is TRUE, TAD borders will not be truncated, and they will
#' span until the end of visible heatmap.
#' 
#' @param highlight_points \strong{Optional}. Not yet implemented.
#' 
#' @param colours \strong{Optional}. Default NULL
#' If tad_ranges is present, colours expects a hexcode value of length 1. But,
#' if tad_colour_col is specified, it expects colours of the same length as 
#' unique tad_ranges$tad_colour_col. 
#' 
#' @param colours_names \strong{Optional}. Default NULL
#' If present, will be assigned to colours. Else, will inherit unique 
#' tad_colour_col. If tad_colour_col is also absent, will revert to a placehold 
#' column name.
#' 
#' @param palette \strong{Required}. Default NULL
#' One of the RColorbrewer or viridis colour palettes
#' 
#' @param col_direction \strong{Optional}. Default 1
#' If -1, the colour scale will be reversed.
#' 
#' @param extrapolate_on \strong{Optional}. Default NULL
#' If present, colours from the palette will be extrapolated between lightest
#' and darkest to create the gradient. This value cannot be more than 100.
#' 
#' @param width \strong{Optional}. Default 10cm
#' Width of the output file units. 
#' 
#' @param height \strong{Optional}. Default 6cm
#' Height of the output file in units.
#' 
#' @param units \strong{Optional}. Default cm
#' Defines the units of the output file width and height.
#' 
#' @param legend_key_width \strong{Optional}. Default unit(3,"cm")
#' Defines the legend key width.
#' 
#' @param legend_key_height \strong{Optional}. Default unit(0.5,"cm")
#' Defines the legend key height.
#' 
#' @return If return_object is set to TRUE, the constructed ggplot2 
#' object will be returned. Else TRUE.
#' 
#' @examples
#' FailSafe_log10 <- function(x){
#'      x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
#'      return(log10(x+1))
#' }
#' 
#' Bintable.path <- system.file(file.path("extdata", "Bintable_100kb.bins"), 
#' package = "HiCBricks")
#' 
#' out_dir <- file.path(tempdir(), "vizart_test")
#' dir.create(out_dir)
#' 
#' My_BrickContainer <- Create_many_Bricks(BinTable = Bintable.path, 
#'   bin_delim = " ", output_directory = out_dir, file_prefix = "Test",
#'   experiment_name = "Vignette Test", resolution = 100000,
#'   remove_existing = TRUE)
#' 
#' Matrix_file <- system.file(file.path("extdata", 
#' "Sexton2012_yaffetanay_CisTrans_100000_corrected_chr3R.txt.gz"), 
#' package = "HiCBricks")
#' 
#' Brick_load_matrix(Brick = My_BrickContainer, chr1 = "chr3R", 
#' chr2 = "chr3R", matrix_file = Matrix_file, delim = " ",
#' remove_prior = TRUE, resolution = 100000)
#' 
#' Brick_vizart_plot_heatmap(File = "./chr3R-1-10000000.pdf", 
#' Bricks = list(My_BrickContainer), resolution = 100000, 
#' x_coords = "chr3R:1:10000000", palette = "Reds", 
#' y_coords = "chr3R:1:10000000", FUN = FailSafe_log10, 
#' value_cap = 0.99, width = 10, height = 11, legend_key_width = unit(3,"mm"),
#' legend_key_height = unit(0.3,"cm"))
#' 
Brick_vizart_plot_heatmap <- function (File, Bricks, resolution,
    x_coords, y_coords, invert_coords = FALSE, FUN = NULL, value_cap = NULL,
    distance = NULL, rotate = FALSE, x_axis = TRUE,  x_axis_title = NULL,
    y_axis = TRUE, y_axis_title = NULL, title = NULL, legend_title = NULL,
    return_object = FALSE, x_axis_num_breaks = 5, y_axis_num_breaks = 5,
    palette, col_direction = 1, extrapolate_on = NULL,
    x_axis_text_size = 10, y_axis_text_size = 10, text_size = 10,
    legend_title_text_size = 8, legend_text_size = 8, title_size = 10,
    tad_ranges = NULL, group_col = NULL, tad_colour_col = NULL, colours = NULL,
    colours_names = NULL, cut_corners = FALSE, highlight_points = NULL,
    width = 10, height = 6, line_width = 0.5, units = "cm",
    legend_key_width = unit(3, "cm"), legend_key_height = unit(0.5,"cm")){
    if (invert_coords == TRUE){
        new_x_coords <- y_coords
        new_y_coords <- x_coords
        x_coords <- new_x_coords
        y_coords <- new_y_coords
    }
    if (!is.list(Bricks)) {
        stop("Bricks expects an argument of type list.",
             " Please refer to the vignette to understand the parameter.")
    }
    Matrix.df <- HiCBricks:::Get_one_or_two_brick_regions(Bricks = Bricks,
        resolution = resolution, x_coords = x_coords, y_coords = y_coords,
        distance = distance, value_cap = value_cap, FUN = FUN)
    if (nrow(Matrix.df) == 0) {
        stop("The matrix was empty!")
    }
    list_of_coords <- list(x_coords = x_coords, y_coords = y_coords)
    Parsed_string <- HiCBricks:::._Parse_genomic_coordinates(list_of_coords)
    x.coord.parsed <- Parsed_string[["x_coords"]]
    y.coord.parsed <- Parsed_string[["y_coords"]]
    x.coord.breaks <- HiCBricks:::make_axis_coord_breaks(from = min(Matrix.df$row),
        to = max(Matrix.df$row), how.many = x_axis_num_breaks,
        two.sample = FALSE)
    x_axis.coord.labs <- HiCBricks:::Make_axis_labels(Brick = Bricks[[1]],
        resolution = resolution, chr = x.coord.parsed["chr"],
        positions = x.coord.breaks)
    two.sample <- (rotate & length(Bricks) == 2)
    y.coord.breaks <- HiCBricks:::make_axis_coord_breaks(from = min(Matrix.df$col),
        to = max(Matrix.df$col), how.many = y_axis_num_breaks,
        two.sample = two.sample)
    y_axis.coord.labs <- HiCBricks:::Make_axis_labels(Brick = Bricks[[1]],
        resolution = resolution, chr = y.coord.parsed["chr"],
        positions = abs(y.coord.breaks))
    Colours <- HiCBricks:::Make_colours(palette = palette, extrapolate_on = extrapolate_on,
        direction = col_direction)
    two.sample <- (length(Bricks) == 2)
    Matrix.df$rescale <- HiCBricks:::rescale_values_for_colours(Object = Matrix.df,
        two.sample = two.sample)
    Value.dist <- HiCBricks:::make_colour_breaks(Object = Matrix.df, how.many = length(Colours),
        two.sample = two.sample)
    Legend.breaks.list <- HiCBricks:::get_legend_breaks(Object = Matrix.df,
        how.many = 5, value_cap = value_cap, colours = Colours,
        two.sample = two.sample)
    Colour.breaks <- Legend.breaks.list[["col.breaks"]]
    Colour.labs <- Legend.breaks.list[["col.labs"]]
    Colours <- Legend.breaks.list[["cols"]]
    if (rotate) {
        y.coord.breaks <- y.coord.breaks - min(y.coord.breaks)
        x.coord.breaks <- x.coord.breaks - min(x.coord.breaks)
        if (length(Bricks) == 2) {
            Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0, ]
            Lower.tri.map <- Matrix.df[Matrix.df$dist <= 0, ]
            Lower.tri.map$dist <- abs(Lower.tri.map$dist)
            Upper.rotated.map <- RotateHeatmap(Matrix = Upper.tri.map,
                value.var = "rescale", upper = TRUE)
            Lower.rotated.map <- RotateHeatmap(Matrix = Lower.tri.map,
                value.var = "rescale", upper = FALSE)
            Entire.rotated.map <- rbind(Upper.rotated.map, Lower.rotated.map)
            y.coord.breaks <- y.coord.breaks/2
            y.coord.breaks <- c(rev(y.coord.breaks) * -1, y.coord.breaks)
            y_axis.coord.labs <- c(rev(y_axis.coord.labs), y_axis.coord.labs)
        } else {
            Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0, ]
            Entire.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Upper.tri.map,
                value.var = "rescale", upper = TRUE)
            y.coord.breaks <- y.coord.breaks/2
        }
    }
    Brick_theme <- HiCBricks:::Get_heatmap_theme(x_axis = x_axis, y_axis = y_axis,
        text_size = text_size, x_axis_text_size = x_axis_text_size,
        y_axis_text_size = y_axis_text_size, legend_title_text_size = legend_title_text_size,
        legend_text_size = legend_text_size, title_size = title_size,
        legend_key_width = legend_key_width, legend_key_height = legend_key_height)
    Labels <- HiCBricks:::Get_heatmap_titles(title = title, x_axis_title = x_axis_title,
        y_axis_title = y_axis_title, legend_title = legend_title,
        x_coords = x_coords, y_coords = y_coords, rotate = rotate)
    Boundaries.obj <- NULL
    if (!is.null(tad_ranges)) {
        Boundaries.obj <- HiCBricks:::Format_boundaries_normal_heatmap(Bricks = Bricks,
            resolution = resolution, Ranges = tad_ranges, group_col = group_col,
            cut_corners = cut_corners, colour.col = tad_colour_col,
            colours = colours, colours_names = colours_names,
            region.chr = x.coord.parsed["chr"], region.start = as.numeric(x.coord.parsed["start"]),
            region.end = as.numeric(x.coord.parsed["end"]), distance = distance,
            rotate = rotate)
    }
    if (rotate) {
        ids <- xcoords <- ycoords <- NULL
        if (invert_coords == TRUE) {
            ThePlot <- ggplot(Entire.rotated.map, aes(x = ycoords,
                y = xcoords))
        } else {
            ThePlot <- ggplot(Entire.rotated.map, aes(x = xcoords,
                y = ycoords))
        }
        ThePlot <- ThePlot + geom_polygon(aes(fill = values,
            group = ids))
        xlims <- c(0, max(Entire.rotated.map[, "xcoords"]))
        ylims <- c(min(Entire.rotated.map[, "ycoords"]), max(Entire.rotated.map[,
            "ycoords"]))
        y.coord.breaks <- seq(ceiling(min(Entire.rotated.map[,
            "ycoords"])), ceiling(max(Entire.rotated.map[, "ycoords"])),
            length.out = y_axis_num_breaks)
        y_axis.coord.labs <- y.coord.breaks * 2
    } else {
        Matrix.df$row <- Matrix.df$row - 0.5
        Matrix.df$col <- Matrix.df$col - 0.5
        if (invert_coords == TRUE) {
            ThePlot <- ggplot(Matrix.df, aes(x = col, y = row))
        } else {
            ThePlot <- ggplot(Matrix.df, aes(x = row, y = col))
        }
        ThePlot <- ThePlot + geom_tile(aes(fill = rescale))
        xlims <- c(min(Matrix.df$row) - 0.5, max(Matrix.df$row) +
            0.5)
        ylims <- c(min(Matrix.df$col) - 0.5, max(Matrix.df$col) +
            0.5)
    }
    if (!is.null(tad_ranges)) {
        if (invert_coords == TRUE) {
            line.group <- y <- x <- NULL
            ThePlot <- ThePlot + geom_line(data = Boundaries.obj,
                aes(x = y, y = x, group = line.group, colour = colours),
                size = line_width)
        } else {
            line.group <- x <- y <- NULL
            ThePlot <- ThePlot + geom_line(data = Boundaries.obj,
                aes(x = x, y = y, group = line.group, colour = colours),
                size = line_width)
        }
        ThePlot <- ThePlot + scale_colour_manual(values = colours)
    }
    if (invert_coords == TRUE) {
        ThePlot <- ThePlot + scale_x_continuous(limits = ylims, expand = c(0,
            0), breaks = y.coord.breaks, labels = y_axis.coord.labs)
        ThePlot <- ThePlot + scale_y_continuous(limits = xlims, expand = c(0,
            0), breaks = x.coord.breaks, labels = x_axis.coord.labs)
    } else {
        ThePlot <- ThePlot + scale_x_continuous(limits = xlims, expand = c(0,
            0), breaks = x.coord.breaks, labels = x_axis.coord.labs)
        ThePlot <- ThePlot + scale_y_continuous(limits = ylims, expand = c(0,
            0), breaks = y.coord.breaks, labels = y_axis.coord.labs)
    }
    ThePlot <- ThePlot + scale_fill_gradientn(legend_title, values = Value.dist,
        breaks = Colour.breaks, labels = Colour.labs, colors = Colours)
    ThePlot <- ThePlot + Brick_theme
    if (invert_coords == TRUE) {
        ThePlot <- ThePlot + labs(title = Labels["title"], x = Labels["y_axis"],
            y = Labels["x_axis"])
    } else {
        ThePlot <- ThePlot + labs(title = Labels["title"], x = Labels["x_axis"],
            y = Labels["y_axis"])
    }
    ggsave(filename = File, plot = ThePlot, width = width, height = height,
        units = units)
    if (return_object) {
        return(ThePlot)
    } else {
        return(TRUE)
    }
}
