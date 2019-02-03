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
#' @param File \strong{Required}
#' A character vector containing the output filename to write.
#' 
#' @param Bricks \strong{Required}
#' A character vector of length 1 (in case of one sample heatmaps) or 2 (in case
#' of two sample heatmaps) specifying the names of the Brick stores from where
#' to fetch the data.
#' 
#' @param x.coords \strong{Required}
#' A character vector of length 1 specifying the coordinates from where to fetch
#' the data.
#' 
#' @param y.coords \strong{Required}
#' A character vector of length 1 specifying the coordinates from where to fetch
#' the data.
#' 
#' @param FUN \strong{Optional}. Default NULL
#' If any sort of transformations should be applied to the data before plotting.
#' Such as, log10 or log2 transformations. 
#' 
#' @param value.cap \strong{Optional}. Default NULL
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
#' @param x.axis \strong{Optional}. Default TRUE
#' If FALSE, the x-axis will be removed (ticks, x-axis labels and title).
#'  
#' @param x.axis.title \strong{Optional}. Default NULL
#' If present, will be the \emph{x-axis} title. Else defaults to the provided
#' x.coords
#' 
#' @param y.axis \strong{Optional}. Default TRUE
#' If FALSE, the y-axis will be removed (ticks, y-axis labels and title).
#' 
#' @param y.axis.title \strong{Optional}. Default NULL
#' If present, will be the \emph{y-axis} title. Else defaults to the provided
#' y.coords
#' 
#' @param title \strong{Optional}. Default NULL
#' If present, will be the \emph{plot} title. Else defaults to the provided
#' x.coords vs y.coords
#' 
#' @param legend.title \strong{Optional}. Default NULL
#' If present will be the title of the legend. Else defaults to "Signal".
#' 
#' @param return.object \strong{Optional}. Default FALSE
#' If present the ggplot object will be returned
#' 
#' @param x.axis.num.breaks \strong{Optional}. Default 5
#' Number of ticks on the x axis
#' 
#' @param y.axis.num.breaks \strong{Optional}. Default 5
#' Number of ticks on the y axis
#' 
#' @param x.axis.text.size \strong{Optional}. Default 10
#' x-axis text size
#' 
#' @param y.axis.text.size \strong{Optional}. Default 10
#' y-axis text size
#' 
#' @param text.size \strong{Optional}. Default 10
#' text size of text elements in the plot.
#' 
#' @param legend.title.text.size \strong{Optional}. Default 8
#' text size of the legend title
#' 
#' @param legend.text.size \strong{Optional}. Default 8
#' text size of the legend text
#' 
#' @param title.size \strong{Optional}. Default 10
#' text size of the title
#' 
#' @param tad.ranges \strong{Optional}. Default NULL
#' A GenomicRanges object specifying the start and end coordinates of TADs to be
#' plotted on the heatmap.
#' 
#' @param group.col \strong{Optional}. Default NULL
#' Name of the column which will be used to categorize TADs as belonging to 
#' either the first or the second Brick stores. This must be a numeric value
#' ranging from 1 to 2. If NULL, TADs will be plotted on both Hi-C maps.
#' 
#' @param tad.colour.col \strong{Optional}. Default NULL
#' tad.colour.col takes as value the column name in the tad.ranges object 
#' corresponding to the column which should be used to define different TAD
#' categories.
#' 
#' @param line.width \strong{Optional}. Default 0.5
#' When plotting TADs set the width of the plotted lines
#' 
#' @param cut.corners \strong{Optional}. Default FALSE
#' if cut.corners is TRUE, TAD borders will not be truncated, and they will
#' span until the end of visible heatmap.
#' 
#' @param highlight.points \strong{Optional}. Not yet implemented.
#' 
#' @param colours \strong{Optional}. Default NULL
#' If tad.ranges is present, colours expects a hexcode value of length 1. But,
#' if tad.colour.col is specified, it expects colours of the same length as 
#' unique tad.ranges$tad.colour.col. 
#' 
#' @param colours.names \strong{Optional}. Default NULL
#' If present, will be assigned to colours. Else, will inherit unique 
#' tad.colour.col. If tad.colour.col is also absent, will revert to a placehold 
#' column name.
#' 
#' @param palette \strong{Required}. Default NULL
#' One of the RColorbrewer or viridis colour palettes
#' 
#' @param col.direction \strong{Optional}. Default 1
#' If -1, the colour scale will be reversed.
#' 
#' @param extrapolate.on \strong{Optional}. Default NULL
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
#' @param legend.key.width \strong{Optional}. Default unit(3,"cm")
#' Defines the legend key width.
#' 
#' @param legend.key.height \strong{Optional}. Default unit(0.5,"cm")
#' Defines the legend key height.
#' 
#' @return If return.object is set to TRUE, the constructed ggplot2 
#' object will be returned. Else TRUE.
#' 
#' @examples
#' FailSafe_log10 <- function(x){
#'      x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
#'      return(log10(x+1))
#' }
#' 
#' Brick.file <- system.file("extdata", "test.hdf", package = "HiCBricks")
#' Brick_vizart_plot_heatmap(File = "./chr19-5000000-10000000.pdf", 
#' Bricks = Brick.file, x.coords = "chr19:5000000:10000000", palette = "Reds",
#' y.coords = "chr19:5000000:10000000", FUN = FailSafe_log10, 
#' value.cap = 0.99, width = 10, height = 11, legend.key.width = unit(3,"mm"),
#' legend.key.height = unit(0.3,"cm"))
#' 
Brick_vizart_plot_heatmap = function(File, Bricks, 
    x.coords, y.coords, FUN = NULL, value.cap = NULL, 
    distance = NULL, rotate = FALSE, x.axis = TRUE, x.axis.title = NULL, 
    y.axis = TRUE, y.axis.title = NULL, title = NULL, legend.title = NULL, 
    return.object=FALSE, x.axis.num.breaks = 5, y.axis.num.breaks = 5, 
    palette, col.direction = 1, extrapolate.on = NULL, 
    x.axis.text.size = 10, y.axis.text.size = 10, text.size = 10,
    legend.title.text.size = 8, legend.text.size = 8, title.size = 10,
    tad.ranges = NULL, group.col = NULL, tad.colour.col = NULL, colours = NULL,
    colours.names = NULL, cut.corners = FALSE, highlight.points = NULL, 
    width = 10, height = 6, line.width = 0.5, units = "cm", 
    legend.key.width = unit(3,"cm"), legend.key.height = unit(0.5,"cm")){

    Matrix.df <- Get_one_or_two_brick_regions(Bricks = Bricks, 
        x.coords = x.coords, y.coords = y.coords, distance = distance,
        value.cap = value.cap, FUN = FUN)
    if(nrow(Matrix.df)==0){
        stop("The matrix was empty!")
    }
    list_of_coords <- list("x.coords" = x.coords, "y.coords" = y.coords)

    Parsed_string <- ._Parse_genomic_coordinates(list_of_coords)
    x.coord.parsed <- Parsed_string[["x.coords"]]
    y.coord.parsed <- Parsed_string[["y.coords"]]
    x.coord.breaks <- make_axis_coord_breaks(from = min(Matrix.df$row), 
        to = max(Matrix.df$row), how.many = x.axis.num.breaks, 
        two.sample = FALSE)
    x.axis.coord.labs <- Make_axis_labels(Brick = Bricks[1], 
        chr = x.coord.parsed['chr'], positions = x.coord.breaks)

    two.sample <- (rotate & length(Bricks)==2)
    y.coord.breaks <- make_axis_coord_breaks(from = min(Matrix.df$row), 
        to = max(Matrix.df$row), how.many = x.axis.num.breaks, 
        two.sample = two.sample)
    y.axis.coord.labs <- Make_axis_labels(Brick = Bricks[1], 
        chr = y.coord.parsed['chr'], positions = abs(y.coord.breaks))

    Colours <- Make_colours(palette = palette, 
        extrapolate.on = extrapolate.on, direction = col.direction)
    # go from min val to mid to max val
    two.sample <- (length(Bricks)==2)
    Matrix.df$rescale <- rescale_values_for_colours(
        Object = Matrix.df, two.sample = two.sample)
    Value.dist <- make_colour_breaks(Object = Matrix.df, 
        how.many = length(Colours), two.sample = two.sample)
    
    Legend.breaks.list <- get_legend_breaks(Object = Matrix.df, 
        how.many = 5, value.cap = value.cap, colours = Colours,
        two.sample = two.sample)
    Colour.breaks <- Legend.breaks.list[["col.breaks"]]
    Colour.labs <- Legend.breaks.list[["col.labs"]]
    Colours <- Legend.breaks.list[["cols"]]
    if(rotate){
        y.coord.breaks <- y.coord.breaks - min(y.coord.breaks)
        x.coord.breaks <- x.coord.breaks - min(x.coord.breaks)
        if(length(Bricks)==2){
            Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0,]
            Lower.tri.map <- Matrix.df[Matrix.df$dist <= 0,]
            Lower.tri.map$dist <- abs(Lower.tri.map$dist)
            Upper.rotated.map <- RotateHeatmap(Matrix=Upper.tri.map, 
                value.var="rescale", upper = TRUE)
            Lower.rotated.map <- RotateHeatmap(Matrix=Lower.tri.map, 
                value.var="rescale", upper = FALSE) 
            Entire.rotated.map <- rbind(Upper.rotated.map,Lower.rotated.map)
            y.coord.breaks <- y.coord.breaks/2
            y.coord.breaks <- c(rev(y.coord.breaks)*-1,y.coord.breaks)
            y.axis.coord.labs <- c(rev(y.axis.coord.labs),y.axis.coord.labs)
        }else{
            Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0,]
            Entire.rotated.map <- RotateHeatmap(Matrix=Upper.tri.map, 
                value.var="rescale", upper = TRUE)
            y.coord.breaks <- y.coord.breaks/2
        }
    }

    # require(ggplot2)
    Brick_theme <- Get_heatmap_theme(x.axis=x.axis, y.axis=y.axis, 
        text.size = text.size, x.axis.text.size = x.axis.text.size, 
        y.axis.text.size = y.axis.text.size, 
        legend.title.text.size = legend.title.text.size, 
        legend.text.size = legend.text.size, 
        title.size = title.size, legend.key.width = legend.key.width, 
        legend.key.height =legend.key.height)
    Labels <- Get_heatmap_titles(title = title, x.axis.title = x.axis.title, 
        y.axis.title = y.axis.title, legend.title = legend.title, 
        x.coords = x.coords, y.coords = y.coords, rotate = rotate)
    Boundaries.obj <- NULL
    if(!is.null(tad.ranges)){
        Boundaries.obj <- Format_boundaries_normal_heatmap(Bricks = Bricks, 
            Ranges = tad.ranges, group.col = group.col, 
        cut.corners = cut.corners, colour.col = tad.colour.col, 
        colours = colours, colours.names = colours.names, 
        region.chr = x.coord.parsed['chr'], 
        region.start = as.numeric(x.coord.parsed['start']), 
        region.end = as.numeric(x.coord.parsed['end']), distance = distance,
        rotate = rotate)
    }

    if(rotate){
        ids <- xcoords <- ycoords <- NULL
        ThePlot <- ggplot(Entire.rotated.map, aes(x = xcoords, y = ycoords))
        ThePlot <- ThePlot + geom_polygon(aes(fill = values, group = ids))
        xlims <- c(0,max(Entire.rotated.map[,"xcoords"]))
        ylims <- c(min(Entire.rotated.map[,"ycoords"]),
            max(Entire.rotated.map[,"ycoords"]))
        y.coord.breaks <- seq(ceiling(min(Entire.rotated.map[,"ycoords"])),
            ceiling(max(Entire.rotated.map[,"ycoords"])),
            length.out = y.axis.num.breaks)
        y.axis.coord.labs <- y.coord.breaks*2
    }else{
        Matrix.df$row <- Matrix.df$row - 0.5
        Matrix.df$col <- Matrix.df$col - 0.5
        ThePlot <- ggplot(Matrix.df, aes(x = row, y = col))
        ThePlot <- ThePlot + geom_tile(aes(fill = rescale))
        xlims <- c(min(Matrix.df$row) - 0.5,max(Matrix.df$row) + 0.5)
        ylims <- c(min(Matrix.df$col) - 0.5,max(Matrix.df$col) + 0.5)
    }
    if(!is.null(tad.ranges)){
        line.group <- x <- y <- NULL
        ThePlot <- ThePlot + geom_line(data = Boundaries.obj, 
            aes(x = x, y = y, group = line.group, colour = colours), 
            size = line.width) 
        ThePlot <- ThePlot + scale_colour_manual(values = colours)
    }
    ThePlot <- ThePlot + scale_x_continuous(limits = xlims, expand = c(0,0),
        breaks = x.coord.breaks, labels = x.axis.coord.labs)
    ThePlot <- ThePlot + scale_y_continuous(limits = ylims, expand = c(0,0),
        breaks = y.coord.breaks, labels = y.axis.coord.labs)
    ThePlot <- ThePlot + scale_fill_gradientn(legend.title,values = Value.dist, 
        breaks = Colour.breaks, labels = Colour.labs, colors = Colours)
    ThePlot<-ThePlot+ Brick_theme
    # return(Entire.rotated.map)
    ThePlot<-ThePlot+labs(title = Labels['title'], x = Labels['x.axis'], 
        y = Labels['y.axis'])
    ggsave(filename = File, plot = ThePlot, width = width, height = height, 
        units = units)
    if(return.object){
        return(ThePlot)
    }else{
        return(TRUE)
    }
}