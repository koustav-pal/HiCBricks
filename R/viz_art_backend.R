._Figure_out_genomic_scale <- function(positions){
    axis.labels <- vapply(positions,function(x){
        if(x > 1000000){
            return(paste(round(x/1000000,1), "mb", sep = ""))
        }
        if(x > 1000){
            return(paste(round(x/1000,1), "kb", sep = ""))
        }
        if(x >= 1){
            return(paste(x/1, "bp", sep = ""))
        }
    }, "")
    return(axis.labels)
}
._Parse_genomic_coordinates <- function(a.list = NULL){
    Reference.object <- GenomicMatrix$new()
    parsed.string.list <- lapply(a.list,function(x){
        Split.string <- Split_genomic_coordinates(Coordinate=x)
        a.vector <- c(Split.string[[1]][1], Split.string[[1]][2], 
            Split.string[[1]][3])
        names(a.vector) <- Reference.object$NonStrandedColNames
        a.vector
    })
    return(parsed.string.list)
}

Get_one_or_two_brick_regions <- function(Bricks = NULL, resolution = NULL, 
    x_coords = NULL, y_coords = NULL, distance = NULL, value_cap = NULL, 
    FUN = NULL){
    Reference.object <- GenomicMatrix$new()
    if(length(Bricks) > 2){
        stop("Higher order polygon layouts have not been implemented yet! 
            So for now we can only do two matrices at a time.\n")
    }
    if(!is.null(value_cap)){
        if(value_cap > 1 | value_cap < 0){
            stop("value_cap must be a value between 0,1 ",
                "representing the quantiles.\n")
        }
    }
    # require(reshape2)
    Matrix.df.list <- list()
    for(i in seq_along(Bricks)){
        Brick <- Bricks[[i]]
        Matrix <- Brick_get_matrix_within_coords(Brick = Brick, 
            x_coords = x_coords, y_coords = y_coords, 
            resolution = resolution, force = TRUE, FUN = FUN)
        Region.position.x <- Brick_return_region_position(Brick = Brick,
            region = x_coords, resolution = resolution)
        Region.position.y <- Brick_return_region_position(Brick = Brick,
            region = y_coords, resolution = resolution)
        if(dim(Matrix)[1] != length(Region.position.x) | 
            dim(Matrix)[2] != length(Region.position.y)){
            stop("Matrix dimensions do not match the expected ",
                "dimensions of the matrix! Please check the ",
                "value transformation!\n")
        }
        Matrix[is.nan(Matrix) | is.infinite(Matrix) | is.na(Matrix)] <- 0
        rownames(Matrix) <- Region.position.x
        colnames(Matrix) <- Region.position.y
        Matrix.df <- melt(Matrix)
        colnames(Matrix.df) <- c("row","col","val")

        if(!is.null(value_cap)){
            capped.val <- quantile(Matrix.df$val,value_cap)
            Matrix.df$val[Matrix.df$val > capped.val] <- capped.val
        }
        Matrix.df$dist <- Matrix.df$col - Matrix.df$row
        Matrix.df$keep <- FALSE
        if(i == 1){
            Matrix.df$keep[Matrix.df$dist >= 0] <- TRUE
        }else{
            Matrix.df$keep[Matrix.df$dist <= 0] <- TRUE
        }
        Matrix.df.list[[i]] <- Matrix.df
    }

    Matrix.df <- do.call(rbind,Matrix.df.list)
    if(length(Bricks)==2){
        Matrix.df <- Matrix.df[Matrix.df$keep,]
        Matrix.df$val[Matrix.df$dist == 0] <- 0
    }
    if(!is.null(distance)){
        Matrix.df <- Matrix.df[Matrix.df$dist <= distance & 
        Matrix.df$dist >= -distance,]
    }
    return(Matrix.df)
}

Make_axis_labels = function(Brick = NULL, chr = NULL, resolution = NULL, 
    positions = NULL){
    Bintable <- Brick_get_bintable(Brick = Brick, chr = chr, 
        resolution = resolution)
    breaks <- end(Bintable[positions])
    breaks[1] <- start(Bintable[positions[1]])
    coord.labs <- ._Figure_out_genomic_scale(breaks)
    return(coord.labs)
}

Make_colours <- function(palette = NULL, extrapolate_on = NULL, direction = 1){
    # require(RColorBrewer)
    # require(viridis)
    viridis.cols <- list("plasma" = plasma, "inferno" = inferno, 
        "magma" = magma, "viridis" = viridis, "cividis" = cividis)
    brewer.names <- rownames(brewer.pal.info)
    brewer.pals <- brewer.pal.info$category 
    brewer.names <- brewer.names[brewer.pals %in% c("div","seq")]
    if(length(palette)!=1){
        stop("palette must be a string of length 1\n")
    }
    if(!(palette %in% brewer.names) & 
        !(palette %in% names(viridis.cols))){
        stop("palette must be a palette from RColorBrewer ",
            "diverging or sequential or Viridis.\n")
    }
    if(palette %in% rownames(brewer.pal.info)){
        Category <- brewer.pal.info$category[
        rownames(brewer.pal.info) == palette]
        Colours <- brewer.pal(name = palette, 
            n = brewer.pal.info$maxcolors[rownames(brewer.pal.info) == palette])
        if(direction == -1){
            Colours <- rev(Colours)
        }
    }else{
        viridis.col.breaks <- 12
        Do.viridis <- TRUE
        viridis.fun <- viridis.cols[[palette]]
        Colours <- viridis.fun(n = viridis.col.breaks, direction = direction)
    }
    if(!is.null(extrapolate_on)){
        if(extrapolate_on > 100){
            stop("I don't think you can actually differentiate ",
                "between more than 100 shades.")
        }
        Colours <- colorRampPalette(Colours)(extrapolate_on)
    }
    return(Colours)
}

RotatePixels <- function(shift = NULL, plot.order = NULL, vector = NULL, 
    distance = NULL, upper = TRUE){
    leftmost.x <- seq((0+shift),length.out=length(vector),by=1)
    leftmost.y <- rep((0+shift),length(vector))
    Order<-seq(plot.order,length.out=length(leftmost.x),by=4)

    uppermid.x <- seq((0.5+shift),length.out=length(vector),by=1)
    uppermid.y <- rep((0.5+shift),length.out=length(vector),by=1)
    Order<-c(Order,seq(min(Order)+1,length.out=length(uppermid.y),by=4))

    rightmost.x <- seq((1+shift),length.out=length(vector),by=1)   
    rightmost.y <- rep((0+shift),length(vector))
    Order<-c(Order,seq(min(Order)+2,length.out=length(uppermid.y),by=4))

    lowermid.x <- leftmost.x + 0.5
    lowermid.y <- leftmost.y - 0.5
    Order <- c(Order,seq(min(Order)+3,length.out=length(uppermid.y),by=4))

    IDs <- factor(paste("Pixel",rightmost.x,distance+1,sep="."),
        levels=paste("Pixel",rightmost.x,distance+1,sep="."))

    if(!upper){
        IDs <- factor(paste("Pixel",rightmost.x,(distance+1)*-1,sep="."),
            levels=paste("Pixel",rightmost.x,(distance+1)*-1,sep="."))
    }

    IDsRep <- rep(IDs,4)

    horizontal.coords <- data.frame(xcoords = c(leftmost.x, uppermid.x, 
        rightmost.x, lowermid.x),
        ycoords = c(leftmost.y, uppermid.y, rightmost.y, lowermid.y), 
        ids = IDs, Distance = distance, Order=Order)
    if(!upper){
        horizontal.coords[,"ycoords"] <- horizontal.coords[,"ycoords"] * -1
    }

    horizontal.values <- data.frame(ids=IDs,values=vector)
    horizontal.data <- merge(horizontal.values, horizontal.coords, by=c("ids"))
    horizontal.data <- horizontal.data[order(horizontal.data$Order),]
    return(horizontal.data)
}

RotateHeatmap = function(Matrix=NULL, value.var=NULL, upper = FALSE){
    rotated.distance <- max(Matrix$dist)
    Rotated.df.list <- lapply(seq_len(rotated.distance),function(x){
        Distance <- x - 1
        Vector <- Matrix[Matrix$dist==x,value.var]
        Shift <- 0.5 * Distance
        Dataframe <- RotatePixels(shift = Shift, plot.order = 1, 
            vector = Vector, distance = Distance, upper = upper)
    })
    Rotated.df <- do.call(rbind,Rotated.df.list)
    return(Rotated.df)
}

make_boundaries_for_heatmap <- function(Object = NULL, region.start = NULL, 
    region.end = NULL, distance = NULL, cut_corners = FALSE){
    if(is.null(distance)){
        distance <- region.end - region.start
    }
    Unique.groups <- unique(Object[,"groups"])
    Object_split <- split(Object, paste(Object$groups, Object$dom.names, 
        Object$colours, sep = ":"))
    Group.list <- lapply(Object_split,function(current_domain){
        current_domain <- unique(current_domain)
        Brick.x <- unique(current_domain$groups)
        colours <- current_domain$colours[
        current_domain[,"type"] == "start"]
        domain_name <- current_domain$dom.names[
        current_domain[,"type"] == "start"]
        Start <- current_domain[
        current_domain[,"type"] == "start", "position"]
        End <- current_domain[
        current_domain[,"type"] == "end", "position"]
        if(Start < region.start & End >= region.start){
            # message("Here")
            Start <- region.start
            if((End - Start) > distance){
                Start <- Start + distance
            }
            Coord.list <- list(x1 = c(Start,End),
                y1 = c(End,End))
            Groups <- rep(paste(domain_name,Brick.x,c(1,2),sep = "."), 
                each = 2)
        }else if(Start >= region.start & End <= region.end){
            # message("2nd Here")
            Coord.list <- list(x1 = c(Start - 1,Start - 1,End), 
                y1 = c(Start - 1,End,End))
            Groups <- rep(paste(domain_name,Brick.x,sep = "."), 
                each = 3)
            My.end <- End
            My.Start <- Start
            if((End - Start) > distance){
                My.end <- Start + distance
                My.Start <- End - distance
            }
            Coord.list <- list(x1 = c(Start - 1,Start - 1,
                My.Start - 1,End), 
                y1 = c(Start - 1,My.end,End,End))
            Groups <- rep(paste(domain_name,Brick.x,
                c(1,2),
                sep = "."), each = 2)
        }else if(Start <= region.end & End > region.end){
            # message("3rd Here")
            End <- region.end
            if((End - Start) > distance){
                End <- Start + distance
            }
            Coord.list <- list(x1 = c(Start - 1,Start - 1),
            y1 = c(Start - 1,End))
            Groups <- rep(paste(domain_name,Brick.x,sep = "."),2)
        }
        if(Brick.x == 2){
            Coord.list <- rev(Coord.list)
        }
        Line <- data.frame(x = Coord.list[[1]], 
            y = Coord.list[[2]], colours = colours,
            line.group = Groups, group = paste("Group",Brick.x,sep = "."))
        if(length(Unique.groups)==1){
            Coord.list <- rev(Coord.list)
            Line.2 <- data.frame(x = Coord.list[[1]],
            y = Coord.list[[2]], colours = colours,
            line.group = Groups, group = paste("Group",Brick.x+1,sep = "."))
            Line <- rbind(Line,Line.2)
        }
        Line
    })
    Group.df <- do.call(rbind,Group.list)
}

make_boundaries_for_rotated_heatmap <- function(Object = NULL, 
    region.start = NULL, region.end = NULL, distance = NULL, 
    cut_corners = FALSE){
    if(is.null(distance)){
        distance <- region.end - region.start
    }
    Shift.seed <- 0.5
    Span <- region.end - region.start
    Unique.groups <- unique(Object[,"groups"])
    Object_split <- split(Object, paste(Object$groups, Object$dom.names, 
        Object$colours, sep = ":"))
    Group.list <- lapply(Object_split,function(current_domain){
        current_domain <- unique(current_domain)
        domain_name <- unique(current_domain[,"dom.names"])
        Brick.x <- unique(current_domain[,"groups"])
        colours <- current_domain$colours[
        current_domain[,"type"] == "start"]
        Start <- current_domain[
        current_domain[,"type"] == "start", "position"]
        End <- current_domain[
        current_domain[,"type"] == "end", "position"]

        Normalised.start.bin <- Start - region.start
        Normalised.end.bin <- End - region.start

        if(cut_corners){
            Max.dist <- (End - Start)/2
            if(Max.dist > distance){
                Max.dist <- distance/2
            }
        }else{
            Max.dist <- distance/2
        }
        Dist.up <- Max.dist
        if((Normalised.end.bin - (Max.dist*2)) < 0){
            Dist.up <- abs(0 - Normalised.end.bin)/2
        }
        x1.start <- Normalised.end.bin - Dist.up
        y1.start <- ifelse(Brick.x == 1, Dist.up, Dist.up*-1)
        x2.start <- Normalised.end.bin
        y2.start <- 0
        End.line <- data.frame(x=c(x1.start,x2.start),
        y=c(y1.start,y2.start), colours = colours,
        line.group = paste(Brick.x, domain_name, "end", sep = "."), 
        group = paste("Group",Brick.x,sep = "."),
        row.names = NULL)

        Dist.down <- Max.dist
        if((Normalised.start.bin + (Max.dist*2)) > Span){
            Dist.down <- (Span - Normalised.start.bin)/ 2
        }
        x1.end <- Normalised.start.bin
        y1.end <- 0
        x2.end <- Normalised.start.bin + Dist.down
        y2.end <- ifelse(Brick.x == 1, Dist.down, Dist.down*-1)
        Start.line <- data.frame(x=c(x1.end,x2.end),
                y=c(y1.end,y2.end), colours=colours,
                line.group = paste(Brick.x, domain_name, 
                    "start", sep = "."),
                group=paste("Group",Brick.x,sep = "."),row.names = NULL)
        Lines <- rbind(End.line,Start.line)
    })
    Group.df <- do.call(rbind,Group.list)
    return(Group.df)
}

Format_boundaries_normal_heatmap <- function(Bricks = NULL, resolution, 
    Ranges = NULL, group_col = NULL, cut_corners = FALSE, colour.col = NULL, 
    colours = NULL, colours_names = NULL, region.chr = NULL, 
    region.start = NULL, region.end = NULL, distance = NULL, 
    rotate = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(!is.null(group_col)){
        Col.values <- unique(elementMetadata(Ranges)[[group_col]])
        if(length(Col.values) > 2 | !is.numeric(Col.values)){
            stop("group_col values must be numeric ",
                "values corresponding to ",
                "the number of Brick objects ",
                "(max. 2) specified.\n")
        }
    }else{
        group_col <- "pseudogroups"
        Ranges.too <- Ranges
        elementMetadata(Ranges)[[group_col]] <- 1
        elementMetadata(Ranges.too)[[group_col]] <- 2
        Ranges <- c(Ranges,Ranges.too)
    }
    
    if(is.null(colours)){
        stop("colours expects a vector of colours of at least length 1")
    }
    if(is.null(colour.col)){
        colour.col <- "pseudo.colour.col"
        elementMetadata(Ranges)[[colour.col]] <- "My_Group"
    }
    Unique.colour.cols <- unique(elementMetadata(Ranges)[[colour.col]])
    if(length(Unique.colour.cols) != length(colours)){
        stop("colours length must be equal to ",
            "number of unique values present in ",
            colour.col,"\n",
            "Length of colours: ",length(colours),"\n",
            "Length of colour names: ",length(Unique.colour.cols),"\n",
            "Names: ",paste(Unique.colour.cols,collapse=","))
    }
    if(is.null(colours_names)){
        colours_names <- Unique.colour.cols
        names(colours) <- colours_names
    }else{
        if(any(!(Unique.colour.cols %in% colours_names))){
            stop("Provided colours_names had differing values ",
                "from unique values present in colour.cols")
        }
        names(colours) <- colours_names
    }

    chr.ranges <- Ranges[seqnames(Ranges) == region.chr]
    chr.ranges <- chr.ranges[end(chr.ranges) >= region.start & 
    start(chr.ranges) <= region.end]
    region <- paste(region.chr, region.start, region.end, 
        sep = Reference.object$Ranges.separator)
    Region.positions <- Brick_return_region_position(Brick = Bricks[[1]], 
        resolution = resolution, 
        region = region)
    Range.to.df.list <- lapply(seq_along(Bricks),function(Brick.x){
        pos.ranges <- chr.ranges[
        elementMetadata(chr.ranges)[[group_col]] == Brick.x]
        chrs <- as.vector(seqnames(pos.ranges))
        start <- start(pos.ranges)
        end <- end(pos.ranges)
        A.ranges <- Brick_fetch_range_index(Brick = Bricks[[Brick.x]], 
            chr = chrs, start = start, end = end, resolution = resolution)
        Position.list <- A.ranges[seqnames(A.ranges) == region.chr]
        check_if_only_one_ranges <- function(x){
            all(!is.na(Position.list$Indexes[[x]]))
        }
        if(!all(vapply(seq_along(Position.list), 
            check_if_only_one_ranges, TRUE))){
            stop("All ranges did not overlap with the bintable.\n")
        }
        Range.positions.start <- vapply(seq_along(Position.list),
            function(x){(min(Position.list$Indexes[[x]]))},1)
        Range.positions.end <- vapply(seq_along(Position.list),
            function(x){(max(Position.list$Indexes[[x]]))},1)
        Range.positions.names <- vapply(seq_along(Position.list),
            function(x){names(Position.list[x])},"")
        Start.df <- data.frame(dom.names = Range.positions.names, 
            position = Range.positions.start, 
            groups = Brick.x, type = "start", 
            colours = elementMetadata(pos.ranges)[[colour.col]])
        End.df <- data.frame(dom.names = Range.positions.names, 
            position = Range.positions.end, 
            groups = Brick.x, type = "end", 
            colours = elementMetadata(pos.ranges)[[colour.col]])
        All.df <- rbind(Start.df,End.df)
        All.df$dom.names <- as.character(All.df$dom.names)
        All.df
    })
    Range.to.df <- do.call(rbind, Range.to.df.list)
    if(rotate){
        Normal.heatmap.lines <- make_boundaries_for_rotated_heatmap(
            Object = Range.to.df, cut_corners = cut_corners,
            region.start = min(Region.positions), 
            region.end = max(Region.positions), distance = distance)      
    }else{
        Normal.heatmap.lines <- make_boundaries_for_heatmap(
            Object = Range.to.df, region.start = min(Region.positions), 
            region.end = max(Region.positions), distance = distance)       
    }
    return(Normal.heatmap.lines)
}

Get_heatmap_theme <- function(x_axis=TRUE, y_axis=TRUE, 
    x_axis.text = NULL, y_axis.text = NULL, text_size = 10, 
    x_axis_text_size = 10, y_axis_text_size = 10,
    legend_title_text_size = 8, legend_text_size = 8, title_size = 10,
    legend_key_width = unit(3,"cm"), legend_key_height = unit(0.5,"cm")){
    if(!x_axis){
        x_axis.ticks <- element_blank()
        x_axis.text <- element_blank()
    }else{
        x_axis.ticks <-element_line(colour = "#000000")
        x_axis.text <- element_text(colour = "#000000", size = x_axis_text_size)
    }
    if(!y_axis){
        y_axis.ticks <- element_blank()
        y_axis.text <- element_blank()
    }else{
        y_axis.ticks <-element_line(colour = "#000000")
        y_axis.text <- element_text(colour = "#000000", size = y_axis_text_size)
    }
    Brick_theme <- theme_bw() + theme(text = element_text(size=text_size),
                plot.background=element_blank(),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                panel.background = element_blank(),
                axis.title.x=x_axis.text,
                axis.title.y=y_axis.text,
                axis.text.x = x_axis.text,
                axis.text.y = x_axis.text,
                axis.ticks.x = x_axis.ticks,
                axis.ticks.y = y_axis.ticks,
                legend.position="bottom",
                legend.key.height = legend_key_height,
                legend.key.width = legend_key_width,
                legend.title=element_text(size=legend_title_text_size),
                legend.text=element_text(size=legend_text_size),
                plot.title=element_text(size=title_size))
    return(Brick_theme)
}

Get_heatmap_titles <- function(title = NULL, x_axis_title = NULL, 
    y_axis_title = NULL, legend_title = NULL, x_coords = NULL, 
    y_coords = NULL, rotate = NULL){
    if(is.null(legend_title)){
        legend_title <- "Signal"
    }
    if(is.null(x_axis_title)){
        x_axis_title <- paste("Genomic position",x_coords,sep = " ")
    }
    if(is.null(y_axis_title)){
        y_axis_title <- paste("Genomic position",y_coords,sep = " ")
    }
    if(is.null(title)){
        title <- paste(x_coords,y_coords,sep = "-")
    }
    if(rotate){
        y_axis_title <- "distance in bins"
    }
    Labels <- c(x_axis_title, y_axis_title, title, legend_title)
    names(Labels) <- c("x_axis","y_axis","title","legend")
    return(Labels)
}

make_axis_coord_breaks <- function(from = NULL, to = NULL, 
    how.many = NULL, two.sample = FALSE){
    Breaks <- round(seq(from,to, length.out = how.many))
    if(two.sample){
        Breaks <- unique(c(rev(Breaks)*-1,Breaks))
    }
    return(Breaks)
}

rescale_values_for_colours <- function(Object = NULL, two.sample = FALSE){
    # require(scales)
    Object$rescale <- 0
    if(two.sample){
        Object$rescale[Object$dist >= 0] <- rescale(
            Object$val[Object$dist >= 0]*-1,c(0,0.5))
        Object$rescale[Object$dist <= 0] <- rescale(
            (Object$val[Object$dist <= 0]),c(0.5,1))
    }else{
        Object$rescale <- rescale(Object$val,c(0,1))
    }
    return(Object$rescale)
}

make_colour_breaks <- function(Object = NULL, how.many = NULL, 
    two.sample = NULL){
    values <- Object[,'rescale']
    distances <- Object[,'dist']
    if(two.sample){
        Value.dist.1 <- seq(min(values[distances >= 0]), 
            max(values[distances >= 0]), length.out = how.many)
        Value.dist.2 <- seq(min(values[distances <= 0]), 
            max(values[distances <= 0]), length.out = how.many)
        Value.dist <- unique(c(Value.dist.1,Value.dist.2))

    }else{
        Value.dist <- seq(min(values),max(values),length.out = how.many)
    }
    return(Value.dist)
}

get_legend_breaks <- function(Object = NULL, mid.val = 0.5, 
    how.many = 5, value_cap = NULL, colours = NULL, two.sample = NULL){
    Len <- length(values)
    values <- Object[,'rescale']
    distances <- Object[,'dist']
    original.values <- Object[,'val']
    Upper.tri <- distances >= 0
    Lower.tri <- distances <= 0
    if(two.sample){
        Colour.breaks.1 <- seq(min(values),mid.val,length.out = 3)
        Colour.labs.1 <- round(seq(max(original.values[Upper.tri]), 
            min(original.values[Upper.tri]), length.out = 3), 2)
        Colour.breaks.2 <- seq(mid.val,max(values),length.out = 3)
        Colour.labs.2 <- round(seq(min(original.values[Lower.tri]), 
            max(original.values[Lower.tri]),length.out = 3),2)
        Colour.labs <- round(
            c(
                Colour.labs.1, 
                Colour.labs.2[2:length(Colour.labs.2)]),2)
        Colour.breaks <- unique(
            c(Colour.breaks.1,Colour.breaks.2))
        Colours <- c(rev(colours),colours)
        if(!is.null(value_cap)){
            Colour.labs[1] <- paste(">",Colour.labs[1],sep = "")
            Colour.labs[length(Colour.labs)] <- paste(">",
                Colour.labs[length(Colour.labs)],sep = "")
        }
        # message(Colour.labs,"\n")
    }else{
        Colour.breaks <- seq(min(values),max(values),length.out = 5)
        Colour.labs <- format(seq(min(original.values),
            max(original.values),length.out = 5), scientific = TRUE)
        if(!is.null(value_cap)){
            Colour.labs[length(Colour.labs)] <- paste(">",
                Colour.labs[length(Colour.labs)], sep = "")
        }
        Colours <- colours
    }
    return(list("cols" = Colours, 
        "col.breaks" = Colour.breaks, "col.labs" = Colour.labs))
}
