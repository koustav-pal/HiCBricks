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

Get_one_or_two_lego_regions <- function(Legos = NULL, x.coords = NULL, 
    y.coords = NULL, distance = NULL, 
    value.cap = NULL, FUN = NULL){
    Reference.object <- GenomicMatrix$new()
    if(length(Legos) > 2){
        stop("Polygonal layouts have not been implemented yet! 
            So for now we can only do two matrices at a time.\n")
    }
    if(!is.null(value.cap)){
        if(value.cap > 1 | value.cap < 0){
            stop("value.cap must be a value between 0,1 ",
                "representing the quantiles.\n")
        }
    }
    # require(reshape2)
    Matrix.df.list <- list()
    for(i in seq_along(Legos)){
        Lego <- Legos[i]
        Matrix <- Lego_get_matrix_within_coords(Lego = Lego, 
            x.coords = x.coords, y.coords = y.coords, force = TRUE, FUN = FUN)
        Region.position.x <- Lego_return_region_position(Lego = Lego, 
            region = x.coords)
        Region.position.y <- Lego_return_region_position(Lego = Lego, 
            region = y.coords)
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

        if(!is.null(value.cap)){
            capped.val <- quantile(Matrix.df$val,value.cap)
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
    if(length(Legos)==2){
        Matrix.df <- Matrix.df[Matrix.df$keep,]
        Matrix.df$val[Matrix.df$dist == 0] <- 0
    }
    if(!is.null(distance)){
        Matrix.df <- Matrix.df[Matrix.df$dist <= distance & 
        Matrix.df$dist >= -distance,]
    }
    return(Matrix.df)
}

Make_axis_labels = function(Lego = NULL, chr = NULL, positions = NULL){
    Bintable <- Lego_get_bintable(Lego = Lego, chr = chr)
    breaks <- end(Bintable[positions])
    breaks[1] <- start(Bintable[positions[1]])
    coord.labs <- ._Figure_out_genomic_scale(breaks)
    return(coord.labs)
}

Make_colours <- function(palette = NULL, extrapolate.on = NULL, direction = 1){
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
    if(!is.null(extrapolate.on)){
        if(extrapolate.on > 100){
            stop("I don't think you can actually differentiate ",
                "between more than 100 shades.")
        }
        Colours <- colorRampPalette(Colours)(extrapolate.on)
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
    region.end = NULL, distance = NULL, cut.corners = FALSE){
    if(is.null(distance)){
        distance <- region.end - region.start
    }
    cat(region.start,region.end,"\n")

    Unique.groups <- unique(Object[,"groups"])
    Group.list <- lapply(Unique.groups,function(Lego.x){
        cat(Lego.x,"\n")
        Domain <- Object[Object[,"groups"] == Lego.x,]
        Domain.names <- unique(Domain[,"dom.names"])
        Dolly.the.sheep.list <- lapply(Domain.names, function(domain.name){
            current.domain <- Domain[Domain[,"dom.names"]==domain.name,]
            colours <- current.domain$colours[
            current.domain[,"type"] == "start"]
            Start <- current.domain[
            current.domain[,"type"] == "start", "position"]
            End <- current.domain[
            current.domain[,"type"] == "end", "position"]
            if(Start < region.start & End >= region.start){
                Start <- region.start
                if((End - Start) > distance){
                    Start <- Start + distance
                }
                Coord.list <- list(x1 = c(Start,End),
                    y1 = c(End,End))
                Groups <- rep(paste(domain.name,Lego.x,c(1,2),sep = "."), 
                    each = 2)
            }else if(Start >= region.start & End <= region.end){
                Coord.list <- list(x1 = c(Start,Start,End), 
                    y1 = c(Start,End,End))
                Groups <- rep(paste(domain.name,Lego.x,sep = "."), 
                    each = 3)
                My.end <- End
                My.Start <- Start
                if((End - Start) > distance){
                    My.end <- Start + distance
                    My.Start <- End - distance
                }
                Coord.list <- list(x1 = c(Start,Start,My.Start,End), 
                    y1 = c(Start,My.end,End,End))
                Groups <- rep(paste(domain.name,Lego.x,
                    c(1,2),
                    sep = "."), each = 2)
            }else if(Start < region.end & End > region.end){
                End <- region.end
                if((End - Start) > distance){
                    End <- Start + distance
                }
                Coord.list <- list(x1 = c(Start,Start), y1 = c(Start,End))
                Groups <- rep(paste(domain.name,Lego.x,sep = "."),2)
            }
            if(Lego.x == 2){
                Coord.list <- rev(Coord.list)
            }
            Line <- data.frame(x = Coord.list[[1]] - 0.5, 
                y = Coord.list[[2]] + 0.5, colours = colours,
                line.group = Groups, group = paste("Group",Lego.x,sep = "."))
            if(length(Unique.groups)==1){
                Coord.list <- rev(Coord.list)
                Line.2 <- data.frame(x = Coord.list[[1]] - 0.5,
                y = Coord.list[[2]] + 0.5, colours = colours,
                line.group = Groups, group = paste("Group",Lego.x+1,sep = "."))
                Line <- rbind(Line,Line.2)
            }
            Line
        })
        Dolly.the.sheep <- do.call(rbind,Dolly.the.sheep.list)
        Dolly.the.sheep
    })
    Group.df <- do.call(rbind,Group.list)
}

make_boundaries_for_rotated_heatmap <- function(Object = NULL, 
    region.start = NULL, region.end = NULL, distance = NULL, 
    cut.corners = FALSE){
    if(is.null(distance)){
        distance <- region.end - region.start
    }
    Shift.seed <- 0.5
    Span <- region.end - region.start
    Unique.groups <- unique(Object[,"groups"])
    Group.list <- lapply(Unique.groups,function(Lego.x){
        Domain <- Object[Object[,"groups"] == Lego.x,]
        Domain.names <- unique(Domain[,"dom.names"])
        Domain.df.list <- lapply(Domain.names,function(x){
            current.domain <- Domain[Domain[,"dom.names"]==x,]
            colours <- current.domain$colours[
            current.domain[,"type"] == "start"]
            Start <- current.domain[
            current.domain[,"type"] == "start", "position"]
            End <- current.domain[
            current.domain[,"type"] == "end", "position"]

            Normalised.start.bin <- Start - region.start
            Normalised.end.bin <- End - region.start

            if(cut.corners){
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
            y1.start <- Dist.up
            x2.start <- Normalised.end.bin
            y2.start <- 0
            End.line <- data.frame(x=c(x1.start,x2.start),
            y=c(y1.start,y2.start), colours = colours,
            line.group = paste(x, "end", sep = "."), 
            group = paste("Group",Lego.x,sep = "."),
            row.names = NULL)

            Dist.down <- Max.dist
            if((Normalised.start.bin + (Max.dist*2)) > Span){
                Dist.down <- (Span - Normalised.start.bin)/ 2
            }
            x1.end <- Normalised.start.bin
            y1.end <- 0
            x2.end <- Normalised.start.bin + Dist.down
            y2.end <- Dist.down
            Start.line <- data.frame(x=c(x1.end,x2.end),
                    y=c(y1.end,y2.end), colours=colours,
                    line.group = paste(x, "start", sep = "."),
                    group=paste("Group",Lego.x,sep = "."),row.names = NULL)
            Lines <- rbind(End.line,Start.line)
        })
        Domain.df <- do.call(rbind,Domain.df.list)
    })
    Group.df <- do.call(rbind,Group.list)
    return(Group.df)
}

Format_boundaries_normal_heatmap <- function(Legos = NULL, Ranges = NULL, 
    group.col = NULL, cut.corners = FALSE, colour.col = NULL, 
    colours = NULL, colours.names = NULL, region.chr = NULL, 
    region.start = NULL, region.end = NULL, distance = NULL, 
    rotate = FALSE){
    Reference.object <- GenomicMatrix$new()
    if(!is.null(group.col)){
        Col.values <- unique(elementMetadata(Ranges)[[group.col]])
        if(!(length(Col.values) > 2 | !is.numeric(Col.values))){
            stop("group.col values must be numeric ",
                "values of for the two Lego objects.\n")
        }
    }else{
        group.col <- "pseudogroups"
        Ranges.too <- Ranges
        elementMetadata(Ranges)[[group.col]] <- 1
        elementMetadata(Ranges.too)[[group.col]] <- 2
        Ranges <- c(Ranges,Ranges.too)
    }
    
    Colour.check <- as.logical(as.numeric(is.null(colour.col)) * 
        as.numeric(is.null(colours)))
    Colour.check.2 <- length(colour.col)/length(colours)
    if(Colour.check){
        stop("colours expects a vector of length 1 with the colour value")
    }else if(is.nan(Colour.check.2)){
        stop("colours is missing\n")
    }else if(!(Colour.check.2 %in% c(0,1))){
        stop("colours and colour.col have different lengths\n")
    }else {
        colour.col <- "pseudo.colour.col"
        elementMetadata(Ranges)[[colour.col]] <- "My_Group"
    }
    if(is.null(colours.names)){
        colours.names <- unique(elementMetadata(Ranges)[[colour.col]])
        names(colours) <- colours.names
    }else{
        names(colours) <- colours.names
    }

    chr.ranges <- Ranges[seqnames(Ranges) %in% region.chr]
    chr.ranges <- chr.ranges[end(chr.ranges) >= region.start & 
    start(chr.ranges) <= region.end]
    region <- paste(region.chr, region.start, region.end, 
        sep = Reference.object$Ranges.separator)
    Region.positions <- Lego_return_region_position(Lego = Legos[1], 
        region = region)
    Range.to.df.list <- lapply(seq_along(Legos),function(Lego.x){
        pos.ranges <- chr.ranges[
        elementMetadata(chr.ranges)[[group.col]] == Lego.x]
        chrs <- as.vector(seqnames(pos.ranges))
        start <- start(pos.ranges)
        end <- end(pos.ranges)
        A.ranges <- Lego_fetch_range_index(Lego = Legos[Lego.x], chr = chrs, 
            start = start, end = end)
        Position.list <- A.ranges[seqnames(A.ranges) %in% region.chr]
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
            groups = Lego.x, type = "start", 
            colours = elementMetadata(pos.ranges)[[colour.col]])
        End.df <- data.frame(dom.names = Range.positions.names, 
            position = Range.positions.end, 
            groups = Lego.x, type = "end", 
            colours = elementMetadata(pos.ranges)[[colour.col]])
        All.df <- rbind(Start.df,End.df)
        All.df$dom.names <- as.character(All.df$dom.names)
        All.df
    })
    Range.to.df <- do.call(rbind, Range.to.df.list)
    if(rotate){
        Normal.heatmap.lines <- make_boundaries_for_rotated_heatmap(
            Object = Range.to.df, cut.corners = cut.corners,
            region.start = min(Region.positions), 
            region.end = max(Region.positions), distance = distance)      
    }else{
        Normal.heatmap.lines <- make_boundaries_for_heatmap(
            Object = Range.to.df, region.start = min(Region.positions), 
            region.end = max(Region.positions), distance = distance)       
    }
    return(Normal.heatmap.lines)
}

Get_heatmap_theme <- function(x.axis=TRUE, y.axis=TRUE, 
    x.axis.text = NULL, y.axis.text = NULL, text.size = 10, 
    x.axis.text.size = 10, y.axis.text.size = 10,
    legend.title.text.size = 8, legend.text.size = 8, title.size = 10,
    legend.key.width = unit(3,"cm"), legend.key.height = unit(0.5,"cm")){
    if(!x.axis){
        x.axis.ticks <- element_blank()
        x.axis.text <- element_blank()
    }else{
        x.axis.ticks <-element_line(colour = "#000000")
        x.axis.text <- element_text(colour = "#000000", size = x.axis.text.size)
    }
    if(!y.axis){
        y.axis.ticks <- element_blank()
        y.axis.text <- element_blank()
    }else{
        y.axis.ticks <-element_line(colour = "#000000")
        y.axis.text <- element_text(colour = "#000000", size = y.axis.text.size)
    }
    Lego_theme <- theme_bw() + theme(text = element_text(size=text.size),
                plot.background=element_blank(),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                panel.background = element_blank(),
                axis.title.x=x.axis.text,
                axis.title.y=y.axis.text,
                axis.text.x = x.axis.text,
                axis.text.y = x.axis.text,
                axis.ticks.x = x.axis.ticks,
                axis.ticks.y = y.axis.ticks,
                legend.position="bottom",
                legend.key.height = legend.key.height,
                legend.key.width = legend.key.width,
                legend.title=element_text(size=legend.title.text.size),
                legend.text=element_text(size=legend.text.size),
                plot.title=element_text(size=title.size))
    return(Lego_theme)
}

Get_heatmap_titles <- function(title = NULL, x.axis.title = NULL, 
    y.axis.title = NULL, legend.title = NULL, x.coords = NULL, 
    y.coords = NULL, rotate = NULL){
    if(is.null(legend.title)){
        legend.title <- "Signal"
    }
    if(is.null(x.axis.title)){
        x.axis.title <- paste("Genomic position",x.coords,sep = " ")
    }
    if(is.null(y.axis.title)){
        y.axis.title <- paste("Genomic position",y.coords,sep = " ")
    }
    if(is.null(title)){
        title <- paste(x.coords,y.coords,sep = "-")
    }
    if(rotate){
        y.axis.title <- "distance in bins"
    }
    Labels <- c(x.axis.title, y.axis.title, title, legend.title)
    names(Labels) <- c("x.axis","y.axis","title","legend")
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
    how.many = 5, value.cap = NULL, colours = NULL, two.sample = NULL){
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
        if(!is.null(value.cap)){
            Colour.labs[1] <- paste(">",Colour.labs[1],sep = "")
            Colour.labs[length(Colour.labs)] <- paste(">",
                Colour.labs[length(Colour.labs)],sep = "")
        }
        cat(Colour.labs,"\n")
    }else{
        Colour.breaks <- seq(min(values),max(values),length.out = 5)
        Colour.labs <- round(seq(min(original.values),
            max(original.values),length.out = 5),2)
        if(!is.null(value.cap)){
            Colour.labs[length(Colour.labs)] <- paste(">",
                Colour.labs[length(Colour.labs)], sep = "")
        }
        Colours <- colours
    }
    return(list("cols" = Colours, 
        "col.breaks" = Colour.breaks, "col.labs" = Colour.labs))
}