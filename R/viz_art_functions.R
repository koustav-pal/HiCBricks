._Figure_out_genomic_scale <- function(positions){
	axis.labels <- vapply(positions,function(x){
		if(x > 1000000){
			return(paste(round(x/1000000,1), "mb", sep = ""))
		}
		if(x > 1000){
			return(paste(round(x/1000,1), "kb", sep = ""))
		}
		if(x > 1){
			return(paste(x/1, "bp", sep = ""))
		}
	}, "")
	return(axis.labels)
}
._Parse_genomic_coordinates <- function(a.list = NULL){
	Reference.object <- GenomicMatrix$new()
	parsed.string.list <- lapply(a.list,function(x){
		Split.string <- Split_genomic_coordinates(Coordinate=x)
		a.vector <- c(Split.string[[1]][1], Split.string[[1]][2], Split.string[[1]][3])
		names(a.vector) <- Reference.object$NonStrandedColNames
		a.vector
	})
	return(parsed.string.list)
}

Get_one_or_two_lego_regions <- function(Legos = NULL, x.coords = NULL, y.coords = NULL, distance = NULL, 
	value.cap = NULL, FUN = NULL){
	Reference.object <- GenomicMatrix$new()
	if(length(Legos) > 2){
		stop("Polygonal layouts have not been implemented yet! 
			So for now we can only do two matrices at a time.\n")
	}
	if(!is.null(value.cap)){
		if(value.cap > 1 | value.cap < 0){
			stop("value.cap must be a value between 0,1 representing the quantiles.\n")
		}
	}
	# require(reshape2)
	Matrix.df.list <- list()
	for(i in seq_along(Legos)){
		Lego <- Legos[i]
		Matrix <- Lego_get_matrix_within_coords(Lego = Lego, x.coords = x.coords, y.coords = y.coords, FUN = FUN)
		Region.position.x <- Lego_return_region_position(Lego = Lego, region = x.coords)
		Region.position.y <- Lego_return_region_position(Lego = Lego, region = y.coords)
		if(dim(Matrix)[1] != length(Region.position.x) | dim(Matrix)[2] != length(Region.position.y)){
			stop("Matrix dimensions do not match the expected dimensions of the matrix!
				Please check the value transformation!\n")
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
		Matrix.df <- Matrix.df[Matrix.df$dist <= distance & Matrix.df$dist >= -distance,]
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
    viridis.cols <- list("plasma" = plasma, "inferno" = inferno, "magma" = magma, "viridis" = viridis, "cividis" = cividis)
    brewer.names <- rownames(brewer.pal.info)
    brewer.pals <- brewer.pal.info$category 
    brewer.names <- brewer.names[brewer.pals %in% c("div","seq")]
    if(length(palette)!=1){
    	stop("palette must be a string of length 1\n")
    }
    if(!(palette %in% brewer.names) & !(palette %in% names(viridis.cols))){
    	stop("palette must be a palette from RColorBrewer diverging or sequential or Viridis.\n")
    }
    if(palette %in% rownames(brewer.pal.info)){
    	Category <- brewer.pal.info$category[rownames(brewer.pal.info) == palette]
    	Colours <- brewer.pal(name = palette, n = brewer.pal.info$maxcolors[rownames(brewer.pal.info) == palette])
    }else{
    	viridis.col.breaks <- 12
    	Do.viridis <- TRUE
    	viridis.fun <- viridis.cols[[palette]]
    	Colours <- viridis.fun(n = viridis.col.breaks, direction = direction)
    }
	if(!is.null(extrapolate.on)){
		if(extrapolate.on > 100){
			stop("I don't think you can actually differentiate between more than 100 shades.")
		}
		Colours <- colorRampPalette(Colours)(extrapolate.on)
	}
	return(Colours)
}

RotatePixels <- function(shift = NULL, plot.order = NULL, vector = NULL, distance = NULL, upper = TRUE){
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

    IDs <- factor(paste("Pixel",rightmost.x,distance+1,sep="."),levels=paste("Pixel",rightmost.x,distance+1,sep="."))

    if(!upper){
	    IDs <- factor(paste("Pixel",rightmost.x,(distance+1)*-1,sep="."),levels=paste("Pixel",rightmost.x,(distance+1)*-1,sep="."))
    }

    IDsRep <- rep(IDs,4)

    horizontal.coords <- data.frame(xcoords = c(leftmost.x, uppermid.x, rightmost.x, lowermid.x),
        ycoords = c(leftmost.y, uppermid.y, rightmost.y, lowermid.y), ids = IDs, Distance = distance, Order=Order)
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
	Rotated.df.list <- lapply(1:rotated.distance,function(x){
	    Distance <- x - 1
	    Vector <- Matrix[Matrix$dist==x,value.var]
	    Shift <- 0.5 * Distance
	    Dataframe <- RotatePixels(shift = Shift, plot.order = 1, vector = Vector, distance = Distance, upper = upper)
	})
	Rotated.df <- do.call(rbind,Rotated.df.list)
	return(Rotated.df)
}

make_boundaries_for_heatmap <- function(Object = NULL, region.start = NULL, 
	region.end = NULL, distance = NULL, cut.corners = FALSE){
	if(is.null(distance)){
		distance <- region.end - region.start
	}
	Shift.seed <- 0.5
	Unique.groups <- unique(Object[,"groups"])
	Group.list <- lapply(Unique.groups,function(Lego.x){
		Domain <- Object[Object[,"groups"] == Lego.x,]
		Domain.names <- unique(Domain[,"dom.names"])
		Dolly.the.sheep.list <- lapply(Domain.names, function(domain.name){
			current.domain <- Domain[Domain[,"dom.names"]==domain.name,]
			colours <- current.domain$colours[current.domain[,"type"] == "start"]
		    Start <- current.domain[current.domain[,"type"] == "start", "position"]
		    End <- current.domain[current.domain[,"type"] == "end", "position"]
		    if(Start < region.start & End >= region.start){
			    Start <- region.start
			    if((End - Start) > distance){
			    	Start <- Start + distance
			    }
		    	Coord.list <- list(x1 = c(Start,End),y1 = c(End,End))
		    	Groups <- rep(paste(domain.name,Lego.x,c(1,2),sep = "."), each = 2)
		    }else if(Start >= region.start & End <= region.end){
			    Coord.list <- list(x1 = c(Start,Start,End), y1 = c(Start,End,End))
			    Groups <- rep(paste(domain.name,Lego.x,sep = "."), each = 3)
			    My.end <- End
			    My.Start <- Start
			    if((End - Start) > distance){
			    	My.end <- Start + distance
			    	My.Start <- End - distance
			    }
		    	Coord.list <- list(x1 = c(Start,Start,My.Start,End), y1 = c(Start,My.end,End,End))
		    	Groups <- rep(paste(domain.name,Lego.x,c(1,2),sep = "."), each = 2)
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
		    Line <- data.frame(x = Coord.list[[1]] - 0.5, y = Coord.list[[2]] + 0.5, colours = colours,
		    	line.group = Groups, group = paste("Group",Lego.x,sep = "."))
		    Line
		})
		Dolly.the.sheep <- do.call(rbind,Dolly.the.sheep.list)
		Dolly.the.sheep
	})
	Group.df <- do.call(rbind,Group.list)
}

make_boundaries_for_rotated_heatmap <- function(Object = NULL, region.start = NULL, 
	region.end = NULL, distance = NULL, cut.corners = FALSE){
	if(is.null(distance)){
		distance <- region.end - region.start
	}
	Shift.seed <- 0.5
	Unique.groups <- unique(Object[,"groups"])
	Group.list <- lapply(Unique.groups,function(Lego.x){
		Domain <- Object[Object[,"groups"] == Lego.x,]
		Domain.names <- unique(Domain[,"dom.names"])
			Domain.df.list <- lapply(Domain.names,function(x){
				current.domain <- Domain[Domain[,"dom.names"]==domain.name,]
				colours <- current.domain$colours[current.domain[,"type"] == "start"]
			    Start <- current.domain[current.domain[,"type"] == "start", "position"]
			    End <- current.domain[current.domain[,"type"] == "end", "position"]
			    Mid.bin <- Start - region.start
			    Max.dist <- distance/2
			    if(Mid.bin - (Max.dist*2) < 0){
			    	Dist.up <- abs(0 - Mid.bin)/2
			    }

			    Lines
			})
			Domain.df <- do.call(rbind,Domain.df.list)
		Dolly.the.sheep <- do.call(rbind,Dolly.the.sheep.list)
		Dolly.the.sheep
	})
	Group.df <- do.call(rbind,Group.list)
}


Format_boundaries_normal_heatmap <- function(Legos = NULL, Ranges = NULL, group.col = NULL, 
	cut.corners = NULL, colour.col = NULL, colours = NULL, colours.names = NULL, region.chr = NULL, 
	region.start = NULL, region.end = NULL, distance = NULL, rotate = FALSE){
	Reference.object <- GenomicMatrix$new()
	if(!is.null(group.col)){
		Col.values <- unique(mcols(Ranges,group.col,use.names = TRUE)[[1]])
		if(!(length(Col.values) > 2 | class(Col.values)!="numeric")){
			stop("group.col values must be numeric values of for the two Lego objects.\n")
		}
	}else{
		group.col <- "pseudogroups"
		Ranges.too <- Ranges
		mcols(Ranges, group.col, use.names = TRUE) <- 1
		mcols(Ranges.too, group.col, use.names = TRUE) <- 2
		Ranges <- c(Ranges,Ranges.too)
	}
	
	Colour.check <- as.logical(as.numeric(is.null(colour.col)) * as.numeric(is.null(colours)))
	Colour.check.2 <- length(colour.col)/length(colours)
	if(Colour.check){
		stop("colours expects a vector of length 1 with the colour value")
	}else if(is.nan(Colour.check.2)){
		stop("colours is missing\n")
	}else if(!(Colour.check.2 %in% c(0,1))){
		stop("colours and colour.col have different lengths\n")
	}else {
		colour.col <- "pseudo.colour.col"
		mcols(Ranges, colour.col, use.names = TRUE) <- "My_Group"
	}
	if(is.null(colours.names)){
		colours.names <- unique(mcols(Ranges, colour.col, use.names = TRUE))
		names(colours) <- colours.names
	}else{
		names(colours) <- colours.names
	}

	chr.ranges <- Ranges[seqnames(Ranges) %in% region.chr]
	chr.ranges <- chr.ranges[start(chr.ranges) < region.end & end(chr.ranges) > region.start]
	region <- paste(region.chr, region.start, region.end, sep = Reference.object$Ranges.separator)
	Region.positions <- Lego_return_region_position(Lego = Legos[1], region = region)

	Range.to.df.list <- lapply(seq_along(Legos),function(Lego.x){
		pos.ranges <- chr.ranges[mcols(chr.ranges,group.col,use.names = TRUE) == Lego.x]
		chrs <- as.vector(seqnames(pos.ranges))
		start <- start(pos.ranges)
		end <- end(pos.ranges)
		A.list <- Lego_fetch_range_index(Lego = Legos[Lego.x], chr = chrs, start = start, end = end)
		Position.list <- A.list[[region.chr]]
		check_if_only_one_ranges <- function(x){
			!is.na(x[["Indexes"]])
		}
		if(!all(vapply(Position.list, check_if_only_one_ranges, TRUE))){
			stop("All ranges did not overlap with the bintable.\n")
		}
		Range.positions.start <- vapply(Position.list,function(x){min(x[["Indexes"]])},1)
		Range.positions.end <- vapply(Position.list,function(x){max(x[["Indexes"]])},1)
		Range.positions.names <- vapply(Position.list,function(x){names(x[["SubjectInfo"]])},"")
		Start.df <- data.frame(dom.names = Range.positions.names, position = Range.positions.start, 
			groups = Lego.x, type = "start",
			colours = mcols(Ranges, colour.col, use.names = TRUE))
		End.df <- data.frame(dom.names = Range.positions.names, position = Range.positions.end, 
			groups = Lego.x, type = "end",
			colours = mcols(Ranges, colour.col, use.names = TRUE))
		All.df <- rbind(Start.df,End.df)
		All.df$dom.names <- as.character(All.df$dom.names)
	})
	Range.to.df <- do.call(rbind, Range.to.df.list)
	if(rotate){
		Normal.heatmap.lines <- make_boundaries_for_rotated_heatmap(Object = Range.to.df, 
			region.start = region.start, region.end = region.end, distance = distance, rotate = rotate)		
	}else{
		Normal.heatmap.lines <- make_boundaries_for_heatmap(Object = Range.to.df, region.start = region.start, 
			region.end = region.end, distance = distance)		
	}
	Plot.object <- geom_line(data = Normal.heatmap.lines, aes(x = x, y = y, group = line.group, colour = colours))	
	Plot.object <- Plot.object + scale_colour_manual(values = colours)
	return(Plot.object)
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

Get_heatmap_titles <- function(title = NULL, x.axis.title = NULL, y.axis.title = NULL, 
	legend.title = NULL, x.coords = NULL, y.coords = NULL, rotate = NULL){
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

make_axis_coord_breaks <- function(from = NULL, to = NULL, how.many = NULL, two.sample = FALSE){
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
    	Object$rescale[Object$dist >= 0] <- rescale(Object$val[Object$dist >= 0]*-1,c(0,0.5))
    	Object$rescale[Object$dist <= 0] <- rescale((Object$val[Object$dist <= 0]),c(0.5,1))
	}else{
		Object$rescale <- rescale(Object$val,c(0,1))
	}
	return(Object$rescale)
}

make_colour_breaks <- function(Object = NULL, how.many = NULL, two.sample = NULL){
	values <- Object[,'rescale']
	distances <- Object[,'dist']
	if(two.sample){
		Value.dist.1 <- seq(min(values[distances >= 0]), max(values[distances >= 0]), length.out = how.many)
		Value.dist.2 <- seq(min(values[distances <= 0]), max(values[distances <= 0]), length.out = how.many)
		Value.dist <- unique(c(Value.dist.1,Value.dist.2))

	}else{
		Value.dist <- seq(min(values),max(values),length.out = how.many)
	}
	return(Value.dist)
}

get_legend_breaks <- function(Object = NULL, mid.val = 0.5, how.many = 5, value.cap = NULL, colours = NULL, two.sample = NULL){
	Len <- length(values)
	values <- Object[,'rescale']
	distances <- Object[,'dist']
	original.values <- Object[,'val']
	Upper.tri <- distances >= 0
	Lower.tri <- distances <= 0
	if(two.sample){
		Colour.breaks.1 <- seq(min(values),mid.val,length.out = 3)
		Colour.labs.1 <- round(seq(max(original.values[Upper.tri]), min(original.values[Upper.tri]), length.out = 3), 2)
		Colour.breaks.2 <- seq(mid.val,max(values),length.out = 3)
		Colour.labs.2 <- round(seq(min(original.values[Lower.tri]),max(original.values[Lower.tri]),length.out = 3),2)
		Colour.labs <- round(c(Colour.labs.1,Colour.labs.2[2:length(Colour.labs.2)]),2)
		Colour.breaks <- unique(c(Colour.breaks.1,Colour.breaks.2))
		Colours <- c(rev(colours),colours)
    	if(!is.null(value.cap)){
	    	Colour.labs[1] <- paste(">",Colour.labs[1],sep = "")
	    	Colour.labs[length(Colour.labs)] <- paste(">",Colour.labs[length(Colour.labs)],sep = "")
    	}
    	cat(Colour.labs,"\n")
	}else{
    	Colour.breaks <- seq(min(values),max(values),length.out = 5)
    	Colour.labs <- round(seq(min(original.values),max(original.values),length.out = 5),2)
    	if(!is.null(value.cap)){
	    	Colour.labs[length(Colour.labs)] <- paste(">",Colour.labs[length(Colour.labs)], sep = "")
    	}
    	Colours <- colours
	}
	return(list("cols" = Colours, "col.breaks" = Colour.breaks, "col.labs" = Colour.labs))
}

Lego_vizart_plot_heatmap = function(File = NULL, Legos = NULL, x.coords = NULL, y.coords = NULL, 
	FUN = NULL, value.cap = NULL, distance = NULL, rotate = FALSE, x.axis = TRUE, x.axis.title = NULL,
	y.axis = TRUE, y.axis.title = NULL, title = NULL, legend.title = NULL, return.object=FALSE,
	x.axis.num.breaks = 5, y.axis.num.breaks = 5, palette = NULL, col.direction = 1, extrapolate.on = NULL, 
	x.axis.text.size = 10, y.axis.text.size = 10, text.size = 10, legend.title.text.size = 8, legend.text.size = 8, 
	title.size = 10, tad.ranges = NULL, group.col = NULL, tad.colour.col = NULL, colours = NULL, colours.names = NULL,
	cut.corners = NULL, highlight.points = NULL, width = 10, height = 6, units = "cm", legend.key.width = unit(3,"cm"), 
	legend.key.height = unit(0.5,"cm")){

	Matrix.df <- Get_one_or_two_lego_regions(Legos = Legos, x.coords = x.coords, y.coords = y.coords, 
		distance = distance, value.cap = value.cap, FUN = FUN)
	if(nrow(Matrix.df)==0){
		stop("The matrix was empty!")
	}
	list_of_coords <- list("x.coords" = x.coords, "y.coords" = y.coords)

	Parsed_string <- ._Parse_genomic_coordinates(list_of_coords)
    x.coord.parsed <- Parsed_string[["x.coords"]]
    y.coord.parsed <- Parsed_string[["y.coords"]]
    x.coord.breaks <- make_axis_coord_breaks(from = min(Matrix.df$row), to = max(Matrix.df$row), 
    	how.many = x.axis.num.breaks, two.sample = FALSE)
    x.axis.coord.labs <- Make_axis_labels(Lego = Legos[1], chr = x.coord.parsed['chr'], positions = x.coord.breaks)

    two.sample <- (rotate & length(Legos)==2)
    y.coord.breaks <- make_axis_coord_breaks(from = min(Matrix.df$row), to = max(Matrix.df$row), 
    	how.many = x.axis.num.breaks, two.sample = two.sample)
    y.axis.coord.labs <- Make_axis_labels(Lego = Legos[1], chr = y.coord.parsed['chr'], positions = abs(y.coord.breaks))

    Colours <- Make_colours(palette = palette, extrapolate.on = extrapolate.on, direction = col.direction)
    # go from min val to mid to max val
    two.sample <- (length(Legos)==2)
    Matrix.df$rescale <- rescale_values_for_colours(Object = Matrix.df, two.sample = two.sample)
    Value.dist <- make_colour_breaks(Object = Matrix.df, how.many = length(Colours), two.sample = two.sample)
    
    Legend.breaks.list <- get_legend_breaks(Object = Matrix.df, how.many = 5, value.cap = value.cap, colours = Colours,
    	two.sample = two.sample)
    Colour.breaks <- Legend.breaks.list[["col.breaks"]]
    Colour.labs <- Legend.breaks.list[["col.labs"]]
    Colours <- Legend.breaks.list[["cols"]]
    if(rotate){
    	y.coord.breaks <- y.coord.breaks - min(y.coord.breaks)
    	x.coord.breaks <- x.coord.breaks - min(x.coord.breaks)
    	if(length(Legos)==2){
    		Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0,]
    		Lower.tri.map <- Matrix.df[Matrix.df$dist <= 0,]
    		Lower.tri.map$dist <- abs(Lower.tri.map$dist)
    		Upper.rotated.map <- RotateHeatmap(Matrix=Upper.tri.map, value.var="rescale", upper = TRUE)
    		Lower.rotated.map <- RotateHeatmap(Matrix=Lower.tri.map, value.var="rescale", upper = FALSE) 
    		Entire.rotated.map <- rbind(Upper.rotated.map,Lower.rotated.map)
    		y.coord.breaks <- y.coord.breaks/2
    		y.coord.breaks <- c(rev(y.coord.breaks)*-1,y.coord.breaks)
			y.axis.coord.labs <- c(rev(y.axis.coord.labs),y.axis.coord.labs)
    	}else{
    		Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0,]
    		Entire.rotated.map <- RotateHeatmap(Matrix=Upper.tri.map, value.var="rescale", upper = TRUE)
    		y.coord.breaks <- y.coord.breaks/2
    	}
    }

	# require(ggplot2)
	Lego_theme <- Get_heatmap_theme(x.axis=x.axis, y.axis=y.axis, text.size = text.size, 
		x.axis.text.size = x.axis.text.size, y.axis.text.size = y.axis.text.size,
		legend.title.text.size = legend.title.text.size, legend.text.size = legend.text.size, 
		title.size = title.size, legend.key.width = legend.key.width, legend.key.height =legend.key.height)
	Labels <- Get_heatmap_titles(title = title, x.axis.title = x.axis.title, y.axis.title = y.axis.title, 
		legend.title = legend.title, x.coords = x.coords, y.coords = y.coords, rotate = rotate)
    Boundaries.obj <- NULL
    if(rotate){
		if(!is.null(tad.ranges)){
		    Boundaries.obj <- Format_boundaries_rotated_heatmap(Legos = Legos, Ranges = tad.ranges, group.col = group.col, 
			cut.corners = cut.corners, colour.col = tad.colour.col, colours = colours, 
			colours.names = colours.names, region.chr = x.coord.parsed['chr'], 
			region.start = min(Matrix.df[,'row']), region.end = max(Matrix.df[,'row']), distance = distance)
		}
	    ThePlot <- ggplot(Entire.rotated.map, aes(x = xcoords, y = ycoords))
	    ThePlot <- ThePlot + geom_polygon(aes(fill = values, group = ids))
	    xlims <- c(0,max(Entire.rotated.map[,"xcoords"]))
	    ylims <- c(min(Entire.rotated.map[,"ycoords"]),max(Entire.rotated.map[,"ycoords"]))
	    y.coord.breaks <- seq(ceiling(min(Entire.rotated.map[,"ycoords"])),ceiling(max(Entire.rotated.map[,"ycoords"])),length.out = y.axis.num.breaks)
    	y.axis.coord.labs <- y.coord.breaks*2
    }else{
		if(!is.null(tad.ranges)){
		    Boundaries.obj <- Format_boundaries_normal_heatmap(Legos = Legos, Ranges = tad.ranges, group.col = group.col, 
			cut.corners = cut.corners, colour.col = tad.colour.col, colours = colours, 
			colours.names = colours.names, region.chr = x.coord.parsed['chr'], 
			region.start = min(Matrix.df[,'row']), region.end = max(Matrix.df[,'row']), distance = distance)
		}
	    ThePlot <- ggplot(Matrix.df, aes(x = row, y = col))
	    ThePlot <- ThePlot + geom_tile(aes(fill = rescale))
	    ThePlot <- ThePlot + Boundaries.obj
	    xlims <- c(min(Matrix.df$row),max(Matrix.df$row))
	    ylims <- c(min(Matrix.df$col),max(Matrix.df$col))
    }

    ThePlot <- ThePlot + scale_x_continuous(limits = xlims, expand = c(0,0),
    	breaks = x.coord.breaks, labels = x.axis.coord.labs)
    ThePlot <- ThePlot + scale_y_continuous(limits = ylims, expand = c(0,0),
    	breaks = y.coord.breaks, labels = y.axis.coord.labs)
    ThePlot <- ThePlot + scale_fill_gradientn(legend.title,values = Value.dist, 
    	breaks = Colour.breaks, labels = Colour.labs, colors = Colours)
    ThePlot<-ThePlot+ Lego_theme
    # return(Entire.rotated.map)
    ThePlot<-ThePlot+labs(title = Labels['title'], x = Labels['x.axis'], y = Labels['y.axis'])
    ggsave(filename = File, plot = ThePlot, width = width, height = height, units = units)
    if(return.object){
	    return(ThePlot)
    }
}






# GenomicVizard <- R6Class("HiCPlotter",
# 	public = list (
# 		ChromosomeList=NA,
# 		initialize = function(){
# 			require(Viridis)
# 			require(RColorBrewer)
# 			private$Viridis.pal <- c("viridis","magma","plasma","inferno")
# 			private$Viridis.funs <- list()
# 			private$Viridis.funs[["viridis"]] <- viridis
# 			private$Viridis.funs[["magma"]] <- magma
# 			private$Viridis.funs[["plasma"]] <- plasma
# 			private$Viridis.funs[["inferno"]] <- inferno
# 			private$ColorBrewer.seq.pal <- c("BuGn" , "BuPu" , "GnBu" , "OrRd" , "PuBu" , "PuRd" , "RdPu" , "YlGn" , "PuBuGn" , "YlGnBu" , "YlOrBr" , "YlOrRd" , "Blues" , "Greens" , "Greys" , "Oranges" , "Purples" , "Reds")
# 			private$ColorBrewer.seq.pal.lim <- c(9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
# 			names(private$ColorBrewer.seq.pal.lim) <- private$ColorBrewer.seq.pal

# 			private$ColorBrewer.cat.pal <-c("Accent" , "Dark2" , "Paired" , "Pastel1" , "Pastel2" , "Set1" , "Set2" , "Set3")
# 			private$ColorBrewer.cat.pal.lim <- c(8, 8, 12, 9, 8, 9, 8, 12)
# 			names(private$ColorBrewer.cat.pal.lim) <- private$ColorBrewer.cat.pal

# 			private$ColorBrewer.div.pal <-c("BrBG" , "PiYG" , "PRGn" , "PuOr" , "RdBu", "RdGy" , "RdYlBu" , "RdYlGn" , "Spectral")
# 			private$ColorBrewer.div.pal.lim <- c(11, 11, 11, 11, 11, 11, 11, 11, 11)
# 			names(private$ColorBrewer.div.pal.lim) <- private$ColorBrewer.div.pal
# 			Brewer.limits <- list()
# 			Brewer.limits[["sequential"]] <- private$ColorBrewer.seq.pal.lim
# 			Brewer.limits[["categorical"]] <- private$ColorBrewer.cat.pal.lim
# 			Brewer.limits[["diverging"]] <- private$ColorBrewer.div.pal.lim

# 			private$Reference.type <- list()
# 			private$Reference.type[["Heatmap"]] <- list()
# 			private$Reference.type[["Heatmap"]][["sequential"]] <-private$ColorBrewer.seq.pal 
# 			private$Reference.type[["Heatmap"]][["categorical"]] <-private$ColorBrewer.cat.pal
# 			private$Reference.type[["Heatmap"]][["diverging"]] <-private$ColorBrewer.div.pal
# 			private$Reference.type[["TwoSampleHeatmap"]] <- list()
# 			private$Reference.type[["TwoSampleHeatmap"]][["sequential"]] <- private$ColorBrewer.seq.pal
# 			private$Reference.type[["TwoSampleHeatmap"]][["categorical"]] <- private$ColorBrewer.cat.pal
# 			private$Reference.type[["TwoSampleHeatmap"]][["diverging"]] <- private$ColorBrewer.div.pal
# 			private$Reference.type[["TwoSampleHeatmap"]] <- list()
# 			private$Reference.type[["PEScan"]] <- list()
# 			private$Reference.type[["PEScan"]][["sequential"]] <- private$ColorBrewer.seq.pal
# 			private$Reference.type[["PEScan"]][["diverging"]] <- private$ColorBrewer.div.pal
# 			private$Reference.type[["SignalDecay"]] <- list()
# 			private$Reference.type[["SignalDecay"]][["categorical"]] <- private$ColorBrewer.cat.pal
# 		},
# 		PlotSimpleHeatmap = function(GenomicMatrix=NULL,XCoords=NULL,YCoords=NULL,type="within",
# 			palette="inferno",color.interpolate.on="palmax",color.interpolate.to=100,na.fill="#FFFFFF",
# 			midpoint=NULL,facet=FALSE,facet.groups=NULL,facet.order=NULL,signal.cap=NULL,FUN=NULL,
# 			na.to.zero=TRUE,rotate=FALSE,rotated.distance=NULL){
# 			private$Figure.type=private$SetType(Type="Heatmap")
# 			private$Rotate <- rotate
# 			if(!is.null(rotated.distance) & (!(class(rotated.distance) %in% c('numeric','integer')) | length(rotated.distance)!=1)){
# 				stop("If rotated.distance is NOT NULL a numeric integer value of length 1 is expected!\n")
# 			}
# 			# Plots a simple symmetric/assymmetric heatmap
# 			if(facet){
# 				if(class(GenomicMatrix) != "list"){
# 					stop("If facet is TRUE GenomicMatrix must be a list containing objects of class GenomicMatrix")
# 				}
# 				if(is.null(facet.groups)){
# 					if(is.null(names(GenomicMatrix))){
# 						stop("When facet.groups is NULL GenomicMatrix names cannot be NULL")
# 					}
# 					facet.groups=names(GenomicMatrix)
# 				}
# 				private$facet.groups <- facet.groups
# 				if(is.null(facet.order)){
# 					facet.order <- c(1:length(GenomicMatrix))
# 				}
# 				if(!(length(GenomicMatrix)==length(facet.groups) & length(facet.groups)==length(facet.order))){
# 					stop("Length of list doesn't match the length of facet.groups and facet.order")
# 				}
# 				private$facet <- facet
# 				private$facet.order <- facet.order
# 				private$facet.groups <- facet.groups
# 			}else{
# 				private$facet <- FALSE
# 				GenomicMatrix.list <- list(GenomicMatrix)
# 				private$facet.order <- 1
# 				private$facet.groups <- NA
# 			}
# 			Colours <- private$CheckAndCreateColors(palette=palette, on=color.interpolate.on,
# 				to=color.interpolate.to, midpoint=midpoint)
# 			if(is.null(GenomicMatrix) | !(class(GenomicMatrix) %in% c("GenomicMatrix"))){
# 				stop("An object of class GenomicMatrix must be provided")
# 			}
# 			if(!na.to.zero){
# 				warning("NAs will not be set to zero.\n")
# 				tryCatch(class(col2rgb(na.fill)),error=function(e){stop("Could not coerce colour to rgb")})				
# 			}
# 			if(!is.null(FUN)){
# 				Seed <- c(10,20,30)
# 				TryLen <- tryCatch(length(FUN(Seed)),error=function(e){stop("Could not apply FUN\n")})
# 				if(TryLen!=length(Seed)){
# 					stop("FUN transformed values returns a length different from the original vector\n")
# 				}
# 			}
# 			Matrix.list <- lapply(private$facet.order,function(iter){
# 				Object <- GenomicMatrix.list[[iter]]
# 				Current.matrix <- private$RetrieveMatrices(GenomicMatrix=Object,XCoords=XCoords,YCoords=YCoords,
# 					group=private$facet.groups[iter],na.to.zero=na.to.zero,FUN=FUN)
# 				OUT <- Current.matrix
# 				if(private$Rotate){
# 					Matrices <- private$RotateHeatmap(Matrix=Matrices,rotated.distance=rotated.distance,triangle=1,id.seed=iter)
# 					OUT <- Matrices
# 				}
# 				return(OUT)
# 			})
# 			Matrices <- do.call(rbind,Matrix.list)

# 			Matrices$value[is.infinite(Matrices$value)] <- 0
# 			Matrices <- private$SetGroupLevels(x=Matrices)

# 			private$PlottingData <- Matrices
# 			self$BuildHeatmap(axis.x.breaks=)
# 		},
# 		PlotTwoSampleHeatmap = function(GenomicMatrix.upper=NULL,GenomicMatrix.lower=NULL,
# 			XCoords=NULL,YCoords=NULL,palette.package="Viridis",
# 			palette="inferno",facet=FALSE,signal.cap=NULL,FUN=NULL,rotate=FALSE){
			
# 			private$Figure.type="TwoSampleHeatmap"
# 			# Plots a simple symmetric/assymmetric heatmap, but produces 
# 			# a faceted plot highlighting regions: TADs, peaks, regions in HiCData
# 			# colorpalette takes arguments Viridis or RColorBrewer
# 			# If it is Rcolorbrewer the palette must be one of the sequential palettes 
# 		},
# 		Plot.PEScan = function(Matrix=NULL,Window=NULL,
# 			palette.package="Viridis",palette="inferno",signal.cap=NULL,FUN=NULL){
# 			private$Figure.type="PEScan"
# 		},
# 		SetType = function(Type=NULL){
# 			if(!(Type %in% TypeList)){
# 				stop("Provided Type is not in allowed types.\nAllowed types are ",private$TypeList)
# 			}
# 		},
# 		BuildHeatmap = function(Matrix = NULL, axis.x.breaks = NULL, axis.y.breaks = NULL,num.breaks.x=5,num.breaks.y=5,
# 			axis.x.labels = NULL,axis.y.labels = NULL, label.start.x=NULL, label.end.x=NULL, label.start.y=NULL, label.end.y=NULL, 
# 			divide.label.by = NULL, axis.x.title = NULL, axis.y.title = NULL, colours = NULL, colour.breaks = NULL, colourbar.breaks.num = NULL, 
# 			remove.colourbar = FALSE, remove.x.axis = FALSE, remove.y.axis = FALSE, remove.plot.title = TRUE, 
# 			na.to.zero = TRUE, na.fill = "#FFFFFF", FUN = NULL,signal.cap=NULL){
# 			# Central function which builds heatmaps for every other function
# 			# Control statements 
# 			require(reshape2)
# 			require(ggplot2)
# 			require(RColorBrewer)
# 			require(scales)
# 			if(!is.null(Matrix)){
# 				if(class(Matrix)!="matrix"){
# 					stop("Matrix must be of class matrix\n")
# 				}
# 				if(class(Matrix)=="matrix"){
# 					rownames(Matrix) <- 1:nrow(Matrix)
# 					colnames(Matrix) <- 1:ncol(Matrix)
# 					Matrix.melt <- melt(Matrix.melt)
# 					colnames(Matrix.melt) <- c('row','col','value')
# 				}			
# 			}
# 			if(na.to.zero){
# 				Matrix.melt$value[is.na(Matrix.melt$value)] <- 0
# 			}else{
# 				tryCatch(class(col2rgb(na.fill)),error=function(e){stop("Could not coerce colour to rgb")})
# 			}
# 			Matrix.melt$value[is.infinite(Matrix.melt$value)] <- 0
# 			if(!is.null(FUN)){
# 				data.transform <- FUN(Matrix.df$value)
# 				Matrix.df$value <- data.transform
# 			}
# 			NAFilter <- !is.na(Matrix.df$value)
# 			if(!is.null(signal.cap)){
# 				if(signal.cap > 1 | signal.cap < 0.5 | length(signal.cap) > 1){
# 					stop("signal.cap takes only length 1 values between 0.5 and 1.\n 
# 						I'm sorry I consider setting a signal cap any smaller as a bad design choice.")
# 				}
# 				Quantile <- quantile(Matrix.df$value[NAFilter],signal.cap)
# 				Matrix.df$value[Matrix.df$value>Quantile & NAFilter] <- Quantile
# 			}
# 			main.dat <- Matrix.df$value[NAFilter]
# 			if (min(Matrix.df$value[NAFilter]) != 0 | max(Matrix.df$value[NAFilter]) != 1){
# 				Rescale <- rescale(Matrix.df$value,c(0,1))
# 				Matrix.df$value <- Rescale
# 			}
# 			# Insert conditionals for rotation later on
# 			# if()
# 			# {
# 			Matrix.df[,"row.adjust"] <- Matrix.df[,"row"] - 0.5
# 			Matrix.df[,"col.adjust"] <- Matrix.df[,"col"] - 0.5
# 			# }
# 			if(!remove.x.axis){
# 				if(is.null(axis.x.breaks)){
# 					if(is.null(num.breaks.x) | length(num.breaks.x) > 1 | num.breaks.x<=0 | !(class(num.breaks.x) %in% c("numeric","integer"))){
# 						stop("num.breaks.x must be a numeric non zero positive value of length 1")
# 					}
# 					axis.x.breaks <- seq(from=min(Matrix.df[,"row.adjust"])-0.5,to=max(Matrix.df[,"row.adjust"])+0.5),length.out=num.breaks.x)
# 					cat("x axis breaks: ",paste(axis.x.breaks,collapse = ", "),"\n")
# 				}
# 				if(!is.null(axis.x.labels)){
# 					if(is.null(label.start.x) | is.null(label.end.x)){
# 						stop("label.start.x and label.end.x cannot be null when remove.x.axis is FALSE\n")
# 					}
# 					axis.x.labels <- seq(from=label.start.x,to=label.end.x,length.out=num.breaks.x)
# 					if(!is.null(divide.label.by)){
# 						axis.x.labels <- axis.x.labels/divide.label.by
# 					}
# 					cat("x axis labels: ",paste(axis.x.labels,collapse = ", "),"\n")
# 				}
# 			}
# 			# if axis labels are null then don't print axis titles.
# 			if(!remove.y.axis){
# 				if(is.null(axis.y.breaks)){
# 					if(is.null(num.breaks.y) | length(num.breaks.y) > 1 | num.breaks.y<=0 | !(class(num.breaks.y) %in% c("numeric","integer"))){
# 						stop("num.breaks.y must be a numeric non zero positive value of length 1")
# 					}
# 					axis.y.breaks <- seq(from=min(Matrix.df[,"col.adjust"])-0.5,to=max(Matrix.df[,"col.adjust"])+0.5),length.out=num.breaks.y)
# 					cat("y axis breaks: ",paste(axis.x.breaks,collapse = ", "),"\n")
# 				}
# 				if(is.null(axis.y.labels)){
# 					if(is.null(label.start.y) | is.null(label.end.y)){
# 						stop("label.start.y and label.end.y cannot be null when remove.y.axis is FALSE\n")
# 					}
# 					axis.y.labels <- seq(from=label.start.y,to=label.end.y,length.out=num.breaks.y)
# 					if(!is.null(divide.label.by)){
# 						axis.y.labels <- axis.y.labels/divide.label.by
# 					}
# 					cat("y axis breaks: ",paste(axis.y.labels,collapse = ", "),"\n")
# 				}
# 			}
# 			if(is.null(colour.breaks) | is.null(colour)){
# 				col.break.num <- 100
# 				colour.breaks <- seq(from=0,to=1,length.out=col.break.num)
# 				colours <- colorRampPalette(brewer.pal(name = "Greys", n = 9))(col.break.num+1)
# 			}
# 			if(!remove.colourbar){
# 				if(is.null(colourbar.breaks.num)){
# 					colourbar.breaks.num <- 5
# 				}
# 				colourbar.breaks <- seq(from=0,to=1,length.out=colourbar.breaks.num)
# 				colourbar.labels <- seq(from=min(main.dat),to=max(main.dat),length.out=colourbar.breaks.num)
# 			}

# 		}
# 		PlotDecayOfSignal = function(GenomicMatrix=NULL,Distances=NULL,
# 			palette="Pastel1",FUN.transform=NULL,FUN.aggregate=NULL){

# 		},		
# 		AddTADBoundaries = function(Chrom=NULL,Start=NULL,End=NULL,Group=NULL,
# 			palette.package="Viridis",palette="inferno",size=1){
# 		},
# 		AddHighlightRegions = function(Chrom=NULL,Start=NULL,End=NULL,Group=NULL,
# 			palette.package="Viridis",palette="inferno",size=1){

# 		},
# 	),
# 	private = list(
# 		Type=NA,
# 		TypeList=c("TwoSampleHeatmap","RotatedHeatmap",),
# 		Viridis.pal = NA,
# 		ColorBrewer.seq.pal = NA,
# 		ColorBrewer.seq.pal.lim = NA,
# 		ColorBrewer.cat.pal = NA,
# 		ColorBrewer.cat.pal.lim =NA,
# 		ColorBrewer.div.pal = NA,
# 		ColorBrewer.div.pal.lim = NA,
# 		ColourNAs = "#FFFFFF",
# 		Viridis.funs = NA,
# 		Figure.type = NA,
# 		Reference.type = NA,
# 		Brewer.limits = NA,
# 		facet = FALSE,
# 		facet.order = NA,
# 		facet.groups = NA,
# 		Rotate = NA,
# 		PlottingData=NA,
# 		CheckAndCreateColors = function(palette=NULL,on=NULL,to=NULL,midpoint=NULL){
# 			if((palette %in% private$Viridis.pal)){
# 				cat("Viridis based colours will ignore the variable color.interpolate.on\n")
# 				require(Viridis)
# 				WhichFun <- which(private$Viridis.pal %in% palette)
# 				if(!(class(to) %in% c("integer","numeric"))){
# 					stop("number of colors to create must be a numeric integer\n")
# 				}
# 				Colours <- private$Viridis.funs[[WhichFun]](n=to)
# 				return(Colours)
# 			}
# 			CheckPalette(palette=palette,on=on,to=to,midpoint=midpoint)
			
# 			Base.colours <- brewer.pal(n=on, name=palette)
# 			if(!is.null(midpoint)){
# 				Lower.tail <- Base.colours[1:ceiling(length(Base.colours)/2)]
# 				Upper.tail <- Base.colours[ceiling(length(Base.colours)/2):length(Base.colours)]
# 				Lower.colours <- colorRampPalette(Lower.tail)(ceiling(to/2))
# 				Upper.colours <- colorRampPalette(Upper.tail)(ceiling(to/2))
# 				Colours <- unique(c(Lower.colours,Upper.colours))
# 			}else{
# 				Colours <- colorRampPalette(Base.colours)(to)
# 			}
# 			return(Colours) 
# 		},
# 		CheckPalette = function(palette=NULL,on=NULL,to=NULL,midpoint=NULL){
# 			WhichPaletteType <- which(sapply(private$Reference.type[[private$Figure.type]],function(x){
# 							palette %in% x[[1]]
# 						}))
# 			if(!(on %in% c("palmax")) & !(class(on) %in% c("numeric","integer"))){
# 				stop("color.interpolate.on takes either palmax or numeric integer values ","\n")				
# 			}
# 			if(length(WhichPaletteType)==0){
# 				stop("palette does not match brewer types allowed for the particular plot type.\n",
# 					"Allowed palettes for the ",private$Figure.type," are ",
# 					paste(names(private$Reference.type[[private$Figure.type]]),collapse=", "),"\n")
# 			}

# 			NameOfPaletteType <- names(private$Reference.type[[type]][WhichPaletteType])
# 			if(!is.null(midpoint) & WhichPaletteType!="diverging"){
# 				stop("When argument midpoint is supplied, only diverging brewer palettes can be used.")
# 			}
# 			if(!is.null(midpoint) & (on %% 2 == 0)){
# 				stop("when argument midpoint is supplied, on takes only odd values.")
# 			}
# 			WhichPalette <- which(private$Reference.type[[type]][[WhichPaletteType]] %in% palette)
# 			if(on=="palmax"){
# 				on <- private$Brewer.limits[[NameOfPaletteType]][WhichPalette]
# 			} else if(on > private$Brewer.limits[[WhichPalette]][WhichPalette]){
# 				stop("Number of colours available in the brewer palette ",names(private$Brewer.limits[[WhichPalette]][WhichPalette]),"is ",
# 				private$Brewer.limits[[WhichPalette]][WhichPalette],"\n")
# 			}
# 			if(NameOfPaletteType=="diverging" & is.null(midpoint)){
# 				stop("For diverging scales a midpoint must be defined\n")
# 			}
# 		}
# 		RetrieveMatrices = function(GenomicMatrix=NULL,XCoords=NULL,YCoords=NULL,group=NULL,na.to.zero=NULL,FUN=NULL){
# 			if(class(GenomicMatrix)!="GenomicMatrix"){
# 				stop("Provided object must be of class GenomicMatrix! This function is not Iterable")
# 			}
# 			Matrix <- GenomicMatrix$FetchMatrixWithinCoords(XCoords=XCoords,YCoords=YCoords,type=type)
# 			rownames(Matrix) <- c(1:nrow(Matrix))
# 			colnames(Matrix) <- c(1:ncol(Matrix))
# 			require(reshape2)
# 			Matrix.melt <- melt(Matrix)
# 			colnames(Matrix.melt) <- c('row','col','value')
# 			if(na.to.zero){
# 				Matrix.melt$value[is.na(Matrix.melt$value)] <- 0
# 			}
# 			Matrix.melt$value[is.infinite(Matrix.melt$value)] <- 0
# 			if(!is.na(group)){
# 				Matrix.melt$group <- group
# 			}
# 			if(!is.null(FUN)){
# 				Non.NA.filter <- !is.na(Matrix.melt$value)
# 				Replace.values <- FUN(Matrix.melt$value[Non.NA.filter])
# 				Matrix.melt$value[Non.NA.filter] <- Replace.values				
# 			}
# 			return(Matrix.melt)
# 		},
# 		SetGroupLevels = function(x=NULL){
# 			if("group" %in% colnames(x)){
# 				x$group <- factor(x$group, levels=private$facet.groups[private$facet.order])
# 			}
# 			return(x)
# 		},
# 		RotatePixels <- function(shift = NULL, plot.order = NULL, vector = NULL, distance = NULL,triangle=1,id.seed=NULL){
# 		    leftmost.x <- seq((0+shift),length.out=length(vector),by=1)
# 		    leftmost.y <- rep((0+shift),length(vector)) * triangle
# 		    Order<-seq(plot.order,length.out=length(leftmost.x),by=4)

# 		    uppermid.x <- seq((0.5+shift),length.out=length(vector),by=1)
# 		    uppermid.y <- rep((0.5+shift),length.out=length(vector),by=1) * triangle
# 		    Order<-c(Order,seq(min(Order)+1,length.out=length(uppermid.y),by=4))
		  
# 		    rightmost.x <- seq((1+shift),length.out=length(vector),by=1)   
# 		    rightmost.y <- rep((0+shift),length(vector)) * triangle
# 		    Order<-c(Order,seq(min(Order)+2,length.out=length(uppermid.y),by=4))

# 		    lowermid.x <- leftmost.x + 0.5
# 		    lowermid.y <- ((leftmost.y * triangle) - 0.5) * triangle
# 		    Order <- c(Order,seq(min(Order)+3,length.out=length(uppermid.y),by=4))

# 		    IDs <- factor(paste("Pixel",rightmost.x,plot.order,id.seed,sep="."),levels=paste("Pixel",rightmost.x,plot.order,id.seed,sep="."))
# 		    IDsRep <- rep(IDs,4)

# 		    horizontal.coords <- data.frame(xcoords = c(leftmost.x, uppermid.x, rightmost.x, lowermid.x),
# 		        ycoords = c(leftmost.y, uppermid.y, rightmost.y, lowermid.y), ids = IDs, Distance = distance, Order=Order, group = id.seed)

# 		    horizontal.values <- data.frame(ids=IDs,values=vector)

# 		    horizontal.data <- merge(horizontal.values, horizontal.coords, by=c("ids"))
# 		    horizontal.data <- horizontal.data[order(horizontal.data$Order),]
# 		    return(horizontal.data)
# 		},
# 		RotateHeatmap = function(Matrix=NULL,rotated.distance=NULL,triangle=1,id.seed=NULL){
# 			Max.row <- max(Matrix$row)
# 			Max.col <- max(Matrix$col)
# 			Max.extent <- Max.col - Max.row
# 			if(is.null(rotated.distance)){
# 				rotated.distance <- Max.extent
# 			}
# 			if(rotate.distance > Max.extent){
# 				stop("Provided rotated.distance was larger than matrix extent.\n")
# 			}
# 			Matrix$dist <- Matrix[,'col'] - Matrix[,'row']
# 			Rotated.df.list <- lapply(1:rotated.distance,function(x){
# 			    Distance <- x - 1
# 			    Vector <- Matrix$value[Matrix$dist==x]
# 			    Shift <- 0.5 * Distance
# 			    Dataframe <- private$RotatePixels(shift = Shift, plot.order = 1, vector = Vector, distance = Distance,triangle=triangle,id.seed=id.seed)
# 			})
# 			Rotated.df <- do.call(rbind,Rotated.df.list)
# 			return(Rotated.df)
# 		}
# 	)
# )