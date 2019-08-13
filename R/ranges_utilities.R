.prepare_ranges_df_for_add <- function(ranges = NULL){
    Reference.object <- GenomicMatrix$new()
    Ranges.df <- as.data.frame(ranges)
    Which.factor <- which(vapply(seq_len(ncol(Ranges.df)), function(x){
            is.factor(Ranges.df[,x])
        }, TRUE))
    Ranges.df[,Which.factor] <- vapply(Which.factor,function(x){
        as.character(Ranges.df[,x])
    },rep("a",nrow(Ranges.df)))
    Metadata.list <- .prepare_ranges_metadata_list(Ranges.df)
    Ranges.df.coords <- Ranges.df[,c(1,2,3)]
    colnames(Ranges.df.coords) <- Reference.object$NonStrandedColNames
    return(list(ranges = Ranges.df.coords, metadata = Metadata.list))
}

.prepare_ranges_metadata_list <- function(ranges_df = NULL){
    Metadata.Cols <- names(ranges_df)[c(4:ncol(ranges_df))]
    Metadata.list <- lapply(Metadata.Cols,function(x){
        if(is.factor(ranges_df[,x])){
            as.character(ranges_df[,x])
        }else{
            ranges_df[,x]
        }
    })
    names(Metadata.list) <- Metadata.Cols
    return(Metadata.list)
}


.send_message <- function(message, issue_stop = TRUE, issue_warning = FALSE){
    if(issue_stop){
        stop(message)
    }else if(issue_warning){
        warning(message)
    }
}

.check_if_rangekey_exists_and_do_something <- function(Brick, rangekey, 
    resolution = NA, all_resolutions = FALSE, issue_stop = TRUE, 
    issue_warning = FALSE, message){
    if(Brick_rangekey_exists(Brick = Brick, rangekey = rangekey, 
        resolution = resolution, all_resolutions = all_resolutions)){
        .send_message(message = message, issue_stop = issue_stop, 
            issue_warning = issue_warning)
    }
}

.check_if_rangekey_not_exists_and_do_something <- function(Brick, rangekey, 
    resolution = NA, all_resolutions = FALSE, issue_stop = TRUE, 
    issue_warning = FALSE, message){
    if(!Brick_rangekey_exists(Brick = Brick, rangekey = rangekey, 
        resolution = resolution, all_resolutions = all_resolutions)){
        .send_message(message = message, issue_stop = issue_stop, 
            issue_warning = issue_warning)
    }
}

.prepare_ranges_metadata_mcols <- function(Brick = NULL, rangekeys = NULL){
    Reference.object <- GenomicMatrix$new()
    mcol.list <- lapply(rangekeys,function(x){
        Handler <- ._Brick_Get_Something_(
            Group.path = Create_Path(c(Reference.object$hdf.ranges.root, x)),
            Brick = Brick, return.what = "group_handle")
        GroupList <- h5ls(Handler,
            datasetinfo = FALSE,
            recursive = FALSE)[,"name"]
        H5Gclose(Handler)
        as_tibble(data.frame(rangekey = x,
                    m.col = GroupList, 
                    stringsAsFactors = FALSE))
    })
    mcol.df <- do.call(rbind,mcol.list)
    hdf.ranges.protected.names <- Reference.object$hdf.ranges.protected.names()
    mcol.df.filter <- !(mcol.df$m.col %in% hdf.ranges.protected.names)
    mcol.df <- mcol.df[mcol.df.filter,]
    if(nrow(mcol.df)==0){
        mcol.df <- NA
    }
    return(mcol.df)
}