library(terra)
library(ncdf4)
library(parallel)

### roi_file can be a path to a shapefile or a manually created polygon using vect()

terraOptions(memfrac = 0.8) # Fraction of memory to allow terra

tmpdir        <- "/mnt/c/Rwork"
out_dir       <- "/mnt/g/GOSAT_RemoTec/extracted/Australia_NZ"
out_name      <- "/Australia_NZ_GOSATRemoTeC_"
f_list        <- list.files("/mnt/g/GOSAT_RemoTeC", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
roi_file      <- "/mnt/g/Australia_NZL/AUS_NZL_shp/aus_nzl.shp"
land          <- 0 # 0 = land; 1 = ocean
land_var <- "flag_landtype"
qc            <- 0 # 0 = good; 1 = bad
qc_var        <- "xco2_quality_flag"  
notes         <- "This data has been filtered to include only good soundings (Quality Flag = 0) that are land (flag_landtype = 0)"


tmp_create <- function(tmpdir) {
  
  p_tmp_dir <- paste0(tmpdir, "/", as.character(Sys.getpid())) # Process ID
  
  if (!dir.exists(p_tmp_dir)) {
    dir.create(p_tmp_dir, recursive = TRUE)
  }
  
  terraOptions(tempdir = p_tmp_dir)
}
tmp_remove <- function(tmpdir) {
  
  p_tmp_dir <- paste0(tmpdir, "/", as.character(Sys.getpid())) # Process ID
  unlink(p_tmp_dir, recursive = TRUE)
}

clip_nc <- function(input_file, roi_file, out_dir, out_name, land,
                    land_var, qc, qc_var, tmpdir) {
  
  tmp_create(tmpdir)
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  time_s <- Sys.time()
  
  # if needed, vectorize roi shp file for clipping
  if (typeof(roi_file) != "S4"){
    roi <- vect(roi_file)
    roi <- aggregate(roi)
  } else {
    roi <- aggregate(roi_file)
  }
  
  t_data <- nc_open(input_file)
  
  # Get spatial and time
  coords           <- cbind(ncvar_get(t_data, "longitude"), ncvar_get(t_data, "latitude"))
  colnames(coords) <- c("lon", "lat")
  t                <- basename(input_file)
  t                <- substr(t, 12, 17)
  t                <- paste0("20", t)
  t                <- gsub("(\\d{4})(\\d{2})(\\d{2})$","\\1-\\2-\\3",t) # add dashes
  
  # Get variables and transform to vect for clipping to ROI
  df_var                 <- data.frame(lon = ncvar_get(t_data, "longitude"))
  df_var$lat             <- ncvar_get(t_data, "latitude")
  
  df_var$xco2              <- ncvar_get(t_data, "xco2")
  df_var$xco2_uncertainty  <- ncvar_get(t_data, "xco2_uncertainty")
  df_var$xco2_quality_flag <- ncvar_get(t_data, qc_var)
  df_var$land              <- ncvar_get(t_data, land_var)
  
  nc_close(t_data)
  
  # Filter by land and QC
  df_var <- df_var[df_var$land >= land, ]
  # df_var <- df_var[df_var$xco2_quality_flag == qc, ] # data set only contains good quality data

  # Put coords in their own
  coords <- cbind(df_var$lon, df_var$lat)
  df_var <- subset(df_var, select = -c(lon,lat))
  
  # Clip data
  vec      <- vect(coords, atts = df_var, crs = "+proj=longlat +datum=WGS84")
  var_roi  <- intersect(vec, roi)
  
  # If number of soundings > 0, then proceed
  if (nrow(crds(var_roi, df = TRUE)) == 0) {
    message(paste0("File for this date is being skipped as it has 0 soundings for the region: ", t, "\n"))
    
  } else {
    # Build data frame for writing to nc file
    df <- crds(var_roi, df = TRUE)
    
    for (i in 1:length(names(var_roi))) {
      df <- cbind(df, var_roi[[i]])
    }
    
    # kick out
    rm(df_var, vec, var_roi)
    
    invisible(gc())
    
    #### Create NC file ####
    ### Note: When creating point files, use number of points as a dim
    ### rather than lon and lat, and make lon and lat variables.
    ###
    
    # Create dimensions nc file
    elemdim <- ncdim_def("obs", "", seq(1, nrow(df)))
    
    # Dates
    t_num   <- as.numeric(julian(as.Date(t), origin = as.Date("1970-01-01")))
    
    # define variables
    fillvalue     <- -9999
    dlname        <- "time"
    time_def      <- ncvar_def("time", "days since 1970-01-01", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "longitude"
    lon_def       <- ncvar_def("lon", "degrees_east", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "latitude"
    lat_def       <- ncvar_def("lat", "degrees_north", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "xco2"
    xco2_def      <- ncvar_def("xco2", "ppm (1e-6)", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "xco2 uncertainty"
    xco2_un_def   <- ncvar_def("xco2_uncertainty", "ppm (1e-6)", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "xco2_quality_flag"
    xco2_qc_def   <- ncvar_def("xco2_quality_flag", "none", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "flag landtype"
    land_def      <- ncvar_def("flag_landtype", "none", elemdim, fillvalue, dlname, prec = "float")
    
    # create netCDF file and put arrays
    out_f <- paste0(out_dir, out_name, t, ".nc")
    
    ncout <- nc_create(out_f,
                       list(time_def, lon_def, lat_def, xco2_def, xco2_un_def, xco2_qc_def, land_def), 
                       force_v4 = TRUE)
    
    # put variables
    ncvar_put(ncout, time_def, rep(t_num, times = nrow(df)))
    ncvar_put(ncout, lon_def, df$x)
    ncvar_put(ncout, lat_def, df$y)
    ncvar_put(ncout, xco2_def, df$xco2)
    ncvar_put(ncout, xco2_un_def, df$xco2_uncertainty)
    ncvar_put(ncout, xco2_qc_def, df$xco2_quality_flag)
    ncvar_put(ncout, land_def, df$land)
    
    # put additional attributes into dimension and data variables
    ncatt_put(ncout,"lon","axis","X")
    ncatt_put(ncout,"lat","axis","Y")
    ncatt_put(ncout,"time","axis","T")
    
    # add global attributes
    ncatt_put(ncout,0,"title", "GOSAT/ACOS v9 XCO2")
    ncatt_put(ncout,0,"institution", "University of Oklahoma")
    ncatt_put(ncout,0,"source", "Russell Doughty, PhD")
    ncatt_put(ncout,0,"date_created", date())
    ncatt_put(ncout,0,"notes", notes)
    
    # Close input file
    nc_close(ncout)
    
    time_e   <- Sys.time()
    time_dif <- difftime(time_e, time_s)
    
    message(paste0("Saved ", out_f, ". Time elapsed: ", time_dif, "\n"))
  }
  
  tmp_remove(tmpdir)
}

mclapply(f_list, clip_nc, mc.cores = 10, mc.preschedule = FALSE, roi_file = roi_file,
         out_dir = out_dir, out_name = out_name, land = land, land_var = land_var,
         qc = qc, qc_var = qc_var,  tmpdir = tmpdir)