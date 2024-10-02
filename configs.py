
# Set a name for this model run (no protected chacacters, including the date is a good idea)
run_name = "Test_06122023"

# set time range and resolution
res = {
    "ystart": 1800, # Years to model
    "yend": 2023,
    "lonstart": -180, # Longitude range and resolution
    "lonend": 180,
    "lonres": 0.5,
    "latstart": -60, # Latitude range and resolution
    "latend": 90,
    "latres": 0.5,
}

# Set time range and resolution
# Temporal frequency is different at different times, so you can give a list of start, end and frequency values
mod_obs_time = {
    "ystart": [1740,1940], # Years to start: Should be a list
    "yend": [1940,2022], # Years to finish: Should be a list
    "yres": [25,2], # Temporal frequency in each period. Note: (yend - ystart)_i must divide by yres_i!!
}


