# The file that contains observation paths
inputTxt = "swift.txt"

# If RA and DEC are set to RA/DEC = "", the pipeline will look for OBJ_RA and OBJ_DEC keywords to find values for RA/DEC by default.
RA = ""
DEC = ""

# When 'pileup' is set to True, source region will be an annulus, with inner radius = 'ignore_radius' and outer radius = 'sourceRadius' pixels
pileup = False
ignore_radius = 3  # In pixels (1 pixel = 2.36 arcseconds)

# If set to True, script will recalculate source coordinates, even if RA and DEC parameters are given above.
# For WT, brightest pixel will be assumed to be the center. For PC, 'xrtcentroid' command will be run.
# The calculated coordinates, especially when the photon count is low, may be incorrect; thus do not forget to check whether the region is correct.
# Set it to False when: you provide correct RA and DEC above, or want to use default RA and DEC in files, or setting it to True gives wrong results
# Set it to True when: You do not have correct RA and DEC values and default values are not correct either
recalculate_source_center = False

# Radius of the source region
sourceRadius = 20

#Background annulus region for windowed timing mode
innerRadius = 90
outerRadius = 110