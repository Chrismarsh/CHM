#pragma  once
#include "triangle_defines.h"

int is_invalid_mean_elevation_diff(vertex triorg, vertex tridest, vertex triapex, double maxtolerance, const double *gt, const struct masked_array* masked_dem, const GDALRasterBandH* dem);