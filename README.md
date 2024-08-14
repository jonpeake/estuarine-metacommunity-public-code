# Assembly processes in estuarine metacommunities are dependent on spatial scale

This is the collective data and R code for Jonathan Peake's work on fish metacommunity assembly in Florida estuaries.

The data file samples.csv includes data on 183-m haul seine catches and their associated habitat, spatial, temporal, and physical data. Each row corresponds to one haul seine sample. The dataset contains the following columns:

1.  Columns prefixed with "Bio\_" contain catch data by species for each haul seine; the trailing numeric string on each of these columns refers to the NODC code for each species (see the species_codes.csv file for species metadata). Catches are reported in total number of each species.
2.  Bay: denotes the estuary where each seine was conducted. Their abbreviations are as follows:
    1.  AP - Apalachicola Bay
    2.  CK - Cedar Key
    3.  TB - Tampa Bay
    4.  CH - Charlotte Harbor
    5.  JX - northeast Florida
    6.  IR - northern Indian River Lagoon
    7.  TQ - southern Indian River Lagoon
3.  StartDepth: the depth of the seine in meters as taken at the center bag
4.  vis: vertical visibility; calculated as the Secchi depth divided by the starting depth (unitless)
5.  Temperature: measured in degrees Celsius at the center bag
6.  pH: measured at the center bag (unitless)
7.  Salinity: measured in PSU at the center bag
8.  DissolvedOxygen: measured in mg/L at the center bag
9.  Columns prefixed with "Bottom\_" contain presence data of bottom types encountered during a seine sampling event; the trailing single-character string for each of these columns refers to the bottom type. See the benthic_codes.csv file for details.
10. Columns prefixed with "Veg\_" contain proportional data of vegetation encountered during a seine sampling event; the trailing double-character string for each of these columns refers to the vegetation classification. See the vegetation_codes.csv file for details. Values reported represent the percent coverage of each of the vegetation types within the sample area.
11. Columns prefixed with "Bycatch\_" contain catch data of non-target benthic invertebrates, algae, leaf litter, and algae for each seine sample; the double-character string for each of these columns refers to the bycatch classification. See the bycatch_codes.csv file for details. Values reported represent the volume in liters of each bycatch type caught in a given sample.
12. Columns prefixed with "Shore\_" contain presence data of shoreline habitats encompassed by a seine sample; a double-character string refers to the classification of shoreline. A trailing "Over" appended to the shoreline type refers to an overhanging element of the shoreline type. See the shoreline_codes.csv file for details.
13. Latitude: reported in decimal degrees
14. Longitude: reported in decimal degrees
15. Month: month in which the sample was collected (1-12, 1 = January, 12 = December)
16. Year: year in which the sample was collected (1998-2020)

The R files include all code necessary to conduct the analyses presented in the manuscript, commented to indicate the usage of functions and code sections.

For further information, please reach out to Jonathan Peake - jpeake1\@usf.edu