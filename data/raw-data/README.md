---
editor_options: 
  markdown: 
    wrap: 72
---

# README

This folder contains raw data used as inputs of the workflow, including:

-   VertebrateSpecies-list_FoS_AcT_Morph_MovM_Aq_DD_LifeH_Diet_HabP_NHabP_Press.xlsx
    (excel file): Species list (149 birds and 55 mammals) with their
    taxonomy and traits. Column specification:

    -   CLASS: species class
    -   ORDER: species order
    -   TAXREFV16_CDREF: taxref (version 16) code
    -   SPECIES_NAME: scientific species name
    -   SPECIES_NAME_SYNONYM: scientific species synonym name
    -   FORAG.STRAT.X: foraging strategy X (binary variable)
    -   ACT.TIME.X: activity time X (binary variable)
    -   MORPHO.BODYMASS.G: body mass in grams
    -   MOV.MODE.X: movement mode X (binary variable)
    -   AQUATIC: aquatic habitat dependence (binary variable)
    -   DISPERSAL_KM: mean natal dispersal distance in kilometers per
        dispersal event
    -   LIFE.HIST.OFFSPRING_PER_YEAR_N: number of offsprings produced
        per year
    -   DIET.BREADTH: number of diet items eaten
    -   HAB.PREF.BREADTH: number of preferred habitats
    -   NEST.HAB.BREADTH: number of preferred nesting habitats
    -   PRESSURE.X : vulnerability to pressure X (binary variable)

-   EnvironmentalVariables_France_Res1000m (RDS file): Spatialized table
    of environmental variables (554783 pixels of 1km2, each row is a
    pixel x 45 environmental variables, each column is a variable). Tif
    and gpkg files are also provided. Variable description:

    -   climatic.bio01d : annual mean near-surface air temperature in
        degree Celsius per pixel
    -   climatic.bio04d : mean daily precipitations in kg/m2 per pixel
    -   climatic.bio12d : temperature seasonality per pixel
    -   climatic.bio15d : precipitation seasonality per pixel
    -   lin.struct.HedgeProportion : proportion of hedges per pixel
        (length of hedges / pixel area)
    -   lin.struct.LakeProportionImportance1to4: proportion of important
        lakes (1 to 4) per pixel (surface of lakes / pixel area)
    -   lin.struct.RailwayProportion: proportion of railways per pixel
        (length of railways / pixel area)
    -   lin.struct.RoadProportionImportance1: proportion of roads
        (importance 1) per pixel (length of roads / pixel area)
    -   lin.struct.RoadProportionImportance2: proportion of roads
        (importance 2) per pixel (length of roads / pixel area)
    -   lin.struct.RoadProportionImportance3: proportion of roads
        (importance 3) per pixel (length of roads / pixel area)
    -   lin.struct.RoadProportionImportance4: proportion of roads
        (importance 4) per pixel (length of roads / pixel area)
    -   lin.struct.RoadProportionImportance5: proportion of roads
        (importance 5) per pixel (length of roads / pixel area)
    -   lin.struct.StreamProportion: proportion of large streams (\>= 15
        meters) per pixel (length of streams / pixel area)
    -   topo.MeanSlope: pixel mean slope
    -   topo.CValti: coefficient of variation of altitude per pixel
    -   land.compo.WAW_prop: proportion of water and wetness per pixel
        (surface of water and wetness / pixel area)
        (<https://land.copernicus.eu/en/products/high-resolution-layer-water-and-wetness>)
    -   land.compo.X: binary variable indicating whether land category X
        is the land category of the pixel. Code specification:

![](images/TableofLandSystCode.png)

-   Vertebrate-Species-GBIF-INPN-IUCNOccurrenceData_France_Res1000m_2010-2020
    (RDS file): Spatialized table of presences / pseudo-absences (554783
    pixels of 1km2 pixel, each row is a pixel x 204 species, each column
    is the presence/pseudo-absence distribution for a species). Pixel
    order in that table is the same as the pixel order in
    EnvironmentalVariables_France_Res1000m. Note that Muscicapa
    tyrrhenica is treated as Muscicapa striata so there are only 203
    columns.

-   ProtectedAreas folder:

    -   Shapefile of strict protection delineation
        (STRICT_PROTECTIONS.shp).

    -   Matrix of Euclidean distance among strict protections, in meters
        (Distance_btw_StrictProtections, RDS file).

-   Grids folder:

    -   Geopackage and raster files of the reference grid (1km^2^) used
        to project occurrences and environmental variables.

-   SamplingEffort_France_Aves_Res1000_2010-2020.tif: raster of
    occurrence sampling effort for Aves.

-   SamplingEffort_France_Mammalia_Res1000_2010-2020.tif: raster of
    occurrence sampling effort for Mammalia.
