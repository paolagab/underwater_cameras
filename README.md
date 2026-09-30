# underwater_cameras

Code for mapping global marine biodiversity sampling gaps from OBIS, intersecting them with drifting longline fishing effort (Global Fishing Watch), and identifying priority areas and flag states for opportunistic biodiversity monitoring via fishing vessels.

## Overview

Ocean biodiversity monitoring remains geographically uneven, leaving large regions undersampled. This pipeline:

1.  Builds a global 1° ocean grid (Mollweide equal-area projection) from GEBCO bathymetry and Natural Earth coastlines.
2.  Quantifies historical OBIS sampling density and classifies undersampled areas (grid cells with no records, or sampling density in the lowest tercile).
3.  Quantifies drifting longline (LLD) apparent fishing effort from Global Fishing Watch (GFW) and classifies high fishing effort areas (upper tercile), alongside a broader fishing footprint (any recorded activity, regardless of intensity).
4.  Computes the spatial overlap between undersampled areas and both the LLD footprint and high-effort areas, globally and by ocean basin.
5.  Classifies overlap cells by jurisdiction (Exclusive Economic Zones vs. High Seas) and identifies the top flag states operating within each, to highlight fleets and governance pathways for adopting opportunistic biodiversity monitoring.

**External data** (not included in this repository; see `R/data_paths.R` for expected paths):

| Dataset                                      | Source                                                    | Used for                         |
|------------------------|------------------------|------------------------|
| GEBCO bathymetry                             | <https://www.gebco.net>                                   | Ocean mask                       |
| Natural Earth land polygons (10m, 50m)       | <https://www.naturalearthdata.com>                        | Ocean mask, map backgrounds      |
| OBIS full occurrence export (Parquet)        | <https://obis.org/data/access/>                           | Sampling density                 |
| Global Fishing Watch apparent fishing effort | <https://globalfishingwatch.org> (API, via `gfwr`)        | Fishing effort                   |
| Global Oceans and Seas (GOaS v1)             | Flanders Marine Institute, <https://doi.org/10.14284/323> | Ocean basin assignment           |
| World EEZ / High Seas boundaries             | Marine Regions, <https://www.marineregions.org>           | EEZ vs. High Seas classification |

**GFW API token**: scripts that query Global Fishing Watch require a personal access token set as the environment variable `GFW_TOKEN` (see <https://globalfishingwatch.org/our-apis/> to request one).

## Running the pipeline

Scripts are numbered and should be run in order within each folder, and folders in the order listed above (`01_study_area` → `02_obis` → `03_fishing_effort` → `04_overlap`). Each script sources `R/data_paths.R` and, where relevant, `R/utils.R`.

## License

MIT

## Citation

If you use this code, please cite:

PENDIENTE
