![GitHub](https://img.shields.io/github/license/ScientistRachel/multiScalePIVanalysis)
# multiScalePIVanalysis
## PIV-based tools for measuring collective cell migration

This code was developed to measure the collective migration of two-dimensional cell sheets that have been imaged over time using phase contrast microscopy.  As described below, _main_batch_script_pairedEdges.m_ will run a suite of analysis tools over a set of images to create multiple .mat and image files of figures as output.  Individual analysis tools can also be run separately given correctly organized inputs.

The current version of the code was adapted for analysis of the data in:
- Lee, Rachel M., Michele I. Vitolo, Wolfgang Losert, and Stuart S. Martin. “Distinct Roles of Tumor-Associated Mutations in Collective Cell Migration.” BioRxiv, June 5, 2020, 2020.06.04.135178. https://doi.org/10.1101/2020.06.04.135178.

Prior iterations of this analysis were also used in the following publications
- Stuelten, Christina H., Rachel M. Lee, Wolfgang Losert, and Carole A. Parent. “Lysophosphatidic Acid Regulates the Motility of MCF10CA1a Breast Cancer Cell Sheets via Two Opposing Signaling Pathways.” Cellular Signalling 45 (May 2018): 1–11. https://doi.org/10.1016/j.cellsig.2018.01.005.
- Lee, Rachel M., Christina H. Stuelten, Carole A. Parent, and Wolfgang Losert. “Collective Cell Migration over Long Time Scales Reveals Distinct Phenotypes.” Convergent Science Physical Oncology 2, no. 2 (May 19, 2016): 025001. https://doi.org/10.1088/2057-1739/2/2/025001.
- Lee, Rachel M., Douglas H. Kelley, Kerstin N. Nordstrom, Nicholas T. Ouellette, and Wolfgang Losert. “Quantifying Stretching and Rearrangement in Epithelial Sheet Migration.” New Journal of Physics 15, no. 2 (2013): 025036. https://doi.org/10.1088/1367-2630/15/2/025036.

## Dependencies and Installation
All scripts were developed and tested in MATLAB 2019a with the Image Processing Toolbox.  Before running the batch analyis script or edge detection for the first time on a system, please run _mexme_dijkstra.m_ to compile the edge detection algorithm (this file is available in the folder EdgeDetection). 

Some scripts further require additional dependencies:
- To work with images in a format besides .tif (e.g. .mvd2 or .zvi), you will need to download OME's Bio-Format's toolbox for MATLAB.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).  This toolbox is required for _convertMVD2.m_.
- The core PIV calculations are performed by MatPIV 1.6.1, which was developed by Johnan Sveen. His work should be cited in any uses of this code: Sveen, Johan Kristian. “An Introduction to MatPIV v. 1.6.1,” 2004. https://www.duo.uio.no/handle/10852/10196.

## Demo

### Example Images and Output
A set of four examples images is available for download along with an Excel file indicating the microscope stage positions of those images (2.49 GB total): [ExampleImages](https://drive.google.com/drive/folders/1LdO4W2dw7qUMEHlDFs0smVRgTt2Wx2rG?usp=sharing).  These images were taken every 3 min for 12 h and have a pixel scale of 0.582 μm. Example output from running _main_batch_script_pairedEdges.m_ over these images is provided in the [ExampleOutput](https://drive.google.com/drive/folders/1y86h-95zrj-62u4K3dSNmjIOwVE8yPjH?usp=sharing) folder (2.15 GB total).  Running on six cores of a PC with a Intel(R) i5-9600K CPU and 16 GB of RAM, it takes approximately 48 min to run the full analysis suite over the provided example images.

### main_batch_script_pairedEdges
This script will run the full suite of tools on a set of .tif images.  The images are assumed to come in pairs; this script was designed to analyze data where two regions of the same circular cell sheet were imaged over time.  The pairs will be combined for downstream analysis into the same radial coordinate system using microscope stage positions (provided by the user as an Excel sheet or extracted from the microscopy data -- see for example _convertMVD2.m_ for extraction of microscope stage positions from a .mvd2 file format).

#### User Inputs
- `directories`: A cell array of full path directories containing multi-image .tif files for analysis.
- `over_write`: When set to 0, previously analyzed images will be skipped during downstream analysis. When set to 1, those files will be *replaced* with new analysis.
- `overwriteEdge`: When set to 0, previously analyzed images will be skipped during edge detection. When set to 1, those files will be *replaced* with new analysis.
- `user_firstframe` and `user_lastframe`: These set the frames for downstream analysis (can be used to only analyze the first 10 h, for example).
- `r_scale`: The resolution of the images in microns per pixel.
- `t_scale`: The time scale of image acquisition, in minutes per frame
- Filtering parameters: `size_exclude`, `size erode`, and `size_gap` are used in the function _mask_and_filter.m_, which removes non-cell areas from the PIV analysis.  These parameters are robust for images collected at the same resolution and do not regularly need to be changed; if filtering of non-cell areas is not workign well, read the header of _mask_and_filter.m_ for more information.
- `T` and `trange`: These parameters set the deformation time (T) and range of frames over which to calculate the finite time Lyapunov exponents (FTLEs).  T should be set large enough that FTLE values have reached an asymptote; practically for epithelial cell migration this is approximately 2 h.  The trange should maximally extend to user_lastframe-T.
- `pos_type`: When set to 1, downstream analysis expects to find microscope stage positions in an Excel sheet named PosList.  With the OME Bioformats package, setting this parameter to 0 will automatically extract positions from .zvi files.
- `parpoolName`: This parameter allows the user to set a custom parallel profile in  MATLAB.  Set to 'local' to use the default MATLAB profile.
- `frameskip`: During edge detection and filtering checks, frames will be plotted as a check on the analysis.  Often plotting every frame (frameskip = 1) is memory-intensive and unnecessary, so values of frameskip = 25 or 50 are recommended.

#### Batch Processing Structure
1. Images are first loaded into _dot_matpiv5_, a function that applies a multipass PIV analysis to the entire .tif file.  Two iterations of 64 x 64 windows are followed by two iterations of 32 x 32 windows, using 50% overlap.  This is the step of the analysis that relies on the MatPIV toolbox by Johnan Sveen.  This creates a \*dot_piv.mat file containing the PIV flow fields (dot_piv.u_fi and dot_piv.v_fi) as well as quality metrics (dot_piv.snr, dot_piv.pkh).
2. Images are then filtered to remove non-cell areas (if the images are entirely covered by cells, this step should be skipped).  This step uses sobel filtering and morphological operators from the Image Processing Toolbox. This creates a \*dot_piv_filtered.mat file with the same fields as \*dot_piv.mat.  The filtered file is used in all downstream analysis.  At this stage, a quiver visualization of the flow field is created using _plot_piv_quiverSpeed_v2.m_ and saved in the quiverSpeed folder.
3. The leading edge of each cell sheet is then detected using the dijkstra algorithm.  Note that before running this script for the first time, the user should run _mexme_dijkstra.m_ (file available in the EdgeDetection folder) to compile the edge detection algorithm. This step creates a \*edgedat.mat file containing the coordinates of the segmented edge, as well as image overlaying edge detection on the original images (found in the EdgeImposed folder).
4. Next edge 'ramen noodle' plots are created, so named because of their wavy appearance.  These images represent the progression of the segmented leading edge over time.  These files are saved to a LengthStacks folder, with the user_lastframe added to the folder name.
5. Finite-time Lyapunov exponent (FTLE) values are calculated for each PIV flow field, which reflect the exponential sensitivity of the flow field to initial conditions.  These are saved in FTLEall.mat.
6. At this point, the user is asked to watch a series of 'filter check' images that confirm that the PIV filtering and edge detection align well.  The edge detection is calculated on a pixel-by-pixel basis, while PIV is calculated on a 16 x 16 pixel grid, so these two images will never perfectly align, but the goal is to see rough alignment with the two types of segementation.  These images are saved in the FilterCheck folder.
7. Downstream analysis of the PIV flow fields begins with a basic set of analysis that measures speed and directionality, which are saved in the Analytics Data and Analytics Figures folders.
8. This is followed by an analysis of these values over time, as well as histograms of each value over time (saved in the Analytics Figures and PIV Hist folders).
9. Additionally, plots of migration metrics versus location in the cell sheet are created (saved in the PIV Time and R folder).
10. Downstream analysis of the FTLE values proceeds similarly: First basic averages are calculated, followed by histograms, time, and spatial analysis (saved in the FTLE Analytics and FTLE Time and R folders).
11. A coarse graining approach is used to calculate characteristic length and time scales for the flow field (saved in the CoarseGrain folder).
12. Spatial auto-correlations are calculated for several metrics of interest (saved in the Correlation Data folder).

Each step of the batch processing script is accomplished through subfunctions which are organized into subfolders in this repository.  In general, the function called by main_batch_script_pairedEdges parses inputs and images, before then passing the information to a further downstream analysis function.
