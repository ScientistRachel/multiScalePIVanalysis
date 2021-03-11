# multiScalePIVanalysis
## PIV-based tools for measuring collective cell migration

This code was developed to measure the collective migration of two-dimensional cell sheets that have been imaged over time using phase contrast microscopy.  As described below, _main_batch_script_pairedEdges.m_ will run a suite of analysis tools over a set of images to create multiple .mat and image files of figures as output.  Individual analysis tools can also be run separately given correctly organized inputs.

The current version of the code was adapted for analysis of the data in:
- Lee, Rachel M., Michele I. Vitolo, Wolfgang Losert, and Stuart S. Martin. “Distinct Roles of Tumor-Associated Mutations in Collective Cell Migration.” BioRxiv, June 5, 2020, 2020.06.04.135178. https://doi.org/10.1101/2020.06.04.135178.

Prior iterations of this analysis were also used in the following publications
- Stuelten, Christina H., Rachel M. Lee, Wolfgang Losert, and Carole A. Parent. “Lysophosphatidic Acid Regulates the Motility of MCF10CA1a Breast Cancer Cell Sheets via Two Opposing Signaling Pathways.” Cellular Signalling 45 (May 2018): 1–11. https://doi.org/10.1016/j.cellsig.2018.01.005.
- Lee, Rachel M., Christina H. Stuelten, Carole A. Parent, and Wolfgang Losert. “Collective Cell Migration over Long Time Scales Reveals Distinct Phenotypes.” Convergent Science Physical Oncology 2, no. 2 (May 19, 2016): 025001. https://doi.org/10.1088/2057-1739/2/2/025001.
- Lee, Rachel M., Douglas H. Kelley, Kerstin N. Nordstrom, Nicholas T. Ouellette, and Wolfgang Losert. “Quantifying Stretching and Rearrangement in Epithelial Sheet Migration.” New Journal of Physics 15, no. 2 (2013): 025036. https://doi.org/10.1088/1367-2630/15/2/025036.

## Dependencies
All scripts were developed and tested in MATLAB 2019a with the Image Processing Toolbox.  Some scripts further require additional dependencies:
- To work with images in a format besides .tif (e.g. .mvd2 or .zvi), you will need to download OME's Bio-Format's toolbox for MATLAB.  See the OME website for more details: [https://www.openmicroscopy.org/bio-formats/downloads/](https://www.openmicroscopy.org/bio-formats/downloads/).  This toolbox is required for _convertMVD2.m_.
- The core PIV calculations are performed by MATPIV 1.6.1, which was developed by Johnan Sveen. His work should be cited in any uses of this code: Sveen, Johan Kristian. “An Introduction to MatPIV v. 1.6.1,” 2004. https://www.duo.uio.no/handle/10852/10196.

## Demo
A set of four examples images is available for download along with an Excel file indicating the microscope stage positions of those images (2.49 GB total): [ExampleImages](https://drive.google.com/drive/folders/1LdO4W2dw7qUMEHlDFs0smVRgTt2Wx2rG?usp=sharing).  These images were taken every 3 min for 12 h and have a pixel scale of 0.582 μm. Example output from running _main_batch_script_pairedEdges.m_ over these images is provided in the [ExampleOutput](https://drive.google.com/drive/folders/1y86h-95zrj-62u4K3dSNmjIOwVE8yPjH?usp=sharing) folder (2.15 GB total).  Running on six cores of a PC with a Intel(R) i5-9600K CPU and 16 GB of RAM it takes approximately 48 min to run the full analysis suite over the provided example images.
