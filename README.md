[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6572895.svg)](https://doi.org/10.5281/zenodo.6572895)

# Source code for research article "Audiovisual adaptation is expressed in spatial and decisional codes" 
**by Máté Aller, Agoston Mihalik, and Uta Noppeney**

Published in [Nature Communications](https://www.nature.com/articles/s41467-022-31549-0)

Cite: Aller, M., Mihalik, A., & Noppeney, U. (2022). Audiovisual adaptation is expressed in spatial and decisional codes. Nature Communications, 13(1), 3924. https://doi.org/10.1038/s41467-022-31549-0

## Installation
1. Clone or download this repository
2. Download all dependencies into their individual subfolders within audiovisual_adaptation_fMRI_EEG/code/toolbox/
3. Download source data from [here](https://doi.org/10.6084/m9.figshare.19469861.v2) and unpack into audiovisual_adaptation_fMRI_EEG/results/data/

## Usage
1. Navigate to the project root folder (audiovisual_adaptation_fMRI_EEG/) and execute set_path.m. This takes care of setting up the environment. 
2. Execute the scripts found in results/ to reproduce the reported results. 
3. All analysis code can be found in audiovisual_adaptation_fMRI_EEG/code/analysis/

## Dependencies
### Standalone toolboxes
- [fieldtrip](https://www.fieldtriptoolbox.org/download/)
- [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
- [Palamedes](https://www.palamedestoolbox.org/download.html)
- [pcm_toolbox](https://github.com/jdiedrichsen/pcm_toolbox)
- [rsatoolbox_matlab](https://github.com/rsagroup/rsatoolbox_matlab)
- [spm12](https://www.fil.ion.ucl.ac.uk/spm/software/download/)
### Matlab File Exchange functions
- [breakxaxis](https://uk.mathworks.com/matlabcentral/fileexchange/42905-break-x-axis)
- [breakyaxis](https://uk.mathworks.com/matlabcentral/fileexchange/45760-break-y-axis)
- [Colormaps](https://uk.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
- [fdr_bh](https://uk.mathworks.com/matlabcentral/fileexchange/29274-dmgroppe-mass_univariate_erp_toolbox)
- [padcat](https://uk.mathworks.com/matlabcentral/fileexchange/22909-padcat)
- [parfor_progress](https://uk.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor-progress-bar-that-works-with-parfor)
- [shadedErrorBar](https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)
- [suplabel](https://uk.mathworks.com/matlabcentral/fileexchange/7772-suplabel)

## Contributors
The code was mainly written by Máté Aller and Agoston Mihalik. Some of the code in the psychometric and neurometric function fitting was kindly provided by David Meijer. Some of the code in the fMRI analysis have been adapted from The Decoding Toolbox (Hebart, Goergen, & Haynes, 2015).
