#!/bin/bash

# Define environment name
ENV_NAME="paper_fig_env"

# 1. Create environment if it doesn't exist (update if it does)
echo "Ensuring conda environment '$ENV_NAME' is installed..."
conda env update -f environment.yaml --prune

# 2. Verify installation using 'conda run' (avoids shell activation issues)
echo "----------------------------------------------------------------"
echo "Verifying R packages in '$ENV_NAME'..."
echo "----------------------------------------------------------------"

conda run -n $ENV_NAME Rscript -e "
pkgs <- c(
  'ggplot2', 'igraph', 'Hmisc', 'ggraph', 'RColorBrewer', 'scales',
  'circlize', 'vegan', 'data.table', 'mixOmics', 'ComplexHeatmap',
  'pheatmap', 'ggrepel'
)

message('Checking package status...')
installed <- sapply(pkgs, requireNamespace, quietly = TRUE)
if (all(installed)) {
  message('[OK] All required packages are found!')
} else {
  missing <- pkgs[!installed]
  stop('Missing packages: ', paste(missing, collapse = ', '))
}
"

# 3. Final instructions
if [ $? -eq 0 ]; then
    echo "----------------------------------------------------------------"
    echo "Setup Successful!"
    echo "To start using the environment, run:"
    echo "    conda activate $ENV_NAME"
    echo "----------------------------------------------------------------"
else
    echo "----------------------------------------------------------------"
    echo "Verification Failed. Checks above for errors."
    echo "----------------------------------------------------------------"
fi
