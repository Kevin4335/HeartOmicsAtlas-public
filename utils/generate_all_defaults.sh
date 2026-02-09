#!/bin/bash
# generate_all_defaults.sh

cd /home/ubuntu/HeartOmics

echo "Generating ACM_VCM_SAN default plots..."
cd resources-NEW/SAN_ACM_VCM && Rscript ACM_VCM_SAN_default_plots.R && cd ../..

echo "Generating Mini-heart default plots..."
cd resources-NEW/Mini-heart && Rscript mini_heart_default_plots.R && cd ../..

echo "Generating SAN-PCO default plots..."
cd resources-NEW/SAN-PCO && Rscript SAN_PCO_default_plots.R && cd ../..

echo "Generating Spatial default plots..."
cd /home/ubuntu/HeartOmics/resources-NEW/spatial_data && Rscript Spatial_default_plots.R && cd /home/ubuntu/HeartOmics

echo "Generating Multiomics default plots..."
cd /home/ubuntu/HeartOmics/resources-NEW/multi_omics_data && Rscript Multiomics_default_plots.R && cd /home/ubuntu/HeartOmics

echo "All default plots generated!"