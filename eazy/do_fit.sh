#!/bin/bash

CATALOGS=(../data/COSMOS_data_PH_27_11_18.csv ../data/XMM_data_PH_27_11_18.csv)
BANDS=(72 314 315 316 317 318 256 257 258 259 18 19)
EAZY=/home/cschreib/programming/eazy-photoz/bin/eazy

# Construct header
HEADER="# id RA Dec "
for BAND in ${BANDS[@]}; do
    HEADER="${HEADER} F${BAND}"
done
for BAND in ${BANDS[@]}; do
    HEADER="${HEADER} E${BAND}"
done
HEADER="${HEADER} z_spec"

for CATALOG in ${CATALOGS[@]}; do
    echo Preparing fit for ${CATALOG}

    # Make input catalog in format expected by EAzY
    #  - add a header at the start of the catalog to identify the columns
    #  - add source ID at beginning of each line
    #  - convert CSV to space-separated columns
    echo ${HEADER} > eazy.cat
    I=1
    cat ${CATALOG} | sed "s/,/ /g" | while read line; do
        echo "${I} ${line}" >> eazy.cat
        I=$((I+1))
    done

    echo Doing fit for ${CATALOG}

    # Clear output directory
    rm -rf output
    mkdir -p output

    # Run EAzY
    ${EAZY} -p eazy.param > eazy.log

    echo Saving fit for ${CATALOG}

    # Copy outputs
    mkdir -p output/pz
    mv output/*.pz output/pz
    mkdir -p output/seds
    mv output/*.temp_sed output/*.obs_sed output/seds

    OUT_DIR=../results/eazy_$(basename ${CATALOG} .csv)
    rm -rf ${OUT_DIR}
    mv output ${OUT_DIR}
done
