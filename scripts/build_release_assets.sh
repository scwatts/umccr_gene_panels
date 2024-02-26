#!/usr/bin/env bash

release_version=24.02.0


# Set up build directory
build_dp=build/${release_version}

[[ -e ${build_dp}/ ]] && rm -r ${build_dp}/
mkdir -p ${build_dp}/


# Build panel data tarball
data_dn=panel-data-v${release_version}
data_dp=${build_dp}/${data_dn}

mkdir -p ${data_dp}/

rsync -aP somatic_panel/4_panel_data/output/ ${data_dp}/somatic/
rsync -aP germline_panel/3_panel_data/output/ ${data_dp}/germline/
rsync -aP \
  --exclude='fusion_database.tsv' \
  --exclude='hmftools_fusion_data.subset.csv' \
  fusion_panel/1_panel_data/output/ ${data_dp}/fusion/

pushd ${build_dp}/
tar -czvf ${data_dn}.tar.gz ${data_dn}/
popd

rm -r ${data_dp}/


# Panels
# NOTE(SW): fusion panel is currently the somatic panel, so is skipped
cp somatic_panel/3_final_panel/final_panel.tsv ${build_dp}/somatic_panel-v${release_version}.tsv
cp germline_panel/2_final_panel/final_panel.tsv ${build_dp}/germline_panel-v${release_version}.tsv
