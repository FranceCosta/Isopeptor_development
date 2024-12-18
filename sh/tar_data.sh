#!/bin/bash -eu

# Store data in a .tar.gz format to be saved on Zenodo
output_dir="tmp"
current_date=$(date +"%Y%m%d")

biochem_analysis="/nfs/research/agb/research/francesco/projects/20241213_IsopeptorDevelopment_v2/analysis/20240904_biochemAnalysis_v1/output" 
test_templates="/nfs/research/agb/research/francesco/projects/20241213_IsopeptorDevelopment_v2/analysis/20241126_testTemplates_v1/output"
model_development="/nfs/research/agb/research/francesco/projects/20241213_IsopeptorDevelopment_v2/analysis/20241203_model_v3/output"

mkdir -p "$output_dir"

# Archive each directory
declare -A dirs=(
    ["biochem_analysis"]=$biochem_analysis
    ["test_templates"]=$test_templates
    ["model_development"]=$model_development
)

for name in "${!dirs[@]}"; do
    dir="${dirs[$name]}"
    if [[ -d "$dir" ]]; then
        tar -czf "$output_dir/${name}.tar.gz" -C "$(dirname "$dir")" "$(basename "$dir")"
        echo "Archived $dir to $output_dir/${name}.tar.gz"
    else
        echo "Warning: Directory $dir does not exist, skipping."
    fi
done

# Create the final archive containing all individual archives
final_archive="$output_dir/${current_date}_data.tar.gz"
(
    cd "$output_dir" && tar --exclude="${current_date}_data.tar.gz" -czf "${current_date}_data.tar.gz" *.tar.gz
)

echo "Final archive created at $final_archive"

