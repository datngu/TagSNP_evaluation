##### configure.sh: set the current path for loading packages and functions
istr="/path/to"
ostr="$PWD"
cat pipeline_source.sh > TagIt_pipeline.sh
eval "sed -i -e 's#"$istr"#"$ostr"#g' TagIt_pipeline.sh"
echo "Done!"