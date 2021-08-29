##### configure.sh: set the current path for loading packages and functions
istr="/path/to"
ostr="$PWD"
cat pipeline_source.sh > FastTagger_pipeline.sh
eval "sed -i -e 's#"$istr"#"$ostr"#g' FastTagger_pipeline.sh"
echo "Done!"