input_tree=tree.animals.tre.newick #fastTree_trimAI_animals_woBL_25.tre
output_tree=${input_tree%.tre.newick}_woBL.nwk

perl -ne '$_=~s/:[\d\.]+//g; print $_;' $input_tree > $output_tree
