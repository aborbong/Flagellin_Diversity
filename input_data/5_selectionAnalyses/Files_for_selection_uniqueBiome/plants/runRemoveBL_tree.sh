input_tree=plants.tree.tre.newick #fastTree_trimAI_animals_woBL_25.tre
output_tree=${input_tree%.tree.tre.newick}_woBL.nwk

perl -ne '$_=~s/:[\d\.]+//g; print $_;' $input_tree > $output_tree
