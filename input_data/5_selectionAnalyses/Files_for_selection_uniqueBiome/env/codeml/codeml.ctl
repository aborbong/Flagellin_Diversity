      seqfile = ../paml.env.phy            * Path to the alignment file
     treefile = ../tree.env_woBL.nwk           * Path to the tree file
      outfile = out_sites_env_codonFreq0.txt            * Path to the output file
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = 1           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
        model = 0         * Models for ω varying across lineages
	  NSsites = 0 1 2 7 8          * Models for ω varying across sites
    CodonFreq = 7        * Codon frequencies
	  estFreq = 1        * Use observed freqs or estimate freqs by ML
        clock = 0          * Clock model
    fix_omega = 0         * Estimate or fix omega
        omega = 0.5        * Initial or fixed omega

