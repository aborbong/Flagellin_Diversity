#$ -N runPhylogeny
#$ -M aborbon@tuebingen.mpg.de
#$ -m ea
#$ -l h_vmem=50G
#$ -l h_rt=70:00:00
#$ -cwd

#export PATH=/ebio/abt3_projects/software/miniconda3/bin:/ebio/abt3_projects/small_projects/aborbon/my_interproscan/interproscan-5.44-79.0:/ebio/abt3_projects/software/miniconda3/envs/megagta/bin:/ebio/abt3_projects/software/miniconda3/condabin:/ebio/abt3/aborbon/bin/direnv:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/ebio/abt3/aborbon/bin/go/bin:/ebio/abt3_projects/small_projects/aborbon/edirect:/ebio/abt3/aborbon/psipred/bin:/ebio/abt3_projects/small_projects/aborbon/MaxQuant/bin:/ebio/abt3_projects/small_projects/aborbon/shortbred

source ~/.bashrc
conda activate py3_fla_2022

file=seqs.flagellins.presence.absence.faa

mafft --ep 0.123 --auto $file > mafft_${file%.faa}.fasta

trimal -gappyout -in mafft_${file%.faa}.fasta -out trimal_mafft_${file%.faa}.fasta

fasttree -lg trimal_mafft_${file%.faa}.fasta > fasttree_trimal_mafft_${file%.faa}.fasta
