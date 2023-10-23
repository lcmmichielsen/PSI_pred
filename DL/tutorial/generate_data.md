# Throughout the tutorial we only show 0.1 and not 0.0
# Models trained on all exons perform better for the DL models, but 0.1 is faster to train

## Data generation

### Replace 0.1 with 0.0 if interested in all exons 
python -u exon_data.py \
--dir /exports/humgen/lmichielsen/Zenodo/Human/HPC \
--file_glia PSI/PSI_glia_norm.csv \
--file_neur PSI/PSI_neur_norm.csv \
--file_fasta /exports/humgen/idenhond/genomes/GRCh38.primary_assembly.genome.fa \
--out /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01 \
--var_threshold 0.1 \
--RBPs "" 

### To add all RBPs as well
python -u exon_data.py \
--dir /exports/humgen/lmichielsen/Zenodo/Human/HPC \
--file_glia PSI/PSI_glia_norm.csv \
--file_neur PSI/PSI_neur_norm.csv \
--file_fasta /exports/humgen/idenhond/genomes/GRCh38.primary_assembly.genome.fa \
--out /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01 \
--var_threshold 0.1 \
--RBPs "AATF,ABCF1,AKAP1,APOBEC3C,AQR,BCCIP,BUD13,CDC40,CPEB4,CPSF6,CSTF2,CSTF2T,DDX21,DDX24,DDX3X,DDX42,DDX51,DDX52,DDX55,DDX6,DGCR8,DHX30,DKC1,DROSHA,EFTUD2,EIF3D,EIF3G,EIF3H,EIF4G2,EWSR1,EXOSC5,FAM120A,FASTKD2,FKBP4,FMR1,FXR1,FXR2,G3BP1,GEMIN5,GNL3,GPKOW,GRSF1,HLTF,HNRNPA1,HNRNPC,HNRNPL,HNRNPM,HNRNPU,HNRNPUL1,IGF2BP1,IGF2BP2,IGF2BP3,ILF3,KHDRBS1,KHSRP,LARP4,LARP7,LIN28B,LSM11,MATR3,METAP2,NCBP2,NIP7,NOL12,NOLC1,NONO,NSUN2,PABPC4,PABPN1,PCBP1,PCBP2,PPIG,PPIL4,PRPF4,PRPF8,PTBP1,PUM1,PUM2,PUS1,QKI,RBM15,RBM22,RBM5,RPS11,SAFB,SAFB2,SBDS,SERBP1,SF3A3,SF3B1,SF3B4,SLBP,SLTM,SMNDC1,SND1,SRSF1,SRSF7,SRSF9,SSB,STAU2,SUGP2,SUPV3L1,TAF15,TARDBP,TBRG4,TIA1,TIAL1,TRA2A,TROVE2,U2AF1,U2AF2,UCHL5,UTP18,UTP3,WDR3,WDR43,XPO5,YBX3,YWHAG,ZC3H11A,ZNF622,ZRANB2"

### To create a mouse dataset
### Use predefined folds to ensure that mouse and human homologs are in the same testset
python -u exon_data.py \
--dir /exports/humgen/lmichielsen/Zenodo/Mouse/HPC \
--file_glia PSI_glia_downsampled.csv \
--file_neur PSI_neur_downsampled.csv \
--file_folds exon_info_folds.csv \
--file_fasta /exports/humgen/idenhond/genomes/mm10.ml.fa \
--out /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_mouse_seqonly_var01 \
--var_threshold 0.1 \
--RBPs "" 