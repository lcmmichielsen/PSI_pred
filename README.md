# PSI_pred

/athena/tilgnerlab/store/lim4020/HumanBrainHippProject

Downsamples PSI values and embeddings: 
- downsampled_PSI.csv (contains also column 'HumanGene' with the genename and 'variabilityStatus' with the group it belongs to (low/medium/high/cons_low/cons_high/fake))
- downsampled_emb.csv

PSI values:
- variable exons: altInHumans_10_90_cellType_withClassification
- const. exons: consInHumans_5_95_cellType
- fake exons: exons_fake.csv (--> just zeros everywhere)

SpliceAI embeddings: 
- exons_var_emb_start.csv
- exons_var_emb_end.csv
- exons_cons_emb_start.csv
- exons_cons_emb_end.csv
- exons_fake_emb_start.csv
- exons_fake_emb_end.csv


