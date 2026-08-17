# Visual figure verification

The ECFP4 target heatmap renders cleanly at 4500×2400 px. It shows the intended 12×14 compound-by-target matrix, readable annotations, and a continuous 0–0.8 color scale. The strongest visible similarity regions are the Gyr/GyrB columns for BI-1/BI-6 and the Gyr/GyrB/MurC-related columns for parts of the X1V and T2Z series; most values are modest, which supports conservative interpretation.

The TMAP-like reference map renders cleanly at 4200×3000 px. The 12 user compounds are clearly labeled and separated into recognizable chemical series: X1V sulfonamide-benzothiazoles cluster together, BI-1/BI-6 form a neighboring pair, and the T2Z/OX-11 xanthine/heteroaromatic series occupies a separate region. Reference ligands are shown with smaller, semi-transparent markers and target-class legend entries. The caption explicitly states that the map is UMAP on ECFP4/Jaccard distance, avoiding an unsupported claim of exact TMAP implementation.

The final organism-adjusted ranking figure was regenerated after adding a visible target-class legend. It now renders at 3645×3600 px with six organism facets, compound labels, the 0.35 moderate-evidence guide line, and a readable legend. The figure is suitable for report inclusion; target-class color interpretation is now explicit.
