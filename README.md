# SpatialTranscriptomics
Code for troubleshooting, analyzing and displaying SlideSeq data
"Prepare" uses the python version of SCTransform to take UMI read matrices and do dimensionality reduction and making a UMAP.  

"Visualize" takes the UMAP from "Prepare" and combines with the spatial information in the SlideSeq barcode file.  Using plotly and napari one can draw on regions in
either sequencing space or in real space and see where those bead locations lie in the other.  One can also find enriched genes.

"Dash" is similar to "Visualize" except that it runs as a simplified Dash app accessible via the web.

