# Biophysical cartography of the native and human-engineered antibody landscapes quantifies the plasticity of antibody developability

## Overview

Designing effective monoclonal antibody (mAb) therapeutics faces a multi-parameter optimization challenge known as “developability”, which reflects an antibody’s ability to progress through development stages based on its physicochemical properties. While natural antibodies may provide valuable guidance for mAb selection, we lack a comprehensive understanding of natural developability parameter (DP) plasticity (redundancy, predictability, sensitivity) and how the DP landscapes of human-engineered and natural antibodies relate to one another. These gaps hinder fundamental developability profile cartography. To chart natural and engineered DP landscapes, we computed 40 sequence- and 46 structure-based DPs of over two million native and human-engineered single-chain antibody sequences. We found lower redundancy among structure-based compared to sequence-based DPs. Sequence DP sensitivity to single amino acid substitutions varied by antibody region and DP, and structure DP values varied across the conformational ensemble of antibody structures. Sequence DPs were more predictable than structure-based ones across different machine-learning tasks and embeddings, indicating a constrained sequence-based design space. Human-engineered antibodies were localized within the developability and sequence landscapes of natural antibodies, suggesting that human-engineered antibodies explore mere subspaces of the natural one. Our work quantifies the plasticity of antibody developability, providing a fundamental resource for multi-parameter therapeutic mAb design.

<img src="./figures/main/Figure_1_new.png">

## Availability

Datasets are available on GitHub: https://github.com/csi-greifflab/developability_profiling. Structural data, such as predicted models and MD trajectories, are stored on Zenodo194. The raw sequencing data for the human donor is available on the Sequence Read Archive (BioProject number PRJNA1043047). Codes are available on GitHub from: https://github.com/csi-greifflab/developability_profiling. For reference, please cite our paper titled "Biophysical cartography of the native and human-engineered antibody landscapes quantifies the plasticity of antibody developability: [https://www.biorxiv.org/content/10.1101/2023.10.26.563958v2](https://www.nature.com/articles/s42003-024-06561-3). 
## Python environment setup

Recommended **Python version**: **3.9**.

Prepare a virtual environment, using `virtualenv` here (but Conda also serves the purpose):

    pip install virtualenv
    virtualenv adaptyv-dl-light
    source adaptyv-dl-light/bin/activate

Install the requirements:

    pip install -r adaptyv-dl-light.txt

### A note on formats

* `.npy` files are simple Numpy arrays of shape _(samples_num, embedding_dim)_.
* `.npz` files also contain Numpy arrays, but should be read using `numpy.load([name_of_file])['arr_0']`.
