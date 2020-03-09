### Growth factor-mediated coupling between lineage size and cell fate choice underlies robustness of mammalian development

#### N. Saiz, L. Mora-Bitria, S. Rahman, H. George, J.P. Herder, J. García-Ojalvo and A.-K. Hadjantonakis

This repository contains the code to reproduce the modeling figures of the paper by Saiz et al, posted in bioRxiv (preprint [2019.12.27.889006](http://dx.doi.org/10.1101/2019.12.27.889006)). It includes the following:

* [`fig2b.ipynb`](fig2b.ipynb): self-contained Jupyter notebook to generate the phase portrait shown in Figure 2b.
* [`fig2c.ipynb`](fig2c.ipynb): Jupyter notebook to generate the lineage dynamics plot shown in Figure 2c. This notebook compiles and runs the C program `embryo_v1.c`, and plots the correspoding results. It requires the C compiler `cc` and the utility `make`, which are available by default in Linux, or by installing the `Xcode` development environment in Mac OS X.
* [`fig2e.ipynb`](fig2e.ipynb): Jupyter notebook to generate the lineage dynamics plots shown in Figure 2e, using again the C program `embryo_v1.c`.
* [`fig3.ipynb`](fig3.ipynb): Jupyter notebook to generate the ESC addition plots shown in Figure 3, using in this case the C program `embryo_v2.c`.
* [`fig4.ipynb`](fig4.ipynb): Jupyter notebook to generate the ablation plots shown in Figure 4, using again the C program `embryo_v1.c`.
