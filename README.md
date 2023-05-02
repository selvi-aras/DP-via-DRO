# Differential Privacy via Distributionally Robust Optimization
GitHub supplement of the paper _Differential Privacy via Distributionally Robust Optimization_ (`DPviaDRO`) by Aras Selvi, Huikang Liu, and Wolfram Wiesemann (2023).

The preprint is available on [Optimization Online](https://optimization-online.org/2023/04/differential-privacy-via-distributionally-robust-optimization/) and [arXiv](https://arxiv.org/abs/2304.12681). The paper was presented at the [Robust Optimization Webinars](https://youtu.be/HIfNWrQ-NS4).


## Description
The following is a guide to use this repository. The folder "GitHub Appendices" includes the additional experiments and discussions that are supplied with the GitHub repository. We now explain the codes in the following. 
<details>
  <summary> <b> Data Independent Noise </b> </summary>
  
  This folder includes the C++ codes for the data independent noise optimization.
</details>

<details>
<summary> <b> Data Dependent Noise </b> </summary>
  
  This folder includes the C++ codes for the data dependent noise optimization.

</details>

<details>
<summary> <b> Dataset Cook </b> </summary>
  
  This folder has two Python3 Jupyter Notebooks which shows how we cleaned datasets examples: one for naïve Bayes, and another one for proximal coordinate descent.

</details>

<details>
<summary> <b> Visualization </b> </summary>
  
  This folder has examples of Python3 visualizations that were used to generate Figures 3 and 4.

</details>

<details>
<summary> <b> HPC Codes </b> </summary>
  
  This folder includes C++ compilation commands on the Imperial HPC's Linux Terminal as well as an example .PBS jobs file. 

</details>

<details>
<summary> <b> ML Algorithms Julia </b> </summary>
  
  Implementations of ML algorithms (NB, PCD) in Julia. `analytic_gaussian.jl` has the Analytic Gaussian mechanism's Julia implementation (we adopted the original Python implementation in Julia). `compare.jl` has the standard logistic regression (LR) implementation in case one wants to benchmark the non-private optimal LR classifier (uses MOSEK solver). `grid_ell1_prep_main.jl`, `grid_ell1_prep_ub.jl`, `grid_ell1_prep_lb.jl` are the files generating the code for Tables 1,4,5. `interpret_results.jl` includes the analysis (p-values, in/out-sample comparison, etc.) of the DP naïve Bayes methods. `main.jl` is the main DP-NB function including the Gaussian, truncated Laplace, and our optimised data independent mechanism, `main_analytic.jl` is the same for analytic Gaussian (added later on, hence in another function), and `main_data_dependent.jl` is the data dependent noise mechanism version which also prepares files for noise distributions to be optimized in C++ with the help of the function in `read_distributions.jl`. `NB.jl` has helper functions for the naïve Bayes method, including the non-noisy counts/statistics. The noisy versions are computed in the file `sensitivites.jl`. The file `PCD_functions.jl` includes the PCD iteration functions to be called in the iterations of the DP PCD  algorithm. `PCD_histogram.jl` gives the histogram in the GitHub Additional Appendices. `PCD_interpret_results ... .jl` are the files we used to interpret PCD results. `PCD.jl` is the main function that runs the DP PCD algorithm, and PCD_multi.jl is the same which uses the data dependent noise (that was prepared in C++). `sample.jl` includes the functions that allow us to sample noise from various distributions. `train_test.jl` splits datasets into training and test-sets. The files `visualize_bounds ... .jl` were used to generate Figure 2.
</details>


## Final Notes
The following scripts are also available upon request:
- The C++ codes that we supplied are the simple implementation of our cutting-plane algorithms. We have several modifications that were used on the HPC to run parallel jobs.
- We shared how we cleaned datasets for two specific datasets. Similar codes for each of the datasets we used, and the final data itself are available upon request.
- We provided some examples of HPC commands; however, for a specific PBS file, or a specific algorithm's parallelised implementation (to be able to run on the cluster computers), please get in touch with us.

## Thank You
Thank you for your interest in our work. If you found this work useful in your research and/or applications, please star this repository and cite:
```
@article{DPviaDRO,
  title={Differential Privacy via Distributionally Robust Optimization},
  author={Selvi, Aras and Liu, Huikang and Wiesemann, Wolfram},
  journal={arXiv preprint 2304.12681},
  year={2023}
}
```
Please contact Aras (a.selvi19@imperial.ac.uk) if you encounter any issues using the scripts. For any other comment or question, please do not hesitate to contact us:

[Aras Selvi](https://www.arasselvi.com/) _(a.selvi19@imperial.ac.uk)_

[Huikang Liu](https://huikang2019.github.io/) _(liuhuikang@shufe.edu.cn)_

[Wolfram Wiesemann](http://wp.doc.ic.ac.uk/wwiesema/) _(ww@imperial.ac.uk)_
