# Differential Privacy via Distributionally Robust Optimization
GitHub supplement of the paper _Differential Privacy via Distributionally Robust Optimization_ (`DPviaDRO`) by Aras Selvi, Huikang Liu, and Wolfram Wiesemann (2023).

The preprint is available on [Optimization Online](https://optimization-online.org/2023/04/differential-privacy-via-distributionally-robust-optimization/) and [arXiv](https://arxiv.org/abs/2304.12681). The paper was presented at the [Robust Optimization Webinars](https://youtu.be/HIfNWrQ-NS4).


## Description
The following is a guide to use this repository.
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

## Final Notes
The following scripts are also available upon request:
- The C++ codes that we supplied are the simple implementation of our cutting-plane algorithms. We have several modifications that were used on the HPC to run parallel jobs.
- We shared how we cleaned datasets for two specific datasets. Similar codes for each of the datasets we used, and the final data itself are available upon request.

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
