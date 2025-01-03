# Mandai storm forest plots

This is a repository for the data and scripts behind the Mandai storm forest plots described in:

> Yee, A.T.K., H.R. Lai, K.Y. Chong, L. Neo, C.Y. Koh, S.Y. Tan, W.W. Seah, J.W. Loh, R.C.J. Lim, M. van Breugel, H.T.W. Tan. (2019) Short-term responses in a secondary tropical forest after a severe windstorm event. Journal of Vegetation Science 30: [720–731](https://onlinelibrary.wiley.com/doi/abs/10.1111/jvs.12753).

A blog [post](https://vegsciblog.org/2020/05/11/trudging-through-treefalls/) was contributed to the Vegetation Science Blog for this publication.

And also:

> Lai, H.R., K.Y. Chong, A.T.K. Yee. (2022) Ten years after: what we learned from the Mandai storm forest. Nature in Singapore Supplement 1: [e2022121](https://www.science.nus.edu.sg/wp-content/uploads/sites/11/2024/03/NIS_S1_207-218.pdf).

Other publications arising from this dataset include:

- Lai, H.R., K.Y. Chong, A.T.K. Yee, H.T.W. Tan, M. van Breugel. (2020) Functional traits that moderate tropical plant recruitment during post-windstorm secondary succession. Journal of Ecology 108: [1322–1333](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.13347). 
    - Also with a blog [post](https://jecologyblog.com/2020/07/01/functional-traits-that-moderate-tropical-tree-recruitment-during-post%e2%80%90windstorm-secondary-succession/) on the Journal of Ecology blog.
    - Model and processed data can be found in another [GitHub repo](https://github.com/hrlai/Lai_et_al_2020_JEcol_Trait_Env)

- Lai, H.R., K.Y. Chong, A.T.K. Yee, M.M. Mayfield, D.B. Stouffer (2022) The role of higher-order biotic interactions on tropical tree growth. Ecology 103: [e03588](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.3588).
    - Model and processed data can be found in another [GitHub repo](https://github.com/stoufferlab/hoi-trees-public)

# To install

```
devtools::install_github("kwekings/Mandai")
```

The object `Mandai_data` containing the data can then be called (i.e., either after `library(Mandai)` or as `Mandai::Mandai_data`).
