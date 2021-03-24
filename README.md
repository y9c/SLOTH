# SLOTH (Single-cell Lineage on Targeted Hypermutation)

<p align="center">
  <img src="docs/SLOTH.gif?raw=true" alt="SLOTH" title="SLOTH" width="200" height="200">
</p>

---

[![DOI:TBD/TBD](https://zenodo.org/badge/DOI/TBD.svg)](https://doi.org/TBD/TBD)
[![GitHub](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE.md)

Mapping single-cell-resolution cell phylogeny reveals cell population dynamics during organ developments

## PROCEDURE

- **Single Molecule Sequecing Pipeline** (smrt)

  - prepare pacbio data
  - quality stat
  - split read
  - unify sequence
  - call mutation

- **Single Cell Tree Construction** (tree)

  - build tree by maximum likelihood method
  - visualize tree
  - benchmarking cell lineage tracing system (tree/benchmark)

- **Population Dynamics** (pop)

  - CoalescentSimulation
  - NpDynamic
  - SShape
  - TissueSimilarity

- simulation (simu)

## READ MORE

scientific question :
https://celllineage.github.io/SLOTH/index.html

technical details:
https://github.com/CellLineage/SLOTH/wiki

---

## TODO

- [ ] turn consensus part into a independent package
- [ ] write wiki page for the technical details
- [ ] write web page for the scientific question
- [ ] add js framework for renderng the tree

## LICENSE

The content of this project itself is created by **Ye Chang** and licensed under the [Creative Commons Attribution 4.0 International License (CC BY)](https://creativecommons.org/licenses/by/4.0/),
and the underlying source code used to format and display that content is licensed under the [MIT license](LICENSE.md).

[![Creative Commons Attribution 4.0 International License](https://github.com/creativecommons/cc-cert-core/blob/master/images/cc-by-88x31.png 'CC BY')](https://creativecommons.org/licenses/by/4.0/)
