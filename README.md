# SLOTH (Single-cell Lineage on Targeted Hypermutation)

<p align="center">
  <img src="docs/SLOTH.gif?raw=true" alt="SLOTH" title="SLOTH" width="200" height="200">
</p>

---

[![DOI:TBD/TBD](https://zenodo.org/badge/DOI/TBD.svg)](https://doi.org/TBD/TBD)
[![GitHub](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE.md)

Mapping single-cell-resolution cell phylogeny reveals cell population dynamics during organ developments

## PROCEDURE

- **Single Molecule Sequecing Pipeline** ([smrt](./smrt))

  - Adaptors Trimming
  - Circular Consensus Sequence Generation
  - Align to Reference Sequence
  - Annotate Sequence Features
  - Link Sample Barcode Sequence
  - Cluster Molecular Indentifier (UMI) Sequence
  - Group Reads by Sample Barcode and Molecular Indentifier
  - Generate Consensus Sequence for Lineage Barcode
  - Call Mutation
  - Post Quality Filter

- **Single Cell Tree Construction** ([tree](./tree))

  - Build Tree by Maximum Likelihood Method
  - Fix Tree
  - Visualize Tree

- **Benchmarking Cell Lineage Tracing System** ([benchmark](./benchmark))

  - `alemany_whole-organism_2018`
  - `bowling_engineered_2020`
  - `chan_molecular_2019`
  - `chen_efficient_2020`
  - `kalhor_developmental_2018`
  - `lee-six_population_2018`
  - `mckenna_whole-organism_2016`
  - `pei_polylox_2017`
  - `quinn_single-cell_2021`
  - `raj_simultaneous_2018`
  - `spanjaard_simultaneous_2018`

- **Simulation** ([simu](./simu))

  - ?

- **Population Dynamics** ([pop](./pop))

  - ?

## READ MORE

scientific question :
https://celllineage.github.io/SLOTH/index.html

technical details:
https://github.com/CellLineage/SLOTH/wiki

## How to use?

```bash
git clone https://github.com/CellLineage/SLOTH.git
```

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
