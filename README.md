# wildfire-speciation

Repo supporting [Krasovich Southworth et al. 2025 "The Influence of Wildfire Smoke on Ambient PM2.5 Chemical Species Concentrations in the Contiguous US"](https://doi.org/10.1021/acs.est.4c09011) in *Environmental Science & Technology*.

## Abstract

Wildfires significantly contribute to ambient air pollution, yet our understanding of how wildfire smoke influences specific chemicals and their resulting concentration in smoke remains incomplete. We combine 15 years of daily species-specific PM2.5 concentrations from 700 air pollution monitors with satellite-derived ambient wildfire smoke PM2.5, and use a panel regression to estimate wildfire smoke's contribution to the concentrations of 27 different chemical species in PM2.5. Wildfire smoke drives detectable increases in the concentration of 25 out of the 27 species with the largest increases observed for organic carbon, elemental carbon, and potassium. We find that smoke originating from wildfires that burned structures had higher concentrations of copper, lead, zinc, and nickel relative to smoke from fires that did not burn structures. Wildfire smoke is responsible for an increasing share of ambient concentrations of multiple species, some of which are particularly harmful to health. Using a risk assessment approach, we find that wildfire-induced enhancement of carcinogenic species concentrations could cause increases in population cancer risk, but these increases are very small relative to other environmental risks. We demonstrate how combining ground-monitored and satellite-derived data can be used to measure wildfire smoke's influence on chemical concentrations and estimate population exposures at large scales.

## Organization

Results from the paper (all main text and supplementary figures) are in the Dropbox folder. Code to replicate results is in the `scripts` folder. Data are in [Dropbox](https://www.dropbox.com/scl/fo/gm5zt8k3mif88awv55tvj/AA8M9YR-UVURtvI_DHIouio?rlkey=11pyi7wsq8d6npqq9zatw4vrd&st=zljm63oe&dl=0).

```
wildfire-speciation
├── scripts/
│   ├── A00_master_workflow.R   # Master script — runs all analysis scripts
│   ├── 0_packages.R            # Required R packages
│   └── ...                     # Individual analysis scripts
├── LICENSE
└── README.md
```

## How to replicate results

1. Download this repository.
2. Download data from [Dropbox](https://www.dropbox.com/scl/fo/gm5zt8k3mif88awv55tvj/AA8M9YR-UVURtvI_DHIouio?rlkey=11pyi7wsq8d6npqq9zatw4vrd&st=zljm63oe&dl=0).
3. Check `scripts/0_packages.R` and install any packages you do not already have. Libraries are loaded automatically when running the workflow.
4. Set working directory to this downloaded repository's root.
5. Run `scripts/A00_master_workflow.R` to execute all analysis scripts.

## Data

All data used in this project as well as all figures in the paper can be found in [Dropbox](https://www.dropbox.com/scl/fo/gm5zt8k3mif88awv55tvj/AA8M9YR-UVURtvI_DHIouio?rlkey=11pyi7wsq8d6npqq9zatw4vrd&st=zljm63oe&dl=0).

Key data sources include:

- **EPA Chemical Speciation Network (CSN)** and **Interagency Monitoring of Protected Visual Environments (IMPROVE)**: Daily species-specific PM2.5 concentrations from ~700 monitors across the contiguous US.
- **Satellite-derived wildfire smoke PM2.5**: Daily ambient wildfire smoke PM2.5 estimates from [Childs et al. 2022](https://doi.org/10.1021/acs.est.2c02934) ([repo](https://github.com/echolab-stanford/daily-10km-smokePM)).

## Citation

If you use this code or data, please cite:

> Krasovich Southworth, E., Qiu, M., Gould, C. F., Kawano, A., Wen, J., Heft-Neal, S., Kilpatrick Voss, K., Lopez, A., Fendorf, S., Burney, J. A., & Burke, M. (2025). The Influence of Wildfire Smoke on Ambient PM2.5 Chemical Species Concentrations in the Contiguous US. *Environmental Science & Technology*, 59(6), 2961–2973. https://doi.org/10.1021/acs.est.4c09011

BibTeX:

```bibtex
@article{krasovichsouthworth2025wildfire,
  title     = {The Influence of Wildfire Smoke on Ambient {PM2.5} Chemical Species Concentrations in the Contiguous {US}},
  author    = {Krasovich Southworth, Emma and Qiu, Minghao and Gould, Carlos F. and Kawano, Ayako and Wen, Jeff and Heft-Neal, Sam and Kilpatrick Voss, Kara and Lopez, Alandra and Fendorf, Scott and Burney, Jennifer Anne and Burke, Marshall},
  journal   = {Environmental Science \& Technology},
  volume    = {59},
  number    = {6},
  pages     = {2961--2973},
  year      = {2025},
  doi       = {10.1021/acs.est.4c09011}
}
```

## Contact

If you find meaningful errors or have questions or suggestions, please contact Emma Krasovich Southworth at [emmars@stanford.edu](mailto:emmars@stanford.edu).

## License

This project is licensed under the [MIT License](LICENSE).
