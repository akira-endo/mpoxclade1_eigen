# mpoxclade1_eigen
Code accompanying [Murayama & Asakura et al. "Roles of community and sexual contacts as drivers of clade I mpox outbreaks"](https://www.medrxiv.org/content/10.1101/2024.10.15.24315554)

## Instructions
* Codes are reproducible via [Codespaces with Jupyter](https://github.blog/changelog/2022-11-09-using-codespaces-with-jupyterlab-public-beta/) or [Docker Desktop](https://www.docker.com/products/docker-desktop/) through Compose plugin. Docker build may take up to about 30 mins.
* Use [Jupytext](https://jupytext.readthedocs.io/en/latest/) to generate Jupyter Notebooks from the source files in /notebook. Main analysis may take up to about 30 mins.
* If Jupyter produces a kernel error, launch Julia from the terminal and try `]activate .` and `]instantiate` before reopening the notebook.
* See ./devcontainer/Dockerfile and Manifest.toml for environments and software dependencies
