# See here for image contents: https://hub.docker.com/r/jupyter/datascience-notebook/

# Specify docker image from https://hub.docker.com/r/jupyter/datascience-notebook/tags.
# Avoid :latest tag. Tags that includes both linux/amd64 and linux/arm64/v8 are recommended.
FROM jupyter/datascience-notebook:2023-10-20

# We want to run common-debian.sh from here:
# https://github.com/microsoft/vscode-dev-containers/tree/main/script-library#development-container-scripts
# But that script assumes that the main non-root user (in this case jovyan)
# is in a group with the same name (in this case jovyan).  So we must first make that so.
COPY .devcontainer/library-scripts/common-debian.sh /tmp/library-scripts/
USER root
RUN apt-get update \
 && groupadd jovyan \
 && usermod -g jovyan -a -G users jovyan \
 && bash /tmp/library-scripts/common-debian.sh \
 && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts

# Jupyter lab extensions
RUN conda install -c conda-forge jupyterlab-git -y && \
    conda install -c conda-forge jupytext -y && \
    jupyter lab build

# Jupyter lab settings
# Alt + Space to send the highlighted code to console
# auto-close brackets + line numbers in notebooks and console
# turn off autosave
RUN printf ' { \
    "@jupyterlab/shortcuts-extension:shortcuts": {   \ 
        "shortcuts":[{"command": "notebook:run-in-console", "keys": ["Alt Space"], "selector": ".jp-Notebook.jp-mod-editMode"}, \
                     {"command": "runmenu:run", "keys": ["Alt Space"], "selector": "[data-jp-code-runner]"}] \
    }, \
    "@jupyterlab/notebook-extension:tracker": {  \ 
        "codeCellConfig": {  "lineNumbers": true, "lineWrap": true, "autoClosingBrackets": true }  \
    },   \
    "@jupyterlab/console-extension:tracker": {   \
        "promptCellConfig": { "lineNumbers": true, "lineWrap": true, "autoClosingBrackets": true }   \
    },  \
    "@jupyterlab/docmanager-extension:plugin": {   \ 
        "autosave": false   \
    }  \
}' > /opt/conda/share/jupyter/lab/settings/overrides.json 

# [Optional] If your Project.toml rarely change, uncomment this section to add them to the image.


# R packages
RUN sudo apt update && \
	sudo apt install cmake -y # CMake required by factoextra
RUN conda install -c conda-forge r-remotes=2.5.0
RUN R -e 'library(remotes) \
	remotes::install_version("socialmixr", version = "0.3.2",repos = "https://cloud.r-project.org/") \
	remotes::install_version("hhh4contacts", version="0.13.3", repos= "https://cloud.r-project.org/") \
	remotes::install_version("factoextra", version="1.0.7", repos= "https://cloud.r-project.org/") \
	remotes::install_version("countrycode", version="1.6.0", repos= "https://cloud.r-project.org/") \
	'

# Fix the issue of old libstdc++
RUN rm -rf /opt/julia-1.9.3/lib/julia/libstdc++.so.6
RUN ln -s /opt/conda/lib/libstdc++.so.6 /opt/julia-1.9.3/lib/julia/libstdc++.so.6



USER jovyan
