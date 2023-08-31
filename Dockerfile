FROM jupyter/minimal-notebook:python-3.9

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

USER root
RUN apt-get update \
   && apt-get install  -yq --no-install-recommends gcc \
   libglu1-mesa \
   libgl1-mesa-glx \
   libxrender1 \
   libxcursor \
   libxft \
   libxinerama \
   && apt-get clean && rm -rf /var/lib/apt/lists/*
USER ${NB_USER}

COPY . ${HOME}
WORKDIR $HOME

RUN pip install .[jupyter]
RUN pip uninstall vtk -y
RUN pip install --no-cache-dir --extra-index-url https://wheels.vtk.org vtk-osmesa

RUN mamba install -c conda-forge -c set3mah mmg -y

WORKDIR $HOME/examples/jupyter_notebooks

# allow jupyterlab for ipyvtk
ENV JUPYTER_ENABLE_LAB=yes
ENV PYVISTA_USE_IPYVTK=true