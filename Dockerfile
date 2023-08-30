# ARG BASE_IMAGE_TAG=latest
# FROM ghcr.io/pyvista/pyvista:$BASE_IMAGE_TAG

FROM jupyter/minimal-notebook:python-3.9

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

USER root
RUN apt-get update \
   && apt-get install  -yq --no-install-recommends git gcc \
   libglu1-mesa \
   libgl1-mesa-glx \
   libxrender1 \
   #  libfontconfig1 \
   #  libosmesa6 \
   #  xvfb \
   && apt-get clean && rm -rf /var/lib/apt/lists/*
USER ${NB_USER}

COPY . ${HOME}
WORKDIR $HOME

# RUN pip install microgen[jupyter] \
RUN pip install "microgen[jupyter] @ git+https://github.com/3MAH/microgen.git@fix_binder"
# && pip uninstall vtk -y \
# && pip install vtk-osmesa --extra-index-url https://wheels.vtk.org

RUN pip list --format=freeze > requirements.txt

# allow jupyterlab for ipyvtk
# ENV JUPYTER_ENABLE_LAB=yes
# ENV PYVISTA_USE_IPYVTK=true