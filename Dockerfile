FROM jupyter/base-notebook

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

USER root
RUN apt-get update \
   && apt-get install  -yq --no-install-recommends git \
   libglu1-mesa \
   #  libgl1-mesa-glx \
   #  libfontconfig1 \
   #  libxrender1 \
   #  libosmesa6 \
   #  xvfb \
   && apt-get clean && rm -rf /var/lib/apt/lists/*
USER ${NB_USER}


# RUN pip install microgen[jupyter] \
RUN pip install "microgen[jupyter] @ git+https://github.com/3MAH/microgen.git@fix_binder" \
   && pip uninstall vtk -y \
   && pip install vtk-osmesa --extra-index-url https://wheels.vtk.org


WORKDIR $HOME

# allow jupyterlab for ipyvtk
ENV JUPYTER_ENABLE_LAB=yes
ENV PYVISTA_USE_IPYVTK=true