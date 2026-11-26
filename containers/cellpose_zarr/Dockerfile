# Adapted from the cellposesam-zarr recipe by the EMBL Centre for Bioimage Analysis
# (https://git.embl.org/grp-cba/containers).

# ---- Build arguments ----
ARG PYTHON_VERSION=3.12
ARG CELLPOSE_VERSION=4.1.1
ARG ZARR_VERSION=3.1.6
ARG NUMPY_VERSION=2.4.4
# ZarrSmith reads the input regions and writes the label stores (cellpose_sam_zarr.py).
# Pinned to a commit of branch consolidated_translation (read_region, label writing),
# which is not released yet: pin a release once there is one.
ARG ZARRSMITH_REF=52b99a296680cbb147e41700674d83d2a7e47953

FROM python:${PYTHON_VERSION}-slim

ARG CELLPOSE_VERSION
ARG ZARR_VERSION
ARG NUMPY_VERSION
ARG ZARRSMITH_REF

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore

# procps provides `ps`, which Nextflow needs to collect task metrics;
# git lets pip install ZarrSmith from its repository
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    libgl1 \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    numpy==${NUMPY_VERSION} \
    zarr==${ZARR_VERSION} \
    "zarrsmith @ git+https://github.com/embl-mif/zarrsmith.git@${ZARRSMITH_REF}" \
    aiohttp requests \
    cellpose==${CELLPOSE_VERSION}

CMD ["python"]
