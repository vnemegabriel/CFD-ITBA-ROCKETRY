# Aconcagua CFD: OpenFOAM v2412 (ESI) with cfMesh, Open MPI and the tools of
# the workflow.  Build and use: docs/guide, chapter "Installation and setup".
#
#     docker build -t aconcagua-cfd -f container/Dockerfile .
#     docker run --rm -it -v "$PWD":/work -v aconcagua-data:/data aconcagua-cfd
#
# The ESI packages exist for amd64 only; on an ARM computer (Apple silicon)
# build and run with --platform linux/amd64 (emulated, slower).
FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl \
 && curl -s https://dl.openfoam.com/add-debian-repo.sh | bash \
 && apt-get install -y --no-install-recommends \
      openfoam2412-default \
      rsync xz-utils gawk time git python3 python3-matplotlib python3-pip gnuplot-nox \
 && pip3 install --no-cache-dir --break-system-packages jinja2-cli \
 && rm -rf /var/lib/apt/lists/*

# Meshes and runs go to /data (a volume), the repository to /work.
ENV LC_ALL=C.UTF-8
RUN useradd -m cfd && mkdir -p /data /work && chown cfd /data /work \
 && echo 'source /usr/lib/openfoam/openfoam2412/etc/bashrc' >> /home/cfd/.bashrc \
 && ln -s /data/meshes /home/cfd/meshes && ln -s /data/runs /home/cfd/runs
USER cfd
WORKDIR /work/flightModel
CMD ["bash", "-lc", "mkdir -p /data/meshes /data/runs && exec bash -l"]
