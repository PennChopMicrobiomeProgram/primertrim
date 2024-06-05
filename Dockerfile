FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home

COPY . .

# Install environment
RUN mamba create --name primertrim -c conda-forge -c bioconda vsearch

ENV PATH="/opt/conda/envs/primertrim/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "--no-capture-output", "-n", "primertrim", "/bin/bash", "-c"]

RUN pip install /home/

RUN echo "Python: $(python --version), Conda: $(conda --version), Vsearch: $(vsearch -v)" > installed_packages.txt

# Run
CMD "bash"