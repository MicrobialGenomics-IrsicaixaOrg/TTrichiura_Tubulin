FROM nfcore/base:1.9
LABEL authors="Francesc Catala Moll" \
      description="Docker image containing all software requirements for the  MicrobialGenomics/TTrichiura_Tubulin pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# For Bandage: otherwise it complains about missing libGL.so.1
RUN apt-get install -y libgl1-mesa-glx && apt-get clean -y

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/trichiura-1.0a/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name trichiura-1.0a > trichiura-1.0a.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
