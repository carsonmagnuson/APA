# syntax=docker/dockerfile:1.4
FROM continuumio/miniconda3:4.12.0 AS runtime

# ----- channels & faster solver -----
RUN conda install -y -n base -c conda-forge mamba \
 && conda clean -afy

# ----- create the env with only essentials -----
RUN mamba create -y -n appenv -c conda-forge -c bioconda \
        python=3.11 \
        paml         \
        pip          \
 && conda clean -afy

# Use the env for every following RUN / CMD
SHELL ["conda", "run", "-n", "appenv", "/bin/bash", "-c"]

WORKDIR /app
# ----- copy the application code -----
COPY ./backend/ /app/backend/
# ----- install your Python deps with pip -----
RUN conda run -n appenv pip install --no-cache-dir -r /app/backend/requirements.txt

EXPOSE 8000

CMD ["conda", "run", "--no-capture-output", "-n", "appenv", "uvicorn", "backend.api_server:app", "--host", "0.0.0.0", "--port", "8000"]
