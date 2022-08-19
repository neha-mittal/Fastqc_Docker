FROM continuumio/miniconda3
RUN conda update conda
COPY requirements.txt .
RUN pip install -r requirements.txt


WORKDIR /app
COPY . /app
CMD ["python",  "fast_qc1.py"]