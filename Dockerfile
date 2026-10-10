FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt pytest

COPY . .

# Run integrity checks during the build as requested
RUN python -m pipeline.verify_integrity

EXPOSE 8501

ENTRYPOINT ["python", "run.py"]
CMD ["serve"]
