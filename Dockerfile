FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    mafft wget \
    && rm -rf /var/lib/apt/lists/*

# Install MFEprimer
RUN wget -q https://github.com/quwubin/MFEprimer-3.0/releases/download/v4.2.4/mfeprimer-4.2.4-linux-amd64.gz -O /tmp/mfeprimer.gz \
    && gunzip /tmp/mfeprimer.gz \
    && chmod +x /tmp/mfeprimer \
    && mv /tmp/mfeprimer /usr/local/bin/mfeprimer

WORKDIR /app
COPY pyproject.toml README.md ./
COPY src/ src/
RUN pip install --no-cache-dir ".[web]" matplotlib

EXPOSE 7860

ENTRYPOINT ["skipalign"]
CMD ["web"]
