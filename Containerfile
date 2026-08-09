# Używamy oficjalnego obrazu NVIDIA CUDA Developer opartego na Ubuntu 24.04
# CUDA 12.8 to minimalna wersja natywnie wspierająca architekturę Blackwell (sm_120)
FROM nvidia/cuda:12.8.0-devel-ubuntu24.04

# Wyłączenie interaktywnych monitów przy instalacji pakietów
ENV DEBIAN_FRONTEND=noninteractive

# 1. Instalacja wymaganych narzędzi systemowych i kompilatorów
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    python3-dev \
    python3-pip \
    python3-venv \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# 2. Przygotowanie środowiska wirtualnego Pythona i instalacja PyTorch (z obsługą Blackwell/cu128)
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --pre torch --index-url https://download.pytorch.org/whl/nightly/cu128 && \
    pip install --no-cache-dir numpy scipy pandas LightwaveExplorer matplotlib

# Domyślny katalog roboczy w kontenerze
WORKDIR /workspace
