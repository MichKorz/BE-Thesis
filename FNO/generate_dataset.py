#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
generate_dataset.py - Skrypt do masowego generowania symulacji LWE na podstawie
ciągu Haltona oraz modulacji fazy w szereg Taylora (TaylorSLMModulator).
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from typing import List, Tuple

import numpy as np
from scipy.stats import qmc

C_NM_FS = 299.792458


# ==============================================================================
# 1. MODULATOR FAZY (Zoptymalizowany TaylorSLMModulator z buforowaniem w RAM)
# ==============================================================================

class PhaseModulatorBase:
    def __init__(self, input_path: str, output_path: str, **kwargs):
        self.input_path = input_path
        self.output_path = output_path
        # Wczytanie pliku bazowego raz do pamięci RAM, aby uniknąć 30 000 odczytów z dysku
        self._cached_data = np.loadtxt(self.input_path)
        wavelengths_nm = self._cached_data[:, 0]

        self._cached_omega = np.zeros_like(wavelengths_nm)
        valid_mask = wavelengths_nm != 0
        self._cached_omega[valid_mask] = (2.0 * np.pi * C_NM_FS) / wavelengths_nm[valid_mask]

    def _load_speck(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return (
            self._cached_data,
            self._cached_omega,
            self._cached_data[:, 3],
            self._cached_data[:, 4],
        )

    def _write_speck(self, original_data: np.ndarray, new_re: np.ndarray, new_im: np.ndarray) -> None:
        output_data = np.copy(original_data)
        output_data[:, 3] = new_re
        output_data[:, 4] = new_im
        np.savetxt(self.output_path, output_data, fmt="%1.6e")

    def apply_phase(self, vector: List[float]):
        raise NotImplementedError()


class TaylorSLMModulator(PhaseModulatorBase):
    def __init__(
        self,
        input_path: str,
        output_path: str,
        lambda_min: float,
        lambda_max: float,
        num_stripes: int,
        lambda_0: float,
        bounds: List[List[float]],
        **kwargs,
    ):
        super().__init__(input_path, output_path, **kwargs)
        self.lambda_min = float(lambda_min)
        self.lambda_max = float(lambda_max)
        self.num_stripes = int(num_stripes)
        self.lambda_0 = float(lambda_0)
        self.bounds = bounds

    def apply_phase(self, vector: np.ndarray):
        orig_data, omega, re, im = self._load_speck()
        wavelengths_nm = orig_data[:, 0]

        omega_0 = (2.0 * np.pi * C_NM_FS) / self.lambda_0

        in_bounds = (
            (wavelengths_nm >= self.lambda_min)
            & (wavelengths_nm <= self.lambda_max)
            & (wavelengths_nm != 0)
        )

        new_re = np.copy(re)
        new_im = np.copy(im)
        final_phase = np.zeros_like(omega)

        if np.any(in_bounds):
            bins = np.linspace(self.lambda_min, self.lambda_max, self.num_stripes + 1)

            for i in range(self.num_stripes):
                mask = (wavelengths_nm >= bins[i]) & (wavelengths_nm <= bins[i + 1]) & in_bounds
                if np.any(mask):
                    lambda_center = 0.5 * (bins[i] + bins[i + 1])
                    omega_center = (2.0 * np.pi * C_NM_FS) / lambda_center

                    delta_omega = omega_center - omega_0
                    taylor_phase = 0.0

                    for order_idx, coeff in enumerate(vector):
                        n = order_idx + 2
                        taylor_phase += (coeff / math.factorial(n)) * (delta_omega ** n)

                    wrapped_phase = taylor_phase % (2.0 * np.pi)
                    final_phase[mask] = wrapped_phase

            cos_phi = np.cos(final_phase[in_bounds])
            sin_phi = np.sin(final_phase[in_bounds])

            new_re[in_bounds] = re[in_bounds] * cos_phi - im[in_bounds] * sin_phi
            new_im[in_bounds] = re[in_bounds] * sin_phi + im[in_bounds] * cos_phi

        self._write_speck(orig_data, new_re, new_im)
        return wavelengths_nm, final_phase


# ==============================================================================
# 2. SILNIK SYMULATORA LWE
# ==============================================================================

class LWESimulator:
    def __init__(self, binary_path: str, config_path: str, debug_mode: bool = False):
        self.binary_path = str(Path(binary_path).resolve())
        self.config_path = str(Path(config_path).resolve())
        self.debug_mode = bool(debug_mode)

    def run(self) -> bool:
        binary_dir = os.path.dirname(self.binary_path)

        if self.debug_mode:
            print("\n=== [LWE Symulator] Start diagnostyki ===")
            print(f"Katalog roboczy (CWD): {binary_dir}")
            print(f"Binarne uruchomienie : {self.binary_path}")
            print(f"Plik konfiguracyjny  : {self.config_path}")

            result = subprocess.run(
                [self.binary_path, self.config_path],
                cwd=binary_dir,
                check=False,
                capture_output=True,
                text=True,
            )
            print(f"[LWE Symulator] Zakończono (Kod: {result.returncode})")
            if result.stdout and result.stdout.strip():
                print("--- LWE STDOUT ---\n" + result.stdout.strip())
            if result.stderr and result.stderr.strip():
                print("--- LWE STDERR ---\n" + result.stderr.strip())
            print("=== [LWE Symulator] Koniec diagnostyki ===\n")
            return result.returncode == 0
        else:
            try:
                subprocess.run(
                    [self.binary_path, self.config_path],
                    cwd=binary_dir,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                return True
            except subprocess.CalledProcessError as e:
                print(f"\n[Błąd Symulatora LWE] Proces zakończył się błędem (Kod: {e.returncode})")
                return False


# ==============================================================================
# 3. GENERATOR CIĄGU HALTONA I ZARZĄDZANIE KONFIGURACJĄ LWE
# ==============================================================================

def get_halton_vector_for_index(
    index: int, bounds: List[List[float]], seed: int = 42, scramble: bool = True
) -> np.ndarray:
    """
    Zwraca deterministyczny wektor współczynników (c2, c3, c4, ...) dla zadanego
    indeksu `index` zgodnie z jednostajnym ciągiem Haltona.
    """
    bounds_arr = np.asarray(bounds, dtype=float)
    dim = bounds_arr.shape[0]

    sampler = qmc.Halton(d=dim, scramble=scramble, seed=seed)
    if index > 0:
        sampler.fast_forward(index)

    u = sampler.random(n=1)[0]
    lows = bounds_arr[:, 0]
    highs = bounds_arr[:, 1]
    return lows + u * (highs - lows)


def generate_all_halton_vectors(
    total_sims: int, bounds: List[List[float]], seed: int = 42, scramble: bool = True
) -> np.ndarray:
    """
    Generuje tablicę (total_sims, dim) wszystkich wektorów Haltona od 0 do total_sims - 1.
    """
    bounds_arr = np.asarray(bounds, dtype=float)
    dim = bounds_arr.shape[0]

    sampler = qmc.Halton(d=dim, scramble=scramble, seed=seed)
    u = sampler.random(n=total_sims)
    lows = bounds_arr[:, 0]
    highs = bounds_arr[:, 1]
    return lows + u * (highs - lows)


def write_lwe_config_for_index(
    template_lines: List[str],
    target_config_path: Path,
    output_base_path: Path,
    modulated_speck_path: Path,
) -> None:
    """
    Aktualizuje w szablonie LWE:
      - Linia 67: 'Output base path: <ścieżka_bez_rozszerzenia>'
      - Linia 70: 'Field 1 file path: <ścieżka_do_zmodyfikowanego_Speck.dat>'
    i zapisuje gotowy plik konfiguracyjny dla bieżącej symulacji.
    """
    new_lines = []
    replaced_output = False
    replaced_field1 = False

    for line in template_lines:
        if line.startswith("Output base path:"):
            new_lines.append(f"Output base path: {output_base_path.resolve()}\n")
            replaced_output = True
        elif line.startswith("Field 1 file path:"):
            new_lines.append(f"Field 1 file path: {modulated_speck_path.resolve()}\n")
            replaced_field1 = True
        else:
            new_lines.append(line)

    if not replaced_output or not replaced_field1:
        sys.exit(
            "[BŁĄD KRYTYCZNY] W szablonie LWE nie znaleziono klucza "
            "'Output base path:' (linia 67) lub 'Field 1 file path:' (linia 70)!"
        )

    with open(target_config_path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


def is_valid_existing_zip(zip_path: Path) -> bool:
    """
    Sprawdza, czy plik .zip istnieje, jest niepusty i ma poprawną strukturę archiwum ZIP.
    Chroni przed pominięciem uszkodzonego pliku po wcześniejszym przerwaniu pracy.
    """
    if not zip_path.is_file() or zip_path.stat().st_size == 0:
        return False
    return zipfile.is_zipfile(zip_path)


def cleanup_partial_outputs(output_base: Path) -> None:
    """
    Usuwa niekompletne pliki z bieżącego indeksu w razie przerwania (Ctrl+C) lub błędu LWE.
    """
    for ext in [".zip", ".txt", "_Ext.dat", "_spectrum.dat", "_Opt.dat"]:
        candidate = Path(str(output_base) + ext)
        if candidate.is_file():
            try:
                candidate.unlink()
            except OSError:
                pass


# ==============================================================================
# 4. GŁÓWNA PĘTLA GENERUJĄCA
# ==============================================================================

def run_generation_pipeline(config_path: str, inspect_idx: int = None) -> None:
    cfg_file = Path(config_path)
    if not cfg_file.is_file():
        sys.exit(f"[BŁĄD] Plik konfiguracyjny JSON nie istnieje: {config_path}")

    with open(cfg_file, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    # Pobranie parametrów z JSON
    paths_cfg = cfg["paths"]
    samp_cfg = cfg["sampling"]
    mod_cfg = cfg["modulator"]
    debug_mode = cfg.get("debug_mode", False)

    bounds = mod_cfg["bounds"]
    seed = int(samp_cfg.get("halton_seed", 42))
    scramble = bool(samp_cfg.get("scramble", True))
    total_sims = int(samp_cfg["total_simulations"])
    mode = samp_cfg.get("mode", "all").strip().lower()

    # Tryb inspekcji pojedynczego indeksu (bez uruchamiania LWE)
    if inspect_idx is not None:
        vec = get_halton_vector_for_index(inspect_idx, bounds, seed=seed, scramble=scramble)
        print("=" * 68)
        print(f" INSPEKCJA INDEKSU HALTONA: {inspect_idx} (plik: {inspect_idx}.zip)")
        print("=" * 68)
        for i, val in enumerate(vec):
            order = i + 2
            print(f" Współczynnik c_{order} (rząd {order}) : {val:+.6f} fs^{order}  (zakres: {bounds[i]})")
        print("=" * 68)
        return

    # Walidacja ścieżek wejściowych
    lwe_binary = Path(paths_cfg["lwe_binary"]).resolve()
    lwe_template = Path(paths_cfg["lwe_template_config"]).resolve()
    base_speck = Path(paths_cfg["base_speck_file"]).resolve()
    output_dir = Path(paths_cfg["output_simulations_dir"]).resolve()
    work_dir = Path(paths_cfg["working_dir"]).resolve()

    if not lwe_binary.is_file():
        sys.exit(f"[BŁĄD] Nie znaleziono pliku wykonywalnego LWE: {lwe_binary}")
    if not lwe_template.is_file():
        sys.exit(f"[BŁĄD] Nie znaleziono szablonu konfiguracji LWE: {lwe_template}")
    if not base_speck.is_file():
        sys.exit(f"[BŁĄD] Nie znaleziono bazowego pliku Speck.dat: {base_speck}")

    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    with open(lwe_template, "r", encoding="utf-8") as f:
        template_lines = f.readlines()

    # Wyznaczenie listy indeksów dla tego komputera
    if mode == "even":
        target_indices = list(range(0, total_sims, 2))
    elif mode == "odd":
        target_indices = list(range(1, total_sims, 2))
    elif mode == "all":
        target_indices = list(range(0, total_sims))
    else:
        sys.exit(f"[BŁĄD] Nieznany tryb 'mode': '{mode}'. Dozwolone: 'even', 'odd', 'all'.")

    # Wygenerowanie macierzy współczynników Haltona
    all_vectors = generate_all_halton_vectors(total_sims, bounds, seed=seed, scramble=scramble)

    # Pliki robocze (osobne dla każdego trybu even/odd/all, aby nie kolidowały)
    active_speck_path = work_dir / f"Speck_active_{mode}.dat"
    active_lwe_cfg_path = work_dir / f"LWE_active_{mode}.txt"

    modulator = TaylorSLMModulator(
        input_path=str(base_speck),
        output_path=str(active_speck_path),
        lambda_min=mod_cfg["lambda_min"],
        lambda_max=mod_cfg["lambda_max"],
        num_stripes=mod_cfg["num_stripes"],
        lambda_0=mod_cfg["lambda_0"],
        bounds=bounds,
    )

    simulator = LWESimulator(
        binary_path=str(lwe_binary),
        config_path=str(active_lwe_cfg_path),
        debug_mode=debug_mode,
    )

    print("=" * 76)
    print(" START GENERATORA ZBIORU TRENINGOWEGO LWE (CIĄG HALTONA)")
    print("=" * 76)
    print(f" Całkowita docelowa liczba symulacji : {total_sims}")
    print(f" Tryb pracy tego komputera (mode)    : {mode.upper()} ({len(target_indices)} indeksów)")
    print(f" Katalog docelowy wynikowych .zip    : {output_dir}")
    print(f" Zakresy Taylora (c2, c3, c4)        : {bounds}")
    print("=" * 76)

    completed_count = 0
    skipped_count = 0
    start_time = time.time()

    for pos, sim_idx in enumerate(target_indices, start=1):
        output_base = output_dir / str(sim_idx)
        expected_zip = output_dir / f"{sim_idx}.zip"

        # Pomijanie już policzonych i poprawnych plików .zip
        if is_valid_existing_zip(expected_zip):
            skipped_count += 1
            continue

        coeff_vector = all_vectors[sim_idx]

        try:
            t0 = time.time()

            # 1. Nałożenie fazy Taylora i zapis roboczego pliku Speck
            modulator.apply_phase(coeff_vector)

            # 2. Aktualizacja linii 67 (Output base path) i 70 (Field 1 file path) w configu LWE
            write_lwe_config_for_index(
                template_lines=template_lines,
                target_config_path=active_lwe_cfg_path,
                output_base_path=output_base,
                modulated_speck_path=active_speck_path,
            )

            # 3. Uruchomienie symulatora LWE
            success = simulator.run()
            dt_sim = time.time() - t0

            if not success or not is_valid_existing_zip(expected_zip):
                print(
                    f"[OSTRZEŻENIE] Symulacja idx={sim_idx} nie utworzyła poprawnego pliku "
                    f"{expected_zip.name}. Czyszczenie niekompletnych plików..."
                )
                cleanup_partial_outputs(output_base)
                continue

            completed_count += 1
            elapsed = time.time() - start_time
            avg_time = elapsed / completed_count
            remaining_todo = len(target_indices) - pos
            eta_min = (remaining_todo * avg_time) / 60.0

            print(
                f"[{pos}/{len(target_indices)}] idx={sim_idx:05d}.zip | "
                f"c2={coeff_vector[0]:+8.1f}, c3={coeff_vector[1]:+9.1f}, c4={coeff_vector[2]:+10.1f} | "
                f"czas: {dt_sim:.2f}s | ETA: {eta_min:.1f} min"
            )

        except KeyboardInterrupt:
            print(
                f"\n[PRZERWANO PRZEZ UŻYTKOWNIKA] Czyszczenie niedokończonej symulacji idx={sim_idx}..."
            )
            cleanup_partial_outputs(output_base)
            sys.exit(0)

    print("=" * 76)
    print(f" ZAKOŃCZONO! Nowo policzone: {completed_count} | Pominięte (istniejące): {skipped_count}")
    print("=" * 76)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Masowe generowanie symulacji LWE z użyciem ciągu Haltona."
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="./sim_config.json",
        help="Ścieżka do pliku konfiguracyjnego JSON (domyślnie: ./sim_config.json).",
    )
    parser.add_argument(
        "--inspect-index",
        type=int,
        default=None,
        metavar="IDX",
        help="Wyświetla dokładne współczynniki (c2, c3, c4) dla podanego indeksu symulacji i kończy działanie.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_generation_pipeline(config_path=args.config, inspect_idx=args.inspect_index)
