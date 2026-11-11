#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
speckutils - Narzędzie do analizy i manipulacji plikami widmowymi FROG (Speck.dat)
oraz przygotowywania profili impulsów dla symulatora Lightwave Explorer (LWE).
"""

import argparse
import sys
from pathlib import Path
import numpy as np

# ==============================================================================
# STAŁE FIZYCZNE
# ==============================================================================
C_NM_FS = 299.792458  # Prędkość światła w próżni [nm / fs]

# ==============================================================================
# ZMIENNE GLOBALNE (Aktualizowane po wczytaniu pliku oraz po każdej modyfikacji)
# ==============================================================================

# 1. Surowe / robocze dane z pliku (tablice 1D o rozmiarze N_GRID)
FILE_PATH = None
RAW_DATA = None
COL_WAVELENGTH_NM = None      # Kolumna 0: Długość fali [nm]
COL_INTENSITY = None          # Kolumna 1: Intensywność I(lambda) znormalizowana do 1.0
COL_PHASE_RAD = None          # Kolumna 2: Faza spektralna z pliku [rad]
COL_REAL = None               # Kolumna 3: Część rzeczywista pola E(nu)
COL_IMAG = None               # Kolumna 4: Część urojona pola E(nu)
COMPLEX_FIELD = None          # Zespolone pole elektryczne: Re + 1j * Im

# 2. Parametry siatki FFT i transformacji czasowo-częstotliwościowej
N_GRID = 0
IS_POWER_OF_TWO = False
FREQ_THZ = None               # Oś częstotliwości optycznej nu = c / lambda [THz]
OMEGA_RAD_FS = None           # Oś częstości kołowej omega = 2 * pi * nu [rad/fs]
DF_THZ = 0.0                  # Krok częstotliwościowy siatki FFT (Delta nu) [THz]
DOMEGA_RAD_FS = 0.0           # Krok częstości kołowej (Delta omega) [rad/fs]
FREQ_SPAN_THZ = 0.0           # Całkowita rozpiętość okna częstotliwościowego [THz]
DT_FS = 0.0                   # Krok czasowy w dziedzinie czasu (Delta t = 1 / Span) [fs]
TIME_SPAN_FS = 0.0            # Całkowity rozmiar okna czasowego (T_window = N * dt) [fs]
NYQUIST_THZ = 0.0             # Częstotliwość Nyquista (Span / 2) [THz]

# 3. Parametry częstotliwości nośnej (środek siatki FFT)
CENTER_GRID_IDX = 0
CARRIER_WL_NM = 0.0
CARRIER_FREQ_THZ = 0.0

# 4. Podział na domenę fizyczną (dodatnie częstotliwości) i niefizyczną (ujemne)
POS_MASK = None
NEG_MASK = None
N_POS_POINTS = 0
N_NEG_POINTS = 0
SINGULARITY_IDX = None
POS_WL_MIN_NM = 0.0
POS_WL_MAX_NM = 0.0

# 5. Charakterystyka widmowa impulsu
PEAK_IDX = 0
PEAK_WL_NM = 0.0
PEAK_FREQ_THZ = 0.0
PEAK_OMEGA_RAD_FS = 0.0
PEAK_STEP_WL_NM = 0.0         # Lokalny krok długości fali w piku widma [nm]

FIELD_PEAK_IDX = 0
FIELD_PEAK_WL_NM = 0.0
FIELD_PEAK_FREQ_THZ = 0.0

COM_FREQ_THZ = 0.0
COM_WL_NM = 0.0
COM_OMEGA_RAD_FS = 0.0

FWHM_WL_NM = 0.0
FWHM_FREQ_THZ = 0.0
RMS_BW_THZ = 0.0

# 6. Parametry Fazy, Phase Blankingu i weryfikacji Jakobianu
RAW_PHASE_RAD = None
UNWRAPPED_PHASE_RAD = None
FROG_PHASE_OFFSET_RAD = 0.0   # Przesunięcie stałe między kolumną 2 a unwrap(atan2(Im, Re))
PHASE_BLANKING_MASK = None
N_ACTIVE_PHASE_PTS = 0
N_EXACT_ZERO_PTS = 0          # Liczba punktów o dokładnie zerowej amplitudzie (Re=0 i Im=0)
EST_BLANKING_THRESH = 1e-4
ACTIVE_WL_MIN_NM = 0.0
ACTIVE_WL_MAX_NM = 0.0
MAX_PHASE_STEP_ACTIVE = 0.0
JACOBIAN_REL_ERROR = 0.0

# 7. Szacowane parametry w dziedzinie czasu (z IFFT kolumn Re/Im)
FWHM_TIME_ACTUAL_FS = 0.0
FWHM_TIME_TL_FS = 0.0


# ==============================================================================
# FUNKCJE POMOCNICZE I OBLICZENIOWE
# ==============================================================================

def is_power_of_two(n: int) -> bool:
    return (n > 0) and ((n & (n - 1)) == 0)


def parse_power_or_threshold(val: float) -> float:
    """
    Konwertuje argument użytkownika na próg intensywności w zakresie (0, 1).
    Obsługuje:
      - ujemny wykładnik potęgi 10: np. -3 -> 10^(-3) = 1e-3
      - dodatni wykładnik potęgi 10: np. 3 -> 10^(-3) = 1e-3
      - bezpośredni ułamek dziesiętny: np. 1e-3 lub 0.001 -> 1e-3
    """
    if val < 0:
        thresh = 10.0 ** val
    elif val >= 1.0:
        thresh = 10.0 ** (-val)
    else:
        thresh = val

    if not (0.0 < thresh < 1.0):
        sys.exit(
            f"[BŁĄD] Nieprawidłowa wartość progu ({val} -> {thresh}). "
            f"Podaj wykładnik potęgi 10 (np. -3 lub 3 dla 10^-3) lub ułamek (np. 1e-3)."
        )
    return float(thresh)


def calculate_fwhm(x: np.ndarray, y: np.ndarray) -> float:
    """
    Oblicza szerokość w połowie wysokości (FWHM) pomiędzy skrajnymi przekroczeniami
    poziomu 0.5 * max(y) z liniową interpolacją podpikselową.
    """
    if len(x) < 2 or np.max(y) <= 0:
        return 0.0

    if x[0] > x[-1]:
        x = x[::-1]
        y = y[::-1]

    half_max = 0.5 * np.max(y)
    above = y >= half_max
    indices = np.where(above)[0]

    if len(indices) == 0:
        return 0.0

    i_left = indices[0]
    if i_left > 0:
        x0, x1 = x[i_left - 1], x[i_left]
        y0, y1 = y[i_left - 1], y[i_left]
        x_left = x0 + (half_max - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
    else:
        x_left = x[0]

    i_right = indices[-1]
    if i_right < len(x) - 1:
        x0, x1 = x[i_right], x[i_right + 1]
        y0, y1 = y[i_right], y[i_right + 1]
        x_right = x0 + (half_max - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x0
    else:
        x_right = x[-1]

    return float(abs(x_right - x_left))


def update_global_diagnostics() -> None:
    """
    Przelicza wszystkie globalne zmienne diagnostyczne na podstawie aktualnego
    stanu tablic COL_WAVELENGTH_NM, COL_INTENSITY, COL_PHASE_RAD, COL_REAL, COL_IMAG.
    """
    global COMPLEX_FIELD
    global N_GRID, IS_POWER_OF_TWO, FREQ_THZ, OMEGA_RAD_FS, DF_THZ, DOMEGA_RAD_FS
    global FREQ_SPAN_THZ, DT_FS, TIME_SPAN_FS, NYQUIST_THZ
    global CENTER_GRID_IDX, CARRIER_WL_NM, CARRIER_FREQ_THZ
    global POS_MASK, NEG_MASK, N_POS_POINTS, N_NEG_POINTS, SINGULARITY_IDX
    global POS_WL_MIN_NM, POS_WL_MAX_NM
    global PEAK_IDX, PEAK_WL_NM, PEAK_FREQ_THZ, PEAK_OMEGA_RAD_FS, PEAK_STEP_WL_NM
    global FIELD_PEAK_IDX, FIELD_PEAK_WL_NM, FIELD_PEAK_FREQ_THZ
    global COM_FREQ_THZ, COM_WL_NM, COM_OMEGA_RAD_FS
    global FWHM_WL_NM, FWHM_FREQ_THZ, RMS_BW_THZ
    global RAW_PHASE_RAD, UNWRAPPED_PHASE_RAD, FROG_PHASE_OFFSET_RAD
    global PHASE_BLANKING_MASK, N_ACTIVE_PHASE_PTS, N_EXACT_ZERO_PTS
    global EST_BLANKING_THRESH, ACTIVE_WL_MIN_NM, ACTIVE_WL_MAX_NM, MAX_PHASE_STEP_ACTIVE
    global JACOBIAN_REL_ERROR, FWHM_TIME_ACTUAL_FS, FWHM_TIME_TL_FS

    COMPLEX_FIELD = COL_REAL + 1j * COL_IMAG

    # 1. Parametry siatki FFT
    N_GRID = len(COL_WAVELENGTH_NM)
    IS_POWER_OF_TWO = is_power_of_two(N_GRID)

    FREQ_THZ = (C_NM_FS / COL_WAVELENGTH_NM) * 1e3
    OMEGA_RAD_FS = (2.0 * np.pi * C_NM_FS) / COL_WAVELENGTH_NM

    diffs_thz = np.diff(FREQ_THZ)
    DF_THZ = float(np.median(np.abs(diffs_thz)))
    DOMEGA_RAD_FS = float(2.0 * np.pi * DF_THZ * 1e-3)

    FREQ_SPAN_THZ = N_GRID * DF_THZ
    NYQUIST_THZ = 0.5 * FREQ_SPAN_THZ
    DT_FS = 1000.0 / FREQ_SPAN_THZ if FREQ_SPAN_THZ > 0 else 0.0
    TIME_SPAN_FS = N_GRID * DT_FS

    # 2. Częstotliwość nośna (środek tablicy FFT)
    CENTER_GRID_IDX = N_GRID // 2
    CARRIER_WL_NM = float(COL_WAVELENGTH_NM[CENTER_GRID_IDX])
    CARRIER_FREQ_THZ = float(FREQ_THZ[CENTER_GRID_IDX])

    # 3. Podział na częstotliwości dodatnie i ujemne
    POS_MASK = COL_WAVELENGTH_NM > 0
    NEG_MASK = COL_WAVELENGTH_NM < 0
    N_POS_POINTS = int(np.sum(POS_MASK))
    N_NEG_POINTS = int(np.sum(NEG_MASK))

    sign_changes = np.where(np.diff(np.sign(COL_WAVELENGTH_NM)) != 0)[0]
    SINGULARITY_IDX = int(sign_changes[0]) if len(sign_changes) > 0 else None

    if N_POS_POINTS > 0:
        POS_WL_MIN_NM = float(np.min(COL_WAVELENGTH_NM[POS_MASK]))
        POS_WL_MAX_NM = float(np.max(COL_WAVELENGTH_NM[POS_MASK]))

    # 4. Charakterystyka widmowa (w domenie fizycznej lambda > 0)
    field_mod_sq = np.abs(COMPLEX_FIELD) ** 2
    N_EXACT_ZERO_PTS = int(np.sum(field_mod_sq == 0.0))

    PEAK_IDX = int(np.argmax(COL_INTENSITY))
    PEAK_WL_NM = float(COL_WAVELENGTH_NM[PEAK_IDX])
    PEAK_FREQ_THZ = float(FREQ_THZ[PEAK_IDX])
    PEAK_OMEGA_RAD_FS = float(OMEGA_RAD_FS[PEAK_IDX])
    if 0 < PEAK_IDX < N_GRID - 1:
        PEAK_STEP_WL_NM = float(
            0.5 * abs(COL_WAVELENGTH_NM[PEAK_IDX + 1] - COL_WAVELENGTH_NM[PEAK_IDX - 1])
        )
    else:
        PEAK_STEP_WL_NM = 0.0

    FIELD_PEAK_IDX = int(np.argmax(field_mod_sq))
    FIELD_PEAK_WL_NM = float(COL_WAVELENGTH_NM[FIELD_PEAK_IDX])
    FIELD_PEAK_FREQ_THZ = float(FREQ_THZ[FIELD_PEAK_IDX])

    if N_POS_POINTS > 1:
        pos_freq = FREQ_THZ[POS_MASK]
        pos_field_sq = field_mod_sq[POS_MASK]

        weight_sum = np.sum(pos_field_sq)
        if weight_sum > 0:
            COM_FREQ_THZ = float(np.sum(pos_freq * pos_field_sq) / weight_sum)
            COM_WL_NM = float((C_NM_FS * 1e3) / COM_FREQ_THZ)
            COM_OMEGA_RAD_FS = float(2.0 * np.pi * COM_FREQ_THZ * 1e-3)
            RMS_BW_THZ = float(
                np.sqrt(np.sum(((pos_freq - COM_FREQ_THZ) ** 2) * pos_field_sq) / weight_sum)
            )

        valid_wl_mask = POS_MASK & (COL_WAVELENGTH_NM < 20000.0)
        if np.sum(valid_wl_mask) > 2:
            FWHM_WL_NM = calculate_fwhm(
                COL_WAVELENGTH_NM[valid_wl_mask], COL_INTENSITY[valid_wl_mask]
            )
        FWHM_FREQ_THZ = calculate_fwhm(pos_freq, pos_field_sq)

    # 5. Analiza fazy, Phase Blankingu i Jakobianu
    RAW_PHASE_RAD = np.arctan2(COL_IMAG, COL_REAL)
    UNWRAPPED_PHASE_RAD = np.unwrap(RAW_PHASE_RAD)

    PHASE_BLANKING_MASK = COL_PHASE_RAD != 0.0
    N_ACTIVE_PHASE_PTS = int(np.sum(PHASE_BLANKING_MASK))

    if N_ACTIVE_PHASE_PTS > 0:
        active_wls = COL_WAVELENGTH_NM[PHASE_BLANKING_MASK]
        ACTIVE_WL_MIN_NM = float(np.min(active_wls))
        ACTIVE_WL_MAX_NM = float(np.max(active_wls))

        FROG_PHASE_OFFSET_RAD = float(
            COL_PHASE_RAD[PEAK_IDX] - UNWRAPPED_PHASE_RAD[PEAK_IDX]
        )

        min_active_int = float(np.min(COL_INTENSITY[PHASE_BLANKING_MASK]))
        zero_phase_mask = ~PHASE_BLANKING_MASK & (COL_INTENSITY > 0)
        max_blanked_int = (
            float(np.max(COL_INTENSITY[zero_phase_mask]))
            if np.any(zero_phase_mask)
            else min_active_int
        )
        EST_BLANKING_THRESH = 0.5 * (min_active_int + max_blanked_int)

        active_indices = np.where(PHASE_BLANKING_MASK)[0]
        if len(active_indices) > 1:
            MAX_PHASE_STEP_ACTIVE = float(
                np.max(np.abs(np.diff(UNWRAPPED_PHASE_RAD[active_indices])))
            )
        else:
            MAX_PHASE_STEP_ACTIVE = 0.0
    else:
        fallback_mask = POS_MASK & (COL_INTENSITY >= 1e-4)
        if np.any(fallback_mask):
            active_wls = COL_WAVELENGTH_NM[fallback_mask]
            ACTIVE_WL_MIN_NM = float(np.min(active_wls))
            ACTIVE_WL_MAX_NM = float(np.max(active_wls))
            active_indices = np.where(fallback_mask)[0]
            MAX_PHASE_STEP_ACTIVE = (
                float(np.max(np.abs(np.diff(UNWRAPPED_PHASE_RAD[active_indices]))))
                if len(active_indices) > 1
                else 0.0
            )
        else:
            ACTIVE_WL_MIN_NM = 0.0
            ACTIVE_WL_MAX_NM = 0.0
            MAX_PHASE_STEP_ACTIVE = 0.0
        FROG_PHASE_OFFSET_RAD = 0.0

    jac_calc = field_mod_sq / (COL_WAVELENGTH_NM ** 2)
    if np.max(jac_calc[POS_MASK]) > 0:
        jac_calc /= np.max(jac_calc[POS_MASK])
    ref_int = COL_INTENSITY / np.max(COL_INTENSITY[POS_MASK])
    JACOBIAN_REL_ERROR = float(np.max(np.abs(jac_calc[POS_MASK] - ref_int[POS_MASK])))

    # 6. Analiza w dziedzinie czasu (IFFT z zespolonej obwiedni)
    t_axis_fs = (np.arange(N_GRID) - (N_GRID // 2)) * DT_FS

    e_time_actual = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(COMPLEX_FIELD)))
    i_time_actual = np.abs(e_time_actual) ** 2
    if np.max(i_time_actual) > 0:
        peak_t_idx = int(np.argmax(i_time_actual))
        i_time_actual_centered = np.roll(i_time_actual, (N_GRID // 2) - peak_t_idx)
        FWHM_TIME_ACTUAL_FS = calculate_fwhm(t_axis_fs, i_time_actual_centered)
    else:
        FWHM_TIME_ACTUAL_FS = 0.0

    field_tl = np.zeros_like(COMPLEX_FIELD)
    field_tl[POS_MASK] = np.abs(COMPLEX_FIELD[POS_MASK])
    e_time_tl = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(field_tl)))
    i_time_tl = np.abs(e_time_tl) ** 2
    FWHM_TIME_TL_FS = calculate_fwhm(t_axis_fs, i_time_tl)


def load_speck_file(filepath: str) -> None:
    """
    Wczytuje plik Speck.dat do pamięci i wywołuje obliczenie zmiennych globalnych.
    """
    global FILE_PATH, RAW_DATA
    global COL_WAVELENGTH_NM, COL_INTENSITY, COL_PHASE_RAD, COL_REAL, COL_IMAG

    path = Path(filepath)
    if not path.is_file():
        sys.exit(f"[BŁĄD] Plik nie istnieje: {filepath}")

    try:
        data = np.loadtxt(path)
    except Exception as e:
        sys.exit(f"[BŁĄD] Nie udało się wczytać pliku '{filepath}': {e}")

    if data.ndim != 2 or data.shape[1] != 5:
        sys.exit(
            f"[BŁĄD] Nieprawidłowy format pliku. Oczekiwano 5 kolumn, "
            f"otrzymano kształt: {data.shape}"
        )

    FILE_PATH = str(path.resolve())
    RAW_DATA = data.copy()
    COL_WAVELENGTH_NM = data[:, 0].copy()
    COL_INTENSITY = data[:, 1].copy()
    COL_PHASE_RAD = data[:, 2].copy()
    COL_REAL = data[:, 3].copy()
    COL_IMAG = data[:, 4].copy()

    update_global_diagnostics()


# ==============================================================================
# FUNKCJE MANIPULACJI DANYMI (Korzystające ze zmiennych globalnych)
# ==============================================================================

def apply_upsampling(target_n: int) -> None:
    """
    Zagęszcza siatkę FFT pliku Speck.dat do `target_n` punktów (potęga dwójki)
    metodą interpolacji liniowej we współrzędnych biegunowych (moduł + rozwinięta faza).
    """
    global COL_WAVELENGTH_NM, COL_INTENSITY, COL_PHASE_RAD, COL_REAL, COL_IMAG

    if not is_power_of_two(target_n) or target_n < 64:
        low_p2 = 2 ** int(np.floor(np.log2(max(64, target_n))))
        high_p2 = 2 ** int(np.ceil(np.log2(max(64, target_n))))
        sys.exit(
            f"[BŁĄD] Liczba punktów siatki musi być potęgą dwójki (>= 64).\n"
            f"       Podano: {target_n}. Najbliższe poprawne wartości to: {low_p2} lub {high_p2}."
        )

    old_n = N_GRID
    old_df = DF_THZ
    old_step_peak = PEAK_STEP_WL_NM

    if target_n == old_n:
        print(f"[INFO] Siatka ma już zadany rozmiar {target_n} punktów. Pominięto resampling.")
        return

    new_df_thz = DF_THZ * (old_n / target_n)
    new_center_idx = target_n // 2

    indices = np.arange(target_n, dtype=float)
    new_freq_thz = CARRIER_FREQ_THZ + (indices - new_center_idx) * new_df_thz

    zero_freq_mask = np.abs(new_freq_thz) < 1e-9
    if np.any(zero_freq_mask):
        new_freq_thz[zero_freq_mask] = 1e-9

    new_wavelength_nm = (C_NM_FS * 1e3) / new_freq_thz

    orig_amp = np.abs(COMPLEX_FIELD)
    orig_unwrapped_phase = UNWRAPPED_PHASE_RAD

    if FREQ_THZ[0] < FREQ_THZ[-1]:
        xp = FREQ_THZ
        fp_amp = orig_amp
        fp_phase = orig_unwrapped_phase
    else:
        xp = FREQ_THZ[::-1]
        fp_amp = orig_amp[::-1]
        fp_phase = orig_unwrapped_phase[::-1]

    new_amp = np.interp(new_freq_thz, xp, fp_amp)
    new_unwrapped_phase = np.interp(new_freq_thz, xp, fp_phase)

    new_real = new_amp * np.cos(new_unwrapped_phase)
    new_imag = new_amp * np.sin(new_unwrapped_phase)

    new_intensity = (new_real ** 2 + new_imag ** 2) / (new_wavelength_nm ** 2)
    max_int = np.max(new_intensity[new_wavelength_nm > 0])
    if max_int > 0:
        new_intensity /= max_int

    blanking_thresh = EST_BLANKING_THRESH if EST_BLANKING_THRESH > 0 else 1e-4
    new_phase_col = np.where(
        (new_wavelength_nm > 0) & (new_intensity >= blanking_thresh),
        new_unwrapped_phase + FROG_PHASE_OFFSET_RAD,
        0.0,
    )

    COL_WAVELENGTH_NM = new_wavelength_nm
    COL_INTENSITY = new_intensity
    COL_PHASE_RAD = new_phase_col
    COL_REAL = new_real
    COL_IMAG = new_imag

    update_global_diagnostics()

    print("=" * 76)
    print(" ZASTOSOWANO ZAGĘSZCZENIE SIATKI WIDMOWEJ (SPECTRAL UPSAMPLING)")
    print("=" * 76)
    print(f" Zmiana rozmiaru siatki FFT  : {old_n} pkt  -->  {N_GRID} pkt ({N_GRID / old_n:.1f}x)")
    print(f" Krok częstotliwości (df)    : {old_df:.6f} THz  -->  {DF_THZ:.6f} THz")
    print(f" Krok w piku widma (d_lambda): {old_step_peak:.3f} nm  -->  {PEAK_STEP_WL_NM:.3f} nm (przy {PEAK_WL_NM:.1f} nm)")
    print(f" Długość fali w środku (N/2) : {CARRIER_WL_NM:.4f} nm (zachowana bez zmian)")
    print(f" Nowe okno czasowe IFFT      : {TIME_SPAN_FS:.2f} fs (krok dt = {DT_FS:.4f} fs)")
    print("=" * 76)


def apply_zero_pulse_phase(threshold_arg: float) -> None:
    """
    Zeruje fazę spektralną samego impulsu (w punktach gdzie lambda > 0 oraz
    Intensity >= threshold), sprowadzając impuls do granicy transformaty (Transform-Limited),
    bez ruszania losowej fazy szumu tła poza impulsem.
    """
    global COL_PHASE_RAD, COL_REAL, COL_IMAG

    threshold = parse_power_or_threshold(threshold_arg)
    exponent = np.log10(threshold)

    pulse_mask = POS_MASK & (COL_INTENSITY >= threshold)
    n_zeroed = int(np.sum(pulse_mask))

    if n_zeroed == 0:
        sys.exit(f"[BŁĄD] Żaden punkt nie przekracza progu {threshold:.2e} (10^{exponent:.1f}).")

    wl_min_zeroed = float(np.min(COL_WAVELENGTH_NM[pulse_mask]))
    wl_max_zeroed = float(np.max(COL_WAVELENGTH_NM[pulse_mask]))
    old_fwhm_time = FWHM_TIME_ACTUAL_FS

    amp_in_pulse = np.sqrt(COL_REAL[pulse_mask] ** 2 + COL_IMAG[pulse_mask] ** 2)
    COL_REAL[pulse_mask] = amp_in_pulse
    COL_IMAG[pulse_mask] = 0.0
    COL_PHASE_RAD[pulse_mask] = 0.0

    update_global_diagnostics()

    print("=" * 76)
    print(" WYZEROWANO FAZĘ IMPULSU (TRANSFORM-LIMITED PULSE PHASE ZEROING)")
    print("=" * 76)
    print(f" Próg początku impulsu       : 10^({exponent:.2f}) = {threshold:.2e} (względem piku I = 1.0)")
    print(f" Liczba punktów z fazą = 0   : {n_zeroed} pkt (zakres: {wl_min_zeroed:.2f} nm .. {wl_max_zeroed:.2f} nm)")
    print(f" Szum poza impulsem (< próg) : Pozostawiony z oryginalną losową fazą (brak piku w t=0)")
    print(f" Czas trwania impulsu (FWHM) : {old_fwhm_time:.2f} fs  -->  {FWHM_TIME_ACTUAL_FS:.2f} fs "
          f"(Transform-Limit: {FWHM_TIME_TL_FS:.2f} fs)")
    print("=" * 76)


def apply_noise_windowing(
    threshold_arg: float, taper_width_thz: float = 15.0, noise_floor: float = 1e-12
) -> None:
    """
    Tłumi szum tła w skrzydłach widma poniżej progu `threshold` (względem piku I=1.0)
    oraz w niefizycznej domenie ujemnych długości fal (lambda < 0) za pomocą
    gładkiej apodyzacji cosinusowej i miękkiego dna `noise_floor`.
    """
    global COL_INTENSITY, COL_PHASE_RAD, COL_REAL, COL_IMAG

    threshold = parse_power_or_threshold(threshold_arg)

    valid_signal_mask = POS_MASK & (COL_INTENSITY >= threshold)
    signal_indices = np.where(valid_signal_mask)[0]

    if len(signal_indices) < 2:
        sys.exit(
            f"[BŁĄD] Brak wystarczającej liczby punktów powyżej progu {threshold:.2e}. "
            f"Zmniejsz próg filtrowania."
        )

    idx_start = int(signal_indices[0])
    idx_end = int(signal_indices[-1])

    f_low_thz = FREQ_THZ[idx_start]
    f_high_thz = FREQ_THZ[idx_end]
    wl_high_nm = COL_WAVELENGTH_NM[idx_start]
    wl_low_nm = COL_WAVELENGTH_NM[idx_end]

    window = np.full(N_GRID, noise_floor, dtype=float)
    window[idx_start : idx_end + 1] = 1.0

    taper_pts = int(np.ceil(taper_width_thz / DF_THZ)) if taper_width_thz > 0 else 0

    if taper_pts > 0:
        first_pos_idx = SINGULARITY_IDX + 1 if SINGULARITY_IDX is not None else 0
        left_start = max(first_pos_idx, idx_start - taper_pts)
        n_left = idx_start - left_start
        if n_left > 0:
            x_left = np.linspace(0.0, 1.0, n_left + 2)[1:-1]
            taper_l = 0.5 * (1.0 - np.cos(np.pi * x_left))
            window[left_start:idx_start] = noise_floor + (1.0 - noise_floor) * taper_l

        right_end = min(N_GRID - 1, idx_end + taper_pts)
        n_right = right_end - idx_end
        if n_right > 0:
            x_right = np.linspace(0.0, 1.0, n_right + 2)[1:-1]
            taper_r = 0.5 * (1.0 + np.cos(np.pi * x_right))
            window[idx_end + 1 : right_end + 1] = noise_floor + (1.0 - noise_floor) * taper_r

    energy_before = float(np.sum((COL_REAL[POS_MASK] ** 2) + (COL_IMAG[POS_MASK] ** 2)))

    COL_REAL = COL_REAL * window
    COL_IMAG = COL_IMAG * window

    energy_after = float(np.sum((COL_REAL[POS_MASK] ** 2) + (COL_IMAG[POS_MASK] ** 2)))
    retained_pct = (100.0 * energy_after / energy_before) if energy_before > 0 else 0.0

    new_int = (COL_REAL ** 2 + COL_IMAG ** 2) / (COL_WAVELENGTH_NM ** 2)
    max_int = np.max(new_int[POS_MASK])
    if max_int > 0:
        new_int /= max_int
    COL_INTENSITY = new_int

    phase_cutoff = min(threshold, EST_BLANKING_THRESH if EST_BLANKING_THRESH > 0 else threshold)
    COL_PHASE_RAD = np.where(POS_MASK & (COL_INTENSITY >= phase_cutoff), COL_PHASE_RAD, 0.0)

    update_global_diagnostics()

    print("=" * 76)
    print(" ZASTOSOWANO OKIENKOWANIE SZUMU (NOISE WINDOWING / APODIZATION)")
    print("=" * 76)
    print(f" Próg odcięcia intensywności : {threshold:.2e} (względem piku I = 1.0)")
    print(f" Pasmo 100% transmisji okna  : {wl_low_nm:.2f} nm .. {wl_high_nm:.2f} nm "
          f"({f_low_thz:.2f} .. {f_high_thz:.2f} THz)")
    print(f" Szerokość zbocza (Cosine)   : {taper_width_thz:.2f} THz ({taper_pts} pikseli FFT z każdej strony)")
    print(f" Poziom tłumienia poza oknem : {noise_floor:.1e} w amplitudzie ({noise_floor**2:.1e} w energii)")
    print(f" Zachowana energia widmowa   : {retained_pct:.4f}% energii fizycznej")
    print("=" * 76)


def apply_hard_zero_noise(threshold_arg: float) -> None:
    """
    [EKSPERYMENTALNE] Twarde zerowanie szumu (Hard Zeroing):
    Wszystkie punkty, dla których lambda < 0 lub Intensity < threshold,
    otrzymują dokładne 0.0 w kolumnach Intensity, Phase, Real oraz Imaginary,
    przy zachowaniu oryginalnej kolumny 0 (długości fal) i rozmiaru siatki N_GRID.
    """
    global COL_INTENSITY, COL_PHASE_RAD, COL_REAL, COL_IMAG

    threshold = parse_power_or_threshold(threshold_arg)
    exponent = np.log10(threshold)

    keep_mask = POS_MASK & (COL_INTENSITY >= threshold)
    zero_mask = ~keep_mask

    n_kept = int(np.sum(keep_mask))
    n_zeroed = int(np.sum(zero_mask))

    if n_kept < 2:
        sys.exit(
            f"[BŁĄD] Powyżej progu {threshold:.2e} (10^{exponent:.1f}) znajduje się mniej niż 2 punkty."
        )

    energy_before = float(np.sum((COL_REAL[POS_MASK] ** 2) + (COL_IMAG[POS_MASK] ** 2)))
    old_fwhm_time = FWHM_TIME_ACTUAL_FS

    wl_min_kept = float(np.min(COL_WAVELENGTH_NM[keep_mask]))
    wl_max_kept = float(np.max(COL_WAVELENGTH_NM[keep_mask]))

    # Twarde wpisanie dokładnych zer 0.0 poza maską impulsu
    COL_REAL[zero_mask] = 0.0
    COL_IMAG[zero_mask] = 0.0
    COL_INTENSITY[zero_mask] = 0.0
    COL_PHASE_RAD[zero_mask] = 0.0

    # Ponowna normalizacja piku Intensity do 1.0 na wypadek drobnych zmian
    max_int = np.max(COL_INTENSITY)
    if max_int > 0:
        COL_INTENSITY /= max_int

    energy_after = float(np.sum((COL_REAL[POS_MASK] ** 2) + (COL_IMAG[POS_MASK] ** 2)))
    retained_pct = (100.0 * energy_after / energy_before) if energy_before > 0 else 0.0

    update_global_diagnostics()

    print("=" * 76)
    print(" [EKSPERYMENTALNE] TWARDE WYZEROWANIE SZUMU DO CZYSTEGO ZERA (0.0)")
    print("=" * 76)
    print(f" Próg odcięcia szumu         : 10^({exponent:.2f}) = {threshold:.2e} (względem piku I = 1.0)")
    print(f" Zachowane punkty impulsu    : {n_kept} z {N_GRID} pkt ({wl_min_kept:.2f} nm .. {wl_max_kept:.2f} nm)")
    print(f" Wyzerowane punkty (0.0)     : {n_zeroed} z {N_GRID} pkt (Re = 0, Im = 0, Int = 0, Phase = 0)")
    print(f" Zachowana energia widmowa   : {retained_pct:.4f}% energii fizycznej")
    print(f" Czas trwania impulsu (FWHM) : {old_fwhm_time:.2f} fs  -->  {FWHM_TIME_ACTUAL_FS:.2f} fs")
    print("=" * 76)


def save_speck_file(output_filepath: str) -> None:
    """
    Zapisuje zmodyfikowane kolumny do pliku tekstowego w 5-kolumnowym formacie FROG/LWE.
    """
    out_path = Path(output_filepath)
    out_matrix = np.column_stack(
        (COL_WAVELENGTH_NM, COL_INTENSITY, COL_PHASE_RAD, COL_REAL, COL_IMAG)
    )
    try:
        np.savetxt(
            out_path,
            out_matrix,
            fmt=["%.7g", "%.7g", "%.7g", "%.7g", "%.7g"],
            delimiter="  ",
        )
        print(f"[ZAPISANO] Zmodyfikowany plik ({N_GRID} pkt) zapisano jako: {out_path.resolve()}")
    except Exception as e:
        sys.exit(f"[BŁĄD] Nie udało się zapisać pliku '{output_filepath}': {e}")


# ==============================================================================
# OBSŁUGA FLAGI --info / -i
# ==============================================================================

def print_file_info() -> None:
    """
    Wyświetla podsumowanie wszystkich obliczonych zmiennych globalnych.
    """
    sep = "=" * 76
    subsep = "-" * 76

    print(sep)
    print(f" DIAGNOSTYKA PLIKU WIDMOWEGO FROG / LWE: {Path(FILE_PATH).name}")
    print(sep)
    print(f" Ścieżka pliku               : {FILE_PATH}")
    print(f" Rozmiar siatki FFT (N_GRID) : {N_GRID} pkt (Potęga 2: {'TAK' if IS_POWER_OF_TWO else 'NIE'})")

    print(subsep)
    print(" 1. PARAMETRY SIATKI FFT (CZĘSTOTLIWOŚĆ I CZAS)")
    print(subsep)
    print(f" Krok częstotliwości (df)    : {DF_THZ:.6f} THz  (d_omega = {DOMEGA_RAD_FS:.6f} rad/fs)")
    print(f" Całkowita szerokość pasma   : {FREQ_SPAN_THZ:.3f} THz  (Nyquist: +/- {NYQUIST_THZ:.3f} THz)")
    print(f" Krok czasowy (dt)           : {DT_FS:.4f} fs")
    print(f" Szerokość okna czasowego    : {TIME_SPAN_FS:.2f} fs")
    print(f" Punkt środkowy FFT (idx {CENTER_GRID_IDX}) : lambda_0 = {CARRIER_WL_NM:.4f} nm ({CARRIER_FREQ_THZ:.3f} THz)")

    print(subsep)
    print(" 2. ZAKRESY DŁUGOŚCI FALI I OSOBLIWOŚĆ (PRZEJŚCIE PRZEZ 0 THz)")
    print(subsep)
    print(f" Punkty fizyczne (lambda > 0): {N_POS_POINTS} wierszy (indeksy {SINGULARITY_IDX + 1 if SINGULARITY_IDX is not None else 0} .. {N_GRID - 1})")
    print(f" Punkty ujemne   (lambda < 0): {N_NEG_POINTS} wierszy (indeksy 0 .. {SINGULARITY_IDX if SINGULARITY_IDX is not None else 0}) -> ignorowane przez LWE")
    if SINGULARITY_IDX is not None:
        wl_neg = COL_WAVELENGTH_NM[SINGULARITY_IDX]
        wl_pos = COL_WAVELENGTH_NM[SINGULARITY_IDX + 1]
        f_neg = FREQ_THZ[SINGULARITY_IDX]
        f_pos = FREQ_THZ[SINGULARITY_IDX + 1]
        print(f" Skok osobliwości (idx {SINGULARITY_IDX}->{SINGULARITY_IDX+1}): "
              f"{wl_neg:.1f} nm ({f_neg:.3f} THz)  -->  +{wl_pos:.1f} nm (+{f_pos:.3f} THz)")
    print(f" Zakres dodatnich fal        : {POS_WL_MIN_NM:.4f} nm  ..  {POS_WL_MAX_NM:.1f} nm")

    print(subsep)
    print(" 3. PARAMETRY WIDMOWE IMPULSU")
    print(subsep)
    print(f" Pik I(lambda) [kolumna 1]   : {PEAK_WL_NM:.4f} nm  ({PEAK_FREQ_THZ:.3f} THz, krok d_lambda = {PEAK_STEP_WL_NM:.3f} nm) [idx {PEAK_IDX}]")
    print(f" Pik |E(nu)|^2 [Re^2 + Im^2] : {FIELD_PEAK_WL_NM:.4f} nm  ({FIELD_PEAK_FREQ_THZ:.3f} THz) [idx {FIELD_PEAK_IDX}]")
    print(f" Środek ciężkości widma (CoM): {COM_WL_NM:.4f} nm  ({COM_FREQ_THZ:.3f} THz, omega_0 = {COM_OMEGA_RAD_FS:.5f} rad/fs)")
    print(f" Szerokość widmowa FWHM (wl) : {FWHM_WL_NM:.2f} nm")
    print(f" Szerokość widmowa FWHM (nu) : {FWHM_FREQ_THZ:.3f} THz  (RMS Bandwidth: {RMS_BW_THZ:.3f} THz)")

    print(subsep)
    print(" 4. FAZA, PHASE BLANKING I SPÓJNOŚĆ DANYCH (JAKOBIAN)")
    print(subsep)
    print(f" Punkty z niezerową kol. 2   : {N_ACTIVE_PHASE_PTS} z {N_GRID} ({100.0 * N_ACTIVE_PHASE_PTS / N_GRID:.1f}%)")
    print(f" Punkty czystego zera (Re=0) : {N_EXACT_ZERO_PTS} z {N_GRID} ({100.0 * N_EXACT_ZERO_PTS / N_GRID:.1f}%)")
    print(f" Szacowany próg Blankingu    : ~{EST_BLANKING_THRESH:.2e} (względem piku I = 1.0)")
    print(f" Okno impulsu (lambda)       : {ACTIVE_WL_MIN_NM:.2f} nm  ..  {ACTIVE_WL_MAX_NM:.2f} nm")
    print(f" Maks. skok fazy w oknie     : {MAX_PHASE_STEP_ACTIVE:.4f} rad/piksel "
          f"({'OK (< 0.5 rad)' if MAX_PHASE_STEP_ACTIVE < 0.5 else 'UWAGA: Duży gradient fazy!'})")
    print(f" Zgodność z Jakobianem 1/l^2 : Maks. odchyłka = {JACOBIAN_REL_ERROR:.2e} "
          f"({'POTWIERDZONA' if JACOBIAN_REL_ERROR < 1e-4 else 'ZMODYFIKOWANA / INNA'})")

    print(subsep)
    print(" 5. SZACOWANE PARAMETRY CZASOWE (Z IFFT OBWIEDNI)")
    print(subsep)
    print(f" Czas trwania Transform-Limit: {FWHM_TIME_TL_FS:.2f} fs (FWHM)")
    print(f" Czas trwania z aktualną fazą: {FWHM_TIME_ACTUAL_FS:.2f} fs (FWHM)")
    if FWHM_TIME_TL_FS > 0:
        print(f" Stosunek poszerzenia impulsu: {FWHM_TIME_ACTUAL_FS / FWHM_TIME_TL_FS:.2f} x Transform-Limit")
    print(sep)


# ==============================================================================
# GŁÓWNY PUNKT WEJŚCIA (CLI)
# ==============================================================================

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="speckutils",
        description="Narzędzie do analizy i modyfikacji plików widmowych FROG (Speck.dat) dla Lightwave Explorer.",
    )
    parser.add_argument(
        "filepath",
        type=str,
        help="Ścieżka do pliku wejściowego (np. ./Speck.dat)",
    )
    parser.add_argument(
        "-i",
        "--info",
        action="store_true",
        help="Wyświetla szczegółowe informacje diagnostyczne i fizyczne o pliku.",
    )
    parser.add_argument(
        "-u",
        "-n",
        "--upsample",
        "--grid-size",
        dest="upsample",
        type=int,
        metavar="N_POINTS",
        help="Zagęszcza (interpoluje) siatkę FFT do wybranej liczby punktów będącej potęgą dwójki "
             "(np. 1024, 2048, 4096).",
    )
    parser.add_argument(
        "-z",
        "--zero-phase",
        type=float,
        metavar="POWER_OR_THRESH",
        help="Zeruje fazę samego impulsu (tworzy impuls Transform-Limited) powyżej zadanego progu "
             "intensywności (np. -3, 3 lub 1e-3 dla progu 10^-3).",
    )
    parser.add_argument(
        "-f",
        "--filter-noise",
        type=float,
        metavar="POWER_OR_THRESH",
        help="Okienkuje (tłumi) szum poniżej podanego progu intensywności względem piku "
             "(np. -3, 3 lub 1e-3) z użyciem gładkiej apodyzacji cosinusowej i miękkiego dna 1e-12.",
    )
    parser.add_argument(
        "-hz",
        "--hard-zero-noise",
        type=float,
        metavar="POWER_OR_THRESH",
        help="[EKSPERYMENTALNE] Twardo zeruje cały szum poniżej podanego progu intensywności "
             "(np. -3, 3 lub 1e-3), wpisując czyste 0.0 do kolumn Intensity, Phase, Re oraz Im.",
    )
    parser.add_argument(
        "--taper-width",
        type=float,
        default=15.0,
        metavar="THZ",
        help="Szerokość gładkiego zbocza cosinusowego (w THz) dla flagi --filter-noise "
             "(domyślnie: 15.0 THz).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="OUT_FILE",
        help="Ścieżka do zapisu zmodyfikowanego pliku wynikowego (np. ./Speck_out.dat).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    # KROK 1: Zawsze wczytujemy plik i obliczamy wszystkie zmienne globalne na starcie
    load_speck_file(args.filepath)

    # KROK 2: Zagęszczenie siatki (jeśli wybrano)
    if args.upsample is not None:
        apply_upsampling(target_n=args.upsample)

    # KROK 3: Zerowanie fazy samego impulsu (jeśli wybrano)
    if args.zero_phase is not None:
        apply_zero_pulse_phase(threshold_arg=args.zero_phase)

    # KROK 4A: Gładkie okienkowanie szumu z miękkim dnem (jeśli wybrano -f)
    if args.filter_noise is not None:
        apply_noise_windowing(
            threshold_arg=args.filter_noise,
            taper_width_thz=args.taper_width,
        )

    # KROK 4B: Eksperymentalne twarde zerowanie szumu do czystego 0.0 (jeśli wybrano -hz)
    if args.hard_zero_noise is not None:
        apply_hard_zero_noise(threshold_arg=args.hard_zero_noise)

    # KROK 5: Wyświetlenie informacji diagnostycznych (po wszystkich modyfikacjach)
    if args.info:
        print_file_info()

    # KROK 6: Zapis zmodyfikowanego pliku, jeśli podano flagę -o / --output
    if args.output:
        save_speck_file(args.output)


if __name__ == "__main__":
    main()
