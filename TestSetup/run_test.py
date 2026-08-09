import os
import shutil
import subprocess
import json

# Konfiguracja ścieżek
RAMDISK_DIR = "/mnt/ramdisk"
TEST_SETUP_DIR = "/workspace/TestSetup"
SHG_BASE_DIR = "/workspace/TestSetup/SHG_BASE" 
MODULATOR_SCRIPT = "/workspace/utils/phase_modulator.py"
LWE_EXE = "/workspace/LightwaveExplorer/build/LightwaveExplorer"
LWE_WORKING_DIR = "/workspace/LightwaveExplorer/build"
PYTHON_BIN = "/opt/venv/bin/python"

# Pliki wejściowe
INPUT_CONFIG = os.path.join(SHG_BASE_DIR, "SHG_BASE.txt")
INPUT_PULSE = os.path.join(SHG_BASE_DIR, "SHG_BASE_pulse1.dat")

# Pliki na RAM dysku
RAMDISK_CONFIG = os.path.join(RAMDISK_DIR, "SHG_BASE.txt")
RAMDISK_PULSE = os.path.join(RAMDISK_DIR, "SHG_BASE_pulse1.dat")
MODIFIED_PULSE = os.path.join(RAMDISK_DIR, "SHG_MODIFIED_pulse1.dat")
OUTPUT_BASE = os.path.join(RAMDISK_DIR, "SHG_BASE_results")
OUTPUT_GDD = os.path.join(RAMDISK_DIR, "SHG_GDD250_results")

def update_lwe_paths(config_path, new_pulse_path, new_output_path):
    with open(config_path, 'r') as f:
        lines = f.readlines()
    
    with open(config_path, 'w') as f:
        for line in lines:
            if line.startswith("Field 1 file path:"):
                f.write(f"Field 1 file path: {new_pulse_path}\n")
            elif line.startswith("Output base path:"):
                f.write(f"Output base path: {new_output_path}\n")
            else:
                f.write(line)

def run_lwe_simulation(config_path, timeout_sec=15):
    try:
        subprocess.run(
            [LWE_EXE, config_path],
            cwd=LWE_WORKING_DIR,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout_sec,
            check=True
        )
        print(" -> Symulacja zakończona pomyślnie.")
    except subprocess.TimeoutExpired:
        print(" -> [BŁĄD] Symulacja przekroczyła czas i została zabita (Timeout).")
    except subprocess.CalledProcessError as e:
        print(f" -> [BŁĄD] Symulacja LWE uległa awarii. Kod błędu: {e.returncode}")

def apply_dispersion(input_pulse, output_pulse, dispersion_values):
    payload = {
        "input_path": input_pulse,
        "output_path": output_pulse, 
        "coefficients": dispersion_values 
    }
    
    result = subprocess.run(
        [PYTHON_BIN, MODULATOR_SCRIPT],
        input=json.dumps(payload),
        text=True,
        capture_output=True,
        check=False
    )
    
    if result.returncode != 0:
        print("\n--- BŁĄD KRYTYCZNY MODULATORA FAZY ---")
        print(f"STDOUT: {result.stdout.strip()}")
        print(f"STDERR: {result.stderr.strip()}")
        print("--------------------------------------\n")
        raise Exception("Zatrzymano test z powodu awarii phase_modulator.py")
    
    try:
        response = json.loads(result.stdout)
        if response.get("status") == "success":
            print(f" -> Zaaplikowano GDD = {dispersion_values[0]}, TOD = {dispersion_values[1]}. Max przesunięcie: {response.get('max_phase_shift'):.2f} rad")
        else:
            print(f"Błąd zwrócony przez skrypt: {response.get('message')}")
    except json.JSONDecodeError:
        print("Nie udało się sparsować odpowiedzi JSON od modulatora.")

if __name__ == "__main__":
    print("1. Inicjalizacja i kopiowanie plików na RAM dysk...")
    shutil.copy(INPUT_CONFIG, RAMDISK_CONFIG)
    shutil.copy(INPUT_PULSE, RAMDISK_PULSE)
    
    print("2. Uruchomienie symulacji bazowej (GDD = 0)...")
    update_lwe_paths(RAMDISK_CONFIG, RAMDISK_PULSE, OUTPUT_BASE)
    run_lwe_simulation(RAMDISK_CONFIG)
    
    print("3. Kopiowanie i modyfikacja fazy impulsu...")
    shutil.copy(RAMDISK_PULSE, MODIFIED_PULSE)
    apply_dispersion(RAMDISK_PULSE, MODIFIED_PULSE, [50.0, 250.0])
    
    print("4. Uruchomienie symulacji zmodyfikowanej (GDD = 50, TOD = 250)...")
    update_lwe_paths(RAMDISK_CONFIG, MODIFIED_PULSE, OUTPUT_GDD)
    run_lwe_simulation(RAMDISK_CONFIG)
    
    print("Zakończono sekwencję testową.")