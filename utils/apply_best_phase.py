import argparse
import json
import sys
import os

# --- KONFIGURACJA ŚCIEŻKI ---
# Podaj "na twardo" ścieżkę do katalogu, w którym leży phase_modulator.py
PHASE_MODULATOR_DIR = "../modulators/"

# Dodanie ścieżki na początek listy systemowej
if os.path.isdir(PHASE_MODULATOR_DIR):
    sys.path.insert(0, PHASE_MODULATOR_DIR)
else:
    print(f"Błąd: Zdefiniowany katalog '{PHASE_MODULATOR_DIR}' nie istnieje.")
    sys.exit(1)

# Import klas po wstrzyknięciu ścieżki
try:
    from phase_modulator import SLMModulator, TaylorSLMModulator
except ImportError as e:
    print(f"Błąd importu: Nie można znaleźć pliku phase_modulator.py w katalogu {PHASE_MODULATOR_DIR}.")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Aplikuje optymalną fazę do pliku Speck.dat na podstawie JSON.")
    parser.add_argument("json_path", help="Ścieżka do pliku konfiguracyjnego JSON z wynikami.")
    parser.add_argument("speck_in", help="Ścieżka do wejściowego pliku Speck.dat.")
    parser.add_argument("--speck_out", default="Speck_modified.dat", help="Ścieżka zapisu (domyślnie: Speck_modified.dat)")
    
    args = parser.parse_args()

    with open(args.json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    mod_config = data.get("config", {}).get("modulator", {})
    mod_class = mod_config.get("class")
    mod_kwargs = mod_config.get("kwargs", {})
    best_vector = data.get("state", {}).get("best_vector")

    if not mod_class or not best_vector:
        print("Błąd: Plik JSON nie zawiera definicji modulatora lub najlepszego wektora.")
        sys.exit(1)

    mod_kwargs["input_path"] = args.speck_in
    mod_kwargs["output_path"] = args.speck_out

    if mod_class == "TaylorSLMModulator":
        modulator = TaylorSLMModulator(**mod_kwargs)
    elif mod_class == "SLMModulator":
        modulator = SLMModulator(**mod_kwargs)
    else:
        print(f"Błąd: Nieobsługiwana klasa modulatora: {mod_class}")
        sys.exit(1)

    print(f"Użyto modulatora: {mod_class}")
    print(f"Aplikowanie wektora: {best_vector}")
    modulator.apply_phase(best_vector)
    print(f"Zakończono sukcesem. Zmodyfikowany plik zapisano jako: {args.speck_out}")


if __name__ == "__main__":
    main()