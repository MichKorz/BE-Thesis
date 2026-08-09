import json
import numpy as np
import matplotlib.pyplot as plt

# --- KONFIGURACJA ---
JSON_PATH = "configs/run_001.json"  # Ścieżka do Twojego pliku JSON z historią runu
OUTPUT_IMAGE = "slm_mask_1920x1080.png"

# Parametry fizyczne matrycy SLM oraz okna optycznego (muszą pasować do tych z JSON/modulatora)
SLM_WIDTH = 1920   # Liczba pikseli w poziomie dla Pluto-2.1
SLM_HEIGHT = 1080  # Liczba pikseli w pionie
LAMBDA_MIN = 1070.0 # nm
LAMBDA_MAX = 3700.0 # nm

def load_history(json_path):
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data.get("history", [])

def find_best_vector(history):
    """Szuka wektora o najwyższym fitness w całej historii treningu."""
    best_vector = None
    max_fitness = float('-inf')
    
    for gen in history:
        for ev in gen.get("evaluations", []):
            if ev["fitness"] > max_fitness:
                max_fitness = ev["fitness"]
                best_vector = ev["vector"]
                
    return best_vector, max_fitness

def main():
    history = load_history(JSON_PATH)
    if not history:
        print("[Błąd] Brak historii w podanym pliku JSON.")
        return
        
    best_vector, best_fitness = find_best_vector(history)
    if best_vector is None:
        print("[Błąd] Nie udało się odnaleźć żadnego wektora w historii.")
        return
        
    print(f"Znaleziono najlepszy wektor o fitness: {best_fitness:.4e}")
    num_stripes = len(best_vector)
    
    # 1. Tworzenie 1D profilu fazy w zależności od długości fali (tak jak w SLMModulatorze)
    # Tworzymy siatkę odpowiadającą szerokości poziomej SLM (1920 pikseli mapujemy na zakres lambda)
    wavelengths = np.linspace(LAMBDA_MIN, LAMBDA_MAX, SLM_WIDTH)
    phase_profile = np.zeros_like(wavelengths)
    
    bins = np.linspace(LAMBDA_MIN, LAMBDA_MAX, num_stripes + 1)
    for i in range(num_stripes):
        mask = (wavelengths >= bins[i]) & (wavelengths <= bins[i+1])
        if np.any(mask):
            phase_profile[mask] = best_vector[i]
            
    # 2. Przeliczenie fazy [0, 2*pi] na skale szarości [0, 255]
    # Normalizacja do przedziału [0, 1] a następnie przemnożenie przez 255
    normalized_phase = np.mod(phase_profile, 2 * np.pi) / (2 * np.pi)
    gray_levels = np.uint8(normalized_phase * 255)
    
    # 3. Rozciągnięcie profilu 1D na całą wysokość matrycy (1920x1080)
    # Paski na SLM są zazwyczaj pionowymi kolumnami rozciągniętymi wzdłuż osi Y
    slm_image_data = np.tile(gray_levels, (SLM_HEIGHT, 1))
    
    # 4. Zapis do pliku PNG
    plt.figure(figsize=(10, 3))
    plt.imshow(slm_image_data, cmap='gray', vmin=0, vmax=255, aspect='auto')
    plt.colorbar(label='Poziom szarości (0-255)')
    plt.title(f"Wygenerowana maska SLM (Pasków: {num_stripes})")
    plt.xlabel("Piksele poziome SLM (0 - 1920)")
    plt.ylabel("Pion")
    plt.tight_layout()
    plt.savefig("mask_preview.png", dpi=150)
    plt.close()
    
    # Zapis właściwej, czystej bitmapy 1920x1080 bez ramek i osi dla urządzenia
    plt.imsave(OUTPUT_IMAGE, slm_image_data, cmap='gray', vmin=0, vmax=255)
    print(f"Sukces! Zapisano maskę gotową dla SLM do pliku: {OUTPUT_IMAGE}")
    print(f"Zapisano również podgląd z osiami do: mask_preview.png")

if __name__ == "__main__":
    main()