import json
import math
import numpy as np
import matplotlib.pyplot as plt

# --- KONFIGURACJA ---
JSON_PATH = "/workspace/Runs/Bluebird/TaylorFitting/ZGP_DFG-BB-TF-30_agents-640_stripes.json"
OUTPUT_IMAGE = "slm_mask_taylor_1920x1080.png"
PREVIEW_IMAGE = "mask_preview_taylor.png"
CONTINUOUS_PHASE_IMAGE = "continuous_phase_taylor.png"

C_NM_FS = 299.792458
SLM_WIDTH = 1920
SLM_HEIGHT = 1080

def load_data(json_path):
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
        
    mod_kwargs = data.get("config", {}).get("modulator", {}).get("kwargs", {})
    lambda_min = mod_kwargs.get("lambda_min", 1000.0)
    lambda_max = mod_kwargs.get("lambda_max", 3200.0)
    num_stripes = mod_kwargs.get("num_stripes", 80)
    lambda_0 = mod_kwargs.get("lambda_0", 2500.0)
    
    history = data.get("history", [])
    best_vector = None
    max_fitness = float('-inf')
    
    for gen in history:
        for ev in gen.get("evaluations", []):
            if ev["fitness"] > max_fitness:
                max_fitness = ev["fitness"]
                best_vector = ev["vector"]
                
    return best_vector, max_fitness, lambda_min, lambda_max, num_stripes, lambda_0

def main():
    best_vector, best_fitness, lambda_min, lambda_max, num_stripes, lambda_0 = load_data(JSON_PATH)
    
    if best_vector is None:
        print("[Błąd] Nie odnaleziono wektora w pliku JSON.")
        return
        
    print(f"Parametry wczytane z konfiguracji:")
    print(f"  Zakres fali: {lambda_min} nm - {lambda_max} nm")
    print(f"  Długość centralna (lambda_0): {lambda_0} nm")
    print(f"  Liczba pasków: {num_stripes}")
    print(f"Najlepszy fitness: {best_fitness:.4e}")
    
    omega_0 = 2 * np.pi * C_NM_FS / lambda_0
    wavelengths = np.linspace(lambda_min, lambda_max, SLM_WIDTH)
    
    # 1. Wykres i obliczenia dla ciągłej fazy Taylora (przed paskowaniem i modulo)
    omega_continuous = 2 * np.pi * C_NM_FS / wavelengths
    delta_omega_continuous = omega_continuous - omega_0
    
    continuous_phase = np.zeros_like(wavelengths)
    for order_idx, coeff in enumerate(best_vector):
        n = order_idx + 2
        continuous_phase += (coeff / math.factorial(n)) * (delta_omega_continuous ** n)
        
    plt.figure(figsize=(10, 4))
    plt.plot(wavelengths, continuous_phase, color='blue', linewidth=2)
    plt.title("Ciągła faza z szeregu Taylora (przed dyskretyzacją i modulo)")
    plt.xlabel("Długość fali (nm)")
    plt.ylabel("Faza (rad)")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(CONTINUOUS_PHASE_IMAGE, dpi=150)
    plt.close()
    
    # 2. Obliczenie wartości Taylora + Phase Wrapping dla każdego paska
    bins = np.linspace(lambda_min, lambda_max, num_stripes + 1)
    stripe_phases = np.zeros(num_stripes)
    
    for i in range(num_stripes):
        lambda_center = (bins[i] + bins[i+1]) / 2.0
        omega_center = 2 * np.pi * C_NM_FS / lambda_center
        delta_omega = omega_center - omega_0
        
        taylor_phase = 0.0
        for order_idx, coeff in enumerate(best_vector):
            n = order_idx + 2
            taylor_phase += (coeff / math.factorial(n)) * (delta_omega ** n)
            
        stripe_phases[i] = taylor_phase % (2 * np.pi)
        
    # 3. Mapowanie wyliczonych faz na fizyczne 1920 pikseli monitora SLM
    phase_profile = np.zeros_like(wavelengths)
    for i in range(num_stripes):
        mask = (wavelengths >= bins[i]) & (wavelengths <= bins[i+1])
        phase_profile[mask] = stripe_phases[i]
        
    # 4. Przeliczenie fazy [0, 2*pi] na poziomy szarości [0, 255] i rozciągnięcie matrycy na 1080px
    normalized_phase = phase_profile / (2 * np.pi)
    gray_levels = np.uint8(normalized_phase * 255)
    slm_image_data = np.tile(gray_levels, (SLM_HEIGHT, 1))
    
    # 5. Generacja podglądu z naniesionymi osiami
    plt.figure(figsize=(12, 4))
    plt.imshow(slm_image_data, cmap='gray', vmin=0, vmax=255, aspect='auto')
    plt.colorbar(label='Poziom szarości (0-255)')
    plt.title(f"Maska Taylor SLM | Zawinięta faza (Modulo 2π)")
    
    num_ticks = 8
    tick_locs = np.linspace(0, SLM_WIDTH - 1, num_ticks)
    tick_labels = np.round(np.linspace(lambda_min, lambda_max, num_ticks), 1)
    plt.xticks(tick_locs, tick_labels)
    plt.xlabel("Długość fali zmapowana na piksele (nm)")
    plt.ylabel("Pion matrycy (px)")
    plt.tight_layout()
    plt.savefig(PREVIEW_IMAGE, dpi=150)
    plt.close()
    
    # 6. Generacja ostatecznej maski dla urządzenia
    plt.imsave(OUTPUT_IMAGE, slm_image_data, cmap='gray', vmin=0, vmax=255)
    print(f"Sukces! Zapisano gotową maskę do pliku: {OUTPUT_IMAGE}")
    print(f"Zapisano również plik podglądowy: {PREVIEW_IMAGE}")
    print(f"Zapisano wykres ciągłej fazy: {CONTINUOUS_PHASE_IMAGE}")

if __name__ == "__main__":
    main()