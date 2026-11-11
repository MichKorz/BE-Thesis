import matplotlib.pyplot as plt
import numpy as np

class ExperimentVisualizer:
    def __init__(self, save_path="/workspace/training_progress.png"):
        self.save_path = save_path

    def update_plot(self, iteration, max_iterations, history, original_pulse_data, best_modified_pulse_data, base_out_pulse_data, best_out_pulse_data, best_fitness_val, eval_f_min, eval_f_max):
        fig, axs = plt.subplots(2, 2, figsize=(16, 10))
        
        c_nm_s = 299792458.0 * 1e9
        lambda_max_integ = c_nm_s / eval_f_min if eval_f_min > 0 else np.inf
        lambda_min_integ = c_nm_s / eval_f_max if eval_f_max > 0 else 0

        # --- Wykres 1 (Lewy Góra): Symulacja Bazowa (base_out_pulse_data) ---
        ax_tl = axs[0, 0]
        if base_out_pulse_data is not None:
            dlugosc_fali_base = base_out_pulse_data[:, 0]
            re_base = base_out_pulse_data[:, 3]
            im_base = base_out_pulse_data[:, 4]
            
            intensywnosc_base = re_base**2 + im_base**2
            faza_base = np.angle(re_base + 1j * im_base)
            
            prog_szumu_base = 0.01 * np.max(intensywnosc_base)
            faza_base[intensywnosc_base < prog_szumu_base] = np.nan
            
            if np.max(intensywnosc_base) > 0:
                intensywnosc_base = intensywnosc_base / np.max(intensywnosc_base)
                
            kolor_int = 'black'
            ax_tl.plot(dlugosc_fali_base, intensywnosc_base, color=kolor_int, label='Intensywność', linewidth=2)
            
            # --- ZMIANA: Obliczanie zacieniowanego pola ---
            maska_calkowania_base = (dlugosc_fali_base >= lambda_min_integ) & (dlugosc_fali_base <= lambda_max_integ)
            
            pole_base = 0.0
            if np.any(maska_calkowania_base):
                # np.abs zabezpiecza przed ujemnym wynikiem jeśli oś X jest malejąca
                pole_base = np.abs(np.trapezoid(intensywnosc_base[maska_calkowania_base], dlugosc_fali_base[maska_calkowania_base]))
                
            label_base = f"Zakres Całkowania (Pole: {pole_base:.2f})"
            ax_tl.fill_between(dlugosc_fali_base, intensywnosc_base, where=maska_calkowania_base, color='gray', alpha=0.4, label=label_base)
            
            # Włączamy legendę żeby wyświetliła obliczone pole
            ax_tl.legend(loc='upper right', fontsize=10)
            # ----------------------------------------------
            
            ax_tl.set_xlabel("Długość fali (nm)", fontsize=10)
            ax_tl.set_ylabel("Znormalizowana Intensywność", color=kolor_int, fontsize=10)
            ax_tl.tick_params(axis='y', labelcolor=kolor_int)
            ax_tl.set_ylim([0, 1.05])
            ax_tl.set_xlim([500, 20000])
            
            ax_tl_phase = ax_tl.twinx()
            kolor_fazy = 'tab:red'
            ax_tl_phase.plot(dlugosc_fali_base, faza_base, color=kolor_fazy, alpha=0.5, linestyle='--', linewidth=2, label='Faza')
            ax_tl_phase.set_ylabel("Faza (rad)", color=kolor_fazy, fontsize=10)
            ax_tl_phase.tick_params(axis='y', labelcolor=kolor_fazy)
            
            ax_tl.set_title("Symulacja Bazowa (przed optymalizacją)")

        # --- Wykres 3 (Prawy Góra): Najlepszy Wynik Symulacji ---
        ax_tr = axs[0, 1]
        if best_out_pulse_data is not None:
            dlugosc_fali_out = best_out_pulse_data[:, 0]
            re_out = best_out_pulse_data[:, 3]
            im_out = best_out_pulse_data[:, 4]
            
            intensywnosc_out = re_out**2 + im_out**2
            faza_out = np.angle(re_out + 1j * im_out)
            
            prog_szumu_out = 0.01 * np.max(intensywnosc_out)
            faza_out[intensywnosc_out < prog_szumu_out] = np.nan
            
            if np.max(intensywnosc_out) > 0:
                intensywnosc_out = intensywnosc_out / np.max(intensywnosc_out)
                
            kolor_int_out = 'black'
            ax_tr.plot(dlugosc_fali_out, intensywnosc_out, color=kolor_int_out, label='Intensywność', linewidth=2)
            
            # --- ZMIANA: Obliczanie zacieniowanego pola ---
            maska_calkowania_out = (dlugosc_fali_out >= lambda_min_integ) & (dlugosc_fali_out <= lambda_max_integ)
            
            pole_out = 0.0
            if np.any(maska_calkowania_out):
                pole_out = pole_out = np.abs(np.trapezoid(intensywnosc_out[maska_calkowania_out], dlugosc_fali_out[maska_calkowania_out]))
                
            label_out = f"Zakres Całkowania (Pole: {pole_out:.2f})"
            ax_tr.fill_between(dlugosc_fali_out, intensywnosc_out, where=maska_calkowania_out, color='firebrick', alpha=0.4, label=label_out)
            
            # Włączamy legendę żeby wyświetliła obliczone pole
            ax_tr.legend(loc='upper right', fontsize=10)
            # ----------------------------------------------
            
            ax_tr.set_xlabel("Długość fali (nm)", fontsize=10)
            ax_tr.set_ylabel("Znormalizowana Intensywność", color=kolor_int_out, fontsize=10)
            ax_tr.tick_params(axis='y', labelcolor=kolor_int_out)
            ax_tr.set_ylim([0, 1.05])
            ax_tr.set_xlim([500, 20000])
            
            ax_tr_phase = ax_tr.twinx()
            kolor_fazy_out = 'green'
            ax_tr_phase.plot(dlugosc_fali_out, faza_out, color=kolor_fazy_out, alpha=0.5, linestyle='--', linewidth=2, label='Faza')
            ax_tr_phase.set_ylabel("Faza (rad)", color=kolor_fazy_out, fontsize=10)
            ax_tr_phase.tick_params(axis='y', labelcolor=kolor_fazy_out)
            
            ax_tr.set_title(f"Najlepszy Wygenerowany Impuls z LWE (Iteracja {iteration}/{max_iterations})")

        # --- Wykres 2 (Lewy Dół): Historia zbieżności Fitness ---
        ax_bl = axs[1, 0]
        best_fitness_history = []
        for gen_data in history:
            best_in_gen = max([e["fitness"] for e in gen_data["evaluations"]])
            best_fitness_history.append(best_in_gen)
            
        ax_bl.plot(best_fitness_history, marker='o', linestyle='-', color='firebrick')
        ax_bl.set_title(f"Historia Zbieżności | Rekord: {best_fitness_val:.4e}")
        ax_bl.set_xlabel("Iteracja")
        ax_bl.set_ylabel("Wartość Funkcji Celu")
        ax_bl.grid(True)

        # --- Wykres 4 (Prawy Dół): Najlepsza Maska (Wykres fazy) ---
        ax_br = axs[1, 1]
        if original_pulse_data is not None and best_modified_pulse_data is not None:
            wl = original_pulse_data[:, 0]
            
            re_orig = original_pulse_data[:, 3]
            im_orig = original_pulse_data[:, 4]
            orig_complex = re_orig + 1j * im_orig
            
            re_mod = best_modified_pulse_data[:, 3]
            im_mod = best_modified_pulse_data[:, 4]
            mod_complex = re_mod + 1j * im_mod
            
            intensywnosc = re_orig**2 + im_orig**2
            prog_szumu = 0.01 * np.max(intensywnosc)
            
            # Odjęcie fazy bazowej od zmodyfikowanej w zespolonej dziedzinie
            mask_phase = np.angle(mod_complex * np.conj(orig_complex))
            # Ograniczenie wartości do przedziału sprzętowego [0, 2π]
            mask_phase = np.mod(mask_phase, 2 * np.pi)
            
            mask_phase[intensywnosc < prog_szumu] = np.nan
            
            ax_br.step(wl, mask_phase, color='blue', linewidth=2, label='Maska fazowa (rad)')
            ax_br.set_xlabel("Długość fali (nm)", fontsize=10)
            ax_br.set_ylabel("Faza (rad)", color='blue', fontsize=10)
            ax_br.tick_params(axis='y', labelcolor='blue')
            ax_br.set_xlim([500, 4000]) 
            ax_br.set_ylim([-0.2, 2 * np.pi + 0.2])
            ax_br.grid(True, linestyle='--', alpha=0.6)
            
            ax_br.set_title("Najlepsza Maska Fazy na SLM (Modulo 2π)")

        plt.tight_layout()
        plt.savefig(self.save_path, dpi=100)
        plt.close(fig)