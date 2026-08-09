import matplotlib.pyplot as plt
import numpy as np

class ExperimentVisualizer:
    def __init__(self, save_path="/workspace/training_progress.png"):
        self.save_path = save_path

    def update_plot(self, iteration, max_iterations, history, modified_pulse_data, original_pulse_data, best_result_data, eval_f_min, eval_f_max):
        fig, axs = plt.subplots(2, 2, figsize=(16, 10))
        
        # --- Wykres 1 (Lewy Góra): Zmodyfikowany impuls wejściowy na SLM (Speck.dat) ---
        ax_tl = axs[0, 0]
        if modified_pulse_data is not None:
            wl = modified_pulse_data[:, 0]
            re_mod = modified_pulse_data[:, 3]
            im_mod = modified_pulse_data[:, 4]
            
            orig_int = re_mod**2 + im_mod**2
            phase = np.angle(re_mod + 1j * im_mod)
            
            prog_szumu = 0.01 * np.max(orig_int)
            phase[orig_int < prog_szumu] = np.nan
            
            if np.max(orig_int) > 0:
                orig_int = orig_int / np.max(orig_int)
                
            kolor_int_mod = 'black'
            ax_tl.set_xlabel("Długość fali (nm)", fontsize=10)
            ax_tl.set_ylabel("Znormalizowana Intensywność", color=kolor_int_mod, fontsize=10)
            wykres_int_mod = ax_tl.plot(wl, orig_int, color=kolor_int_mod, label="Intensywność", linewidth=2)
            ax_tl.tick_params(axis='y', labelcolor=kolor_int_mod)
            ax_tl.set_xlim([500, 4000])
            ax_tl.set_ylim([0, 1.05])
            ax_tl.grid(True, linestyle='--', alpha=0.6)
            
            ax_tl_phase = ax_tl.twinx()
            kolor_fazy_mod = 'blue'
            ax_tl_phase.set_ylabel("Faza (rad)", color=kolor_fazy_mod, fontsize=10)
            
            indeksy_w_zakresie_mod = (wl >= 500) & (wl <= 4000)
            wykres_faza_mod = ax_tl_phase.step(wl[indeksy_w_zakresie_mod], phase[indeksy_w_zakresie_mod], 
                                               color=kolor_fazy_mod, lw=2, label="Faza po modulacji")
            ax_tl_phase.tick_params(axis='y', labelcolor=kolor_fazy_mod)
            
            linie_mod = wykres_int_mod + wykres_faza_mod
            etykiety_mod = [l.get_label() for l in linie_mod]
            ax_tl.legend(linie_mod, etykiety_mod, loc='upper right')
            
            ax_tl.set_title(f"Iteracja {iteration}/{max_iterations} | Zmodyfikowany impuls (Speck.dat)")

        # --- Wykres 2 (Lewy Dół): Historia zbieżności Fitness ---
        ax_bl = axs[1, 0]
        best_fitness_history = []
        for gen_data in history:
            best_in_gen = max([e["fitness"] for e in gen_data["evaluations"]])
            best_fitness_history.append(best_in_gen)
            
        ax_bl.plot(best_fitness_history, marker='o', linestyle='-', color='firebrick')
        best_fit_val = best_result_data["fitness"] if best_result_data else 0.0
        ax_bl.set_title(f"Historia Zbieżności | Rekord Całki: {best_fit_val:.4e}")
        ax_bl.set_xlabel("Iteracja")
        ax_bl.set_ylabel("Wartość Funkcji Celu")
        ax_bl.grid(True)

        # --- Wykres 3 (Prawy Góra): Oryginalny impuls wejściowy (Speck.dat) ---
        ax_tr = axs[0, 1]
        if original_pulse_data is not None:
            zakres_dlugosci_fali = [500, 4000]
            dlugosc_fali = original_pulse_data[:, 0]
            
            re_comp = original_pulse_data[:, 3]
            im_comp = original_pulse_data[:, 4]
            intensywnosc = re_comp**2 + im_comp**2
            faza = np.angle(re_comp + 1j * im_comp)
            
            prog_szumu = 0.01 * np.max(intensywnosc)
            faza[intensywnosc < prog_szumu] = np.nan
            
            if np.max(intensywnosc) > 0:
                intensywnosc = intensywnosc / np.max(intensywnosc)
                
            kolor_int = 'tab:blue'
            ax_tr.set_xlabel('Długość fali (nm)', fontsize=10)
            ax_tr.set_ylabel('Znormalizowana Intensywność', color=kolor_int, fontsize=10)
            wykres_int = ax_tr.plot(dlugosc_fali, intensywnosc, color=kolor_int, label='Intensywność', linewidth=2)
            ax_tr.tick_params(axis='y', labelcolor=kolor_int)
            ax_tr.set_xlim(zakres_dlugosci_fali)
            ax_tr.set_ylim([0, 1.05])
            ax_tr.grid(True, linestyle='--', alpha=0.6)
            
            os_fazy = ax_tr.twinx()
            kolor_fazy = 'tab:red'
            os_fazy.set_ylabel('Faza (rad)', color=kolor_fazy, fontsize=10)
            
            indeksy_w_zakresie = (dlugosc_fali >= zakres_dlugosci_fali[0]) & (dlugosc_fali <= zakres_dlugosci_fali[1])
            wykres_faza = os_fazy.plot(dlugosc_fali[indeksy_w_zakresie], faza[indeksy_w_zakresie], 
                                       color=kolor_fazy, label='Faza', linestyle='--', linewidth=2)
            os_fazy.tick_params(axis='y', labelcolor=kolor_fazy)
            
            c_nm_s = 299792458.0 * 1e9
            lambda_max_integ = c_nm_s / eval_f_min if eval_f_min > 0 else np.inf
            lambda_min_integ = c_nm_s / eval_f_max if eval_f_max > 0 else 0
            
            maska_calkowania = (dlugosc_fali >= lambda_min_integ) & (dlugosc_fali <= lambda_max_integ)
            ax_tr.fill_between(dlugosc_fali, intensywnosc, where=maska_calkowania, color='gray', alpha=0.3, label='Obszar całkowania')
            
            linie = wykres_int + wykres_faza
            etykiety = [l.get_label() for l in linie]
            ax_tr.legend(linie, etykiety, loc='upper right')
            
            ax_tr.set_title('Oryginalny impuls wejściowy (Speck.dat)')

        def hz_to_nm(freq_array):
            wl_array = np.zeros_like(freq_array)
            valid = freq_array > 0
            wl_array[valid] = (299792458.0 / freq_array[valid]) * 1e9
            return wl_array

        # --- Wykres 4 (Prawy Dół): Najlepszy globalny wynik Spektrum z LWE ---
        ax_br = axs[1, 1]
        if best_result_data and best_result_data["freq"] is not None:
            f = best_result_data["freq"]
            s = best_result_data["spec"]
            p = best_result_data["phase"]
            
            wl_nm = hz_to_nm(f)
            
            p = p.astype(float)
            if np.max(s) > 0:
                p[s < 0.01 * np.max(s)] = np.nan
            
            ax_br.plot(wl_nm, s, color='black')
            mask = (f >= eval_f_min) & (f <= eval_f_max)
            ax_br.fill_between(wl_nm, s, where=mask, color='firebrick', alpha=0.4, label="Zakres Całkowania")
            ax_br.set_xlabel("Długość fali (nm)")
            ax_br.set_ylabel("Intensywność wyjściowa", color='black')
            
            if np.any(mask):
                valid_wl = wl_nm[mask]
                span = valid_wl.max() - valid_wl.min()
                ax_br.set_xlim(max(0, valid_wl.min() - 0.2 * span), valid_wl.max() + 0.2 * span)
            else:
                ax_br.set_xlim(500, 4000)
                
            ax_br_phase = ax_br.twinx()
            ax_br_phase.plot(wl_nm, p, color='green', alpha=0.5, linestyle='--')
            ax_br_phase.set_ylabel("Faza (rad)", color='green')
            ax_br.set_title("Najlepszy Wygenerowany Impuls z LWE (Globalnie)")

        plt.tight_layout()
        plt.savefig(self.save_path, dpi=100)
        plt.close(fig)