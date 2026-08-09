import numpy as np
import matplotlib.pyplot as plt

# 1. Definiowanie granic osi X w nanometrach
zakres_dlugosci_fali = [300, 5000]

# 2. Wczytanie pliku Speck.dat
# Format: 0: Długość fali (nm), 1: Intensywność, 2: Faza (ignorowana), 3: Re, 4: Im
dane = np.loadtxt('Base_0mm_Speck.dat')

dlugosc_fali = dane[:, 0]
intensywnosc = dane[:, 1]
czesc_rzeczywista = dane[:, 3]
czesc_urojona = dane[:, 4]

# 3. Obliczenie fazy z części rzeczywistej i urojonej
# Użycie np.unwrap zapobiega skokom fazy o 2*pi
faza_obliczona = np.unwrap(np.arctan2(czesc_urojona, czesc_rzeczywista))

# 4. Normalizacja intensywności
if np.max(intensywnosc) > 0:
    intensywnosc = intensywnosc / np.max(intensywnosc)

# 5. Inicjalizacja wykresu
fig, os_intensywnosci = plt.subplots(figsize=(10, 6))

# Konfiguracja głównej osi Y (Intensywność)
kolor_int = 'tab:blue'
os_intensywnosci.set_xlabel('Długość fali (nm)', fontsize=12)
os_intensywnosci.set_ylabel('Znormalizowana Intensywność', color=kolor_int, fontsize=12)
wykres_int = os_intensywnosci.plot(dlugosc_fali, intensywnosc, color=kolor_int, label='Intensywność', linewidth=2)
os_intensywnosci.tick_params(axis='y', labelcolor=kolor_int)

# Ustawienie limitów osi X i Y
os_intensywnosci.set_xlim(zakres_dlugosci_fali)
os_intensywnosci.set_ylim([0, 1.05])
os_intensywnosci.grid(True, linestyle='--', alpha=0.6)

# Konfiguracja drugiej osi Y (Faza) współdzielącej oś X
os_fazy = os_intensywnosci.twinx()
kolor_fazy = 'tab:red'
os_fazy.set_ylabel('Faza (rad)', color=kolor_fazy, fontsize=12)

# Rysowanie fazy tylko w zadanym zakresie poprawia czytelność autoskalowania
indeksy_w_zakresie = (dlugosc_fali >= zakres_dlugosci_fali[0]) & (dlugosc_fali <= zakres_dlugosci_fali[1])
wykres_faza = os_fazy.plot(dlugosc_fali[indeksy_w_zakresie], faza_obliczona[indeksy_w_zakresie], 
                           color=kolor_fazy, label='Faza (z Re/Im)', linestyle='--', linewidth=2)
os_fazy.tick_params(axis='y', labelcolor=kolor_fazy)

# Tworzenie wspólnej legendy
linie = wykres_int + wykres_faza
etykiety = [l.get_label() for l in linie]
os_intensywnosci.legend(linie, etykiety, loc='upper right')

plt.title('Widmo impulsu: Znormalizowana intensywność i faza z Re/Im', fontsize=14)
plt.tight_layout()
plt.savefig("Base_0mm_Speck.png")