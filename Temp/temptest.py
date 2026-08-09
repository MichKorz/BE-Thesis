import numpy as np
import LightwaveExplorer as lwe

symulacja = lwe.load("./Base_0mm.zip", True)
srodek_x = symulacja.Nspace // 2

if len(symulacja.Ext_x.shape) == 2:
    E_t = symulacja.Ext_x[:, srodek_x]
else:
    E_t = symulacja.Ext_x[:, srodek_x].flatten()[:symulacja.Ntime]

# 1. USUNIĘCIE FAZY LINIOWEJ: Znalezienie szczytu i przesunięcie do indeksu 0
indeks_szczytu = np.argmax(np.abs(E_t))
E_t_wycentrowane = np.roll(E_t, -indeks_szczytu)

# 2. Transformata na wycentrowanym impulsie
E_f = np.fft.fft(E_t_wycentrowane)

# 3. Odcięcie częstotliwości 0 Hz
E_f_dodatnie = E_f[1:symulacja.Nfreq]
czestotliwosci = symulacja.frequencyVectorSpectrum[1:]

# 4. Poprawna normalizacja zespolona
max_amplituda = np.max(np.abs(E_f_dodatnie))
if max_amplituda > 0:
    E_f_norm = E_f_dodatnie / max_amplituda
else:
    E_f_norm = E_f_dodatnie

# 5. Ekstrakcja parametrów
intensywnosc = np.abs(E_f_norm)**2
faza = np.unwrap(np.angle(E_f_norm))
czesc_rzeczywista = np.real(E_f_norm)
czesc_urojona = np.imag(E_f_norm)

# 6. Konwersja osi poziomej
dlugosc_fali_nm = (299792458.0 / czestotliwosci) * 1e9
dane = np.column_stack((dlugosc_fali_nm, intensywnosc, faza, czesc_rzeczywista, czesc_urojona))

# 7. Sortowanie i zapis
sort_idx = np.argsort(dane[:, 0])
dane_posortowane = dane[sort_idx]
np.savetxt('Base_0mm_Speck.dat', dane_posortowane, fmt='%1.6e', delimiter='\t')