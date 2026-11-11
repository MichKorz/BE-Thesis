# Dokumentacja formatu plików profilu impulsu i plików wyjściowych symulacji (LightwaveExplorer / FROG)

Pliki tekstowe definiujące profile impulsów optycznych w dziedzinie częstotliwości (zwyczajowo nazywane `Speck.dat`, `spec.dat`, `pulse1.dat` itp.)[cite: 2, 3], które są odczytywane przez symulator LightwaveExplorer (LWE) oraz oprogramowanie FROG, posiadają ścisłą, 5-kolumnową strukturę numeryczną[cite: 1, 2, 3]. Reprezentuje ona wynik dyskretnej, szybkiej transformaty Fouriera (FFT) zespolonej obwiedni pola elektrycznego[cite: 2, 3].

Należy podkreślić, że **LWE nie generuje bezpośrednio gotowych plików w tym formacie po zakończeniu symulacji**[cite: 3]. Proces odzyskiwania pełnego widma z fazą wymaga przetworzenia surowych plików binarnych w dziedzinie czasu, które program generuje jako swój standardowy wynik[cite: 3].

---

## 1. Tabela Struktury Danych wejściowych (Speck.dat / pulse.dat)

| Kolumna (Indeks) | Parametr | Wartości w Domenie FFT | Jednostka domyślna |
| --- | --- | --- | --- |
| **0** | Długość fali ($\lambda$)[cite: 2, 3] | Ujemne i Dodatnie (stały krok częstotliwości $\Delta\nu$)[cite: 1, 2, 3] | nm[cite: 2, 3] |
| **1** | Intensywność spektralna $I(\lambda)$[cite: 2, 3] | Dodatnie (przeskalowane przez $1/\lambda^2$ i znormalizowane)[cite: 1, 2, 3] | j.u. (max 1.0)[cite: 1, 2, 3] |
| **2** | Faza spektralna[cite: 2, 3] | Ujemne i Dodatnie (zerowana poniżej progu *Phase Blanking*)[cite: 1, 2, 3] | rad[cite: 2, 3] |
| **3** | Rzeczywista część pola ($\text{Re}$)[cite: 2, 3] | Ujemne i Dodatnie (w dziedzinie częstotliwości $\tilde{E}(\nu)$)[cite: 1, 2, 3] | j.u.[cite: 1] |
| **4** | Urojona część pola ($\text{Im}$)[cite: 2, 3] | Ujemne i Dodatnie (w dziedzinie częstotliwości $\tilde{E}(\nu)$)[cite: 1, 2, 3] | j.u.[cite: 1] |

---

## 2. Pliki wyjściowe generowane przez symulację LWE

Po zakończeniu propagacji, LightwaveExplorer generuje standardowo pakiet trzech głównych plików (lub archiwum `.zip` je zawierające)[cite: 3]. Nazwa każdego pliku oparta jest na nazwie nadanej symulacji (np. `NazwaSymulacji`)[cite: 3].

### A. Plik konfiguracyjny: `[NazwaSymulacji].txt`

* **Format:** Zwykły tekst (ASCII/UTF-8)[cite: 3].
* **Zawartość:** Zapis pełnej konfiguracji numerycznej i fizycznej użytej w symulacji (wartości GDD, krok czasowy, wymiary siatki `Ntime` i `Nspace`, skany pętli *batch*, parametry kryształów)[cite: 3].
* **Jednostki:** Standardowe jednostki fizyczne SI lub przedrostki zdefiniowane w silniku LWE (np. GDD w $\text{fs}^2$ / $10^{-30} \text{ s}^2$, czas w sekundach, dystans w metrach)[cite: 3].
* **Zastosowanie:** Parametry z tego pliku są niezbędne do nadania fizycznego sensu surowym danym z plików binarnych[cite: 3]. Służą do rekonstrukcji osi czasu (`timeVector`), osi przestrzennej (`spaceVector`) oraz osi częstotliwości (`frequencyVectorSpectrum`)[cite: 3]. Bez odczytania kroku czasowego $dt$ (zapisywanego jako `Time step`), pliki `.dat` są bezużytecznymi ciągami liczb[cite: 3].

### B. Pełne pole elektryczne: `[NazwaSymulacji]_Ext.dat`

* **Format:** Binarny, ciągły strumień liczb zmiennoprzecinkowych podwójnej precyzji (`double`, 64-bit)[cite: 3]. Dane są ułożone w strukturze kolumnowej (konwencja Fortran-order, `order='F'` w NumPy), bez żadnych nagłówków tekstu[cite: 3].
* **Zawartość:** Rzeczywiste wartości amplitudy pola elektrycznego w dziedzinie czasu $E(t)$ wyciągnięte na samym końcu zdefiniowanej propagacji przestrzennej[cite: 3].
* **Wymiary:** Macierz wielowymiarowa[cite: 3]. Rozmiar i kształt zależą od parametrów wpisanych w `.txt`[cite: 3]. Tablica zawiera osobne bloki dla polaryzacji X (`Ext_x`) oraz polaryzacji Y (`Ext_y`)[cite: 3]. Podstawowy układ to siatka czasowo-przestrzenna `(Ntime, Nspace)`[cite: 3]. Przy aktywnych pętlach skanowania, dodawane są kolejne wymiary (np. `Nsims`)[cite: 3].
* **Jednostki:** Pole elektryczne w woltach na metr (V/m) w zależności od zastosowanego skalowania, ułożone na dyskretnej siatce punktów czasowych i przestrzennych[cite: 3].
* **Zastosowanie:** Podstawa do pełnej analizy[cite: 3]. Poprzez operację 1D FFT (lub 2D FFT) na tych danych odzyskuje się pełne zespolone widmo, fazę spektralną, dyspersję i tworzy pliki dla algorytmów FROG (jak `Speck.dat`)[cite: 3].

### C. Zintegrowane widmo: `[NazwaSymulacji]_spectrum.dat`

* **Format:** Binarny, ciągły strumień liczb zmiennoprzecinkowych podwójnej precyzji (`double`, 64-bit), układ Fortran-order[cite: 3]. Brak nagłówków[cite: 3].
* **Zawartość:** Natężenie (intensywność) widmowa scałkowana po profilu przestrzennym[cite: 3]. LWE wylicza to natężenie automatycznie z modułu do kwadratu wewnątrz kodu C++[cite: 3]. Tablica podzielona jest na całkowite natężenie (`spectrumTotal`) oraz rzuty na osie polaryzacji (`spectrum_x`, `spectrum_y`)[cite: 3].
* **Wymiary:** Zredukowane względem pliku `_Ext.dat`[cite: 3]. Oś przestrzenna jest zintegrowana (znika), a długość wektora odpowiada liczbie dodatnich częstotliwości FFT (`Nfreq = Ntime / 2 + 1`)[cite: 3]. Układ obejmuje `(Nfreq, 3, Nsims)`[cite: 3].
* **Zastosowanie i Ograniczenia:** Służy do błyskawicznego wykreślenia rozkładu energii i wygenerowanych harmonicznych bez konieczności ręcznego aplikowania transformaty Fouriera[cite: 3]. **Plik ten nie posiada żadnych informacji o fazie spektralnej** (zostały utracone w operacji potęgowania modułu zespolonego) ani o czasie nadejścia impulsu[cite: 3]. Nie nadaje się do odtwarzania parametrów fazowych[cite: 3].

---

## 3. Oś dziedziny w plikach wejściowych i siatka FFT

Największą pułapką w interpretacji wygenerowanego przez FROG pliku `Speck.dat` jest zawartość pierwszej kolumny[cite: 3]. Choć siatka FFT operuje w dziedzinie **stałego kroku częstotliwości ($\Delta\nu = \text{const}$)**[cite: 2], plik zapisuje odpowiadające tym częstotliwościom wartości długości fali w **nanometrach** ($\lambda_k = c / \nu_k$) без wykonywania interpolacji (re-binningu) osi[cite: 1, 2, 3].

* **Pochodzenie wartości ujemnych i asymptota:** Oprogramowanie FROG wykonuje transformatę FFT na zespolonej obwiedni pola wokół częstotliwości nośnej $\nu_0 = c/\lambda_0$[cite: 2]. Jeśli krok czasowy $\Delta t$ jest bardzo mały, szerokość okna częstotliwościowego FFT ($1/\Delta t$) staje się tak duża, że dla skrajnych ujemnych odstrojeń $\Delta\nu_k$ całkowita częstotliwość optyczna spada poniżej zera ($\nu_k = \nu_0 + \Delta\nu_k < 0$). Ponieważ program oblicza długość fali jako $\lambda_k = c / \nu_k$, przejście siatki przez $\nu = 0$ powoduje rozerwanie osi i gwałtowny skok wartości z dużych liczb ujemnych (np. `-254557.4 nm`) bezpośrednio do dużych liczb dodatnich (np. `+1058561 nm`)[cite: 1, 3]. Wiersze o $\lambda < 0$ są niefizycznym artefaktem zbyt szerokiego okna FFT i są automatycznie odrzucane przez solver UPPE w LWE.
* **Porządek sortowania wierszy:** W oryginalnym pliku `Speck.dat` wiersze są uporządkowane rosnąco według **częstotliwości $\nu_k$** (kolejnych binów tablicy FFT)[cite: 1, 2]. Oznacza to, że wartości długości fali $\lambda$ w pierwszej kolumnie **maleją** w ramach gałęzi ujemnej, skaczą przez nieskończoność i ponownie **maleją** w ramach gałęzi dodatniej (np. od `+1058561 nm` do `+337.41 nm`)[cite: 1]. Podczas generowania własnego pliku dla LWE najbezpieczniej jest zapisać wyłącznie fizyczną, dodatnią gałąź częstotliwości ($\lambda > 0$) uporządkowaną monotonicznie (według rosnącej częstotliwości lub rosnącej długości fali).
* **Wymóg konwersji:** Do wszelkich operacji na fazie (np. dodawania dyspersji czy symulacji modulatora SLM) wymagana jest konwersja dodatnich długości fal do częstości kołowej ($\omega$) za pomocą prędkości światła[cite: 3]:

$$\omega = \frac{2\pi c}{\lambda}$$

gdzie $c = 299.792458\text{ nm/fs}$, co zwraca częstość kołową w jednostce $\text{rad/fs}$[cite: 3].

---

## 4. Zespolone Pole Elektryczne, Intensywność i Faza

* **Konstrukcja pola z kolumn 3 i 4:** Zarówno program FROG, jak i LWE przy wczytywaniu pliku konstruują zespolone pole elektryczne **wyłącznie na podstawie dwóch ostatnich kolumn**: części rzeczywistej ($\text{Re}$, kolumna 3) oraz urojonej ($\text{Im}$, kolumna 4)[cite: 2, 3]. Kolumny 1 (`Intensity`) oraz 2 (`Phase`) są ignorowane podczas wczytywania pola do solwera i pełnią w pliku funkcję pomocniczą (do bezpośredniego podglądu wykresów)[cite: 2].
* **Skalowanie Jakobianem w kolumnie Intensywności:** W pliku `Speck.dat` kolumny $\text{Re}$ oraz $\text{Im}$ reprezentują pole zespolone w dziedzinie częstotliwości $\tilde{E}(\nu_k)$[cite: 2]. Natomiast kolumna 1 (`Intensity`) reprezentuje gęstość widmową mocy przeliczoną na stały krok długości fali $I(\lambda_k)$ (jak w spektrometrze siatkowym), znormalizowaną do `1.0` w maksimum[cite: 1, 2]. Z tego powodu między kolumnami **nie zachodzi** prosta równość $I = \text{Re}^2 + \text{Im}^2$, lecz relacja uwzględniająca Jakobian transformacji ($\vert{}d\nu/d\lambda\vert{} = c/\lambda^2$)[cite: 1, 2]:

$$I(\lambda_k) = \text{const} \cdot \frac{\text{Re}_k^2 + \text{Im}_k^2}{\lambda_k^2}$$

* **Progowanie fazy (*Phase Blanking*):** W oryginalnym pliku `Speck.dat` kolumna 2 (`Phase`) zawiera zera we wszystkich wierszach, w których znormalizowana intensywność spada poniżej zadanego progu (domyślnie $10^{-4}$), mimo że kolumny $\text{Re}$ i $\text{Im}$ zachowują tam pełne, niezerowe wartości szumu[cite: 1, 2]. **Nie wolno zerować fazy w kolumnach $\text{Re}$ i $\text{Im}$ dla szumu w skrzydłach widma** (np. przez zastąpienie $\text{Re} = \sqrt{\text{Re}^2 + \text{Im}^2}, \text{Im} = 0$), ponieważ wymusi to koherencję szumu z całego pasma FFT i wygeneruje niefizyczny pik pasożytniczy w $t = 0$. Jeśli konieczne jest oczyszczenie szumu w skrzydłach, należy wyzerować lub stłumić gładkim oknem (np. super-Gaussa) całą amplitudę zespoloną ($\text{Re} = 0$ oraz $\text{Im} = 0$).
* **Modyfikacja fazy (np. symulacja SLM / dyspersji):** Wprowadzanie zmian fazowych $\Delta\phi(\omega)$ wymaga obrócenia wektora $(\text{Re}, \text{Im})$ na płaszczyźnie zespolonej o kąt $\Delta\phi$[cite: 3]:

$$\text{Re}_{\text{nowe}} = \text{Re} \cos(\Delta\phi) - \text{Im} \sin(\Delta\phi)$$

$$\text{Im}_{\text{nowe}} = \text{Re} \sin(\Delta\phi) + \text{Im} \cos(\Delta\phi)$$

---

## 5. Odtwarzanie pliku Speck.dat z pliku _Ext.dat

Z powodu braku fazy w pliku `_spectrum.dat`, odtworzenie struktury wejściowej z wyników LWE odbywa się przez obróbkę pliku `_Ext.dat`[cite: 3].

Kluczowe kroki matematyczne podczas odtwarzania[cite: 3]:

1. **Centrowanie impulsu (Redukcja opóźnienia grupowego):** Pole elektryczne w pliku binarnym rzadko znajduje się dokładnie w punkcie $t=0$[cite: 3]. Bezpośrednie nałożenie FFT na przesunięty w czasie impuls skutkuje nałożeniem na widmo bardzo stromej, szybkozmiennej fazy liniowej $e^{-i\omega t_0}$[cite: 3]. Faza ta powoduje silne oscylacje składowych $\text{Re}$ oraz $\text{Im}$ pomiędzy sąsiednimi punktami siatki, co przy ponownym wczytywaniu i interpolacji na inną siatkę w LWE prowadzi do sztucznego tłumienia amplitudy[cite: 3]. Należy użyć kołowego przesunięcia tablicy numerycznej (`fftshift` / `roll`), aby środek ciężkości impulsu znalazł się pod indeksem odpowiadającym $t=0$[cite: 3].
2. **Transformata Fouriera:** Wykonanie algorytmu 1D FFT (`rfft`) na wycentrowanym w czasie rzeczywistym polu elektrycznym i odrzucenie składowej stałej ($\omega = 0$)[cite: 3].
3. **Zespolona normalizacja:** Podzielenie tablicy zespolonej przez jej maksymalny moduł, co gwarantuje czytelny zakres wartości $\text{Re}$ i $\text{Im}$ (analogicznie do konwencji FROG)[cite: 1, 3].
4. **Ekstrakcja składowych:** Wyliczenie długości fali $\lambda = 2\pi c / \omega$ dla dodatnich częstości, obliczenie intensywności widmowej, uciąglenie fazy (`unwrap`) oraz zapisanie składowych $\text{Re}$ i $\text{Im}$[cite: 3].

---

## 6. Dodatnie częstości w UPPE i Konwencja Fazy

Solver UPPE w programie LightwaveExplorer operuje na sygnale analitycznym, czyli wyłącznie na **dodatnich częstotliwościach optycznych ($\omega > 0$, $\lambda > 0$)**. Przy odtwarzaniu rzeczywistego pola $E(t)$ w dziedzinie czasu symetria hermitowska ($E(-\omega) = E^*(\omega)$) jest wymuszana automatycznie wewnątrz silnika C++ poprzez transformatę odwrotną typu `c2r` (complex-to-real). Dzięki temu podczas modyfikacji pliku `Speck.dat` operujemy wyłącznie na fizycznych, dodatnich częstościach ($\omega > 0$).

Wprowadzając dyspersję za pomocą szeregu Taylora wokół częstości centralnej impulsu $\omega_0 = \frac{2\pi c}{\lambda_0}$ (gdzie $\omega > 0$):

$$\Delta\phi(\omega) = \pm \sum_{n=2}^{N} \frac{c_n}{n!} (\omega - \omega_0)^n$$

gdzie współczynniki $c_n = \left.\frac{\partial^n \phi}{\partial \omega^n}\right\vert{}_{\omega_0}$ oznaczają odpowiednio GDD ($n=2$), TOD ($n=3$), FOD ($n=4$) itd.[cite: 3] Znak ($\pm$) przed sumą zależy od przyjętej konwencji transformaty Fouriera ($\exp(-i\omega t)$ vs $\exp(+i\omega t)$) między plikiem wejściowym a silnikiem LWE i warto go potwierdzić testem pojedynczego impulsu z dodatnim GDD[cite: 3].

---

## 7. Artefakty i zjawiska numeryczne przy modyfikacji widma

* **Brak okienkowania (Windowing) przy nakładaniu fazy:** Nakładanie wielomianu dyspersji lub fazy SLM bezpośrednio na cały wektor danych zmienia fazę również tam, gdzie sygnał stanowi wyłącznie szum tła[cite: 3]. Dopóki modyfikacja polega wyłącznie na obrocie zespolonym ($\text{Re} + i\text{Im})e^{i\Delta\phi}$ bez zmiany modułu, poziom szumu się nie zwiększa, jednak dobrą praktyką jest odfiltrowanie niefizycznych skrzydeł widma gładkim oknem przed symulacją nieliniową.
* **Szybki wzrost fazy na krańcach wysokich częstotliwości ($\lambda \rightarrow 0, \omega \rightarrow \infty$):** Ponieważ wielomian Taylora skaluje się z potęgą odstrojenia częstotliwościowego $(\omega - \omega_0)^n$, dla bardzo krótkich fal (dużych częstości $\omega$ na końcu pliku `Speck.dat`) faza narasta bardzo gwałtownie (natomiast dla długich fal $\lambda \rightarrow \infty$, gdzie $\omega \rightarrow 0$, odstrojenie dąży do skończonej wartości $-\omega_0$). Gdy różnica fazy między sąsiednimi węzłami siatki przekroczy $\sim 0.5\text{ rad}$, liniowa interpolacja składowych $\text{Re}$ i $\text{Im}$ na nową siatkę w LWE spowoduje sztuczną utratę amplitudy. W takim przypadku przed nałożeniem dużej fazy należy najpierw zagęścić siatkę częstotliwościową pliku wejściowego.
* **Aliasing w oknie czasowym:** Jeśli nałożona zewnętrznie dyspersja (np. duże GDD lub faza SLM wpisana do `Speck.dat`) rozciągnie czasowo składowe częstotliwościowe impulsu tak silnie, że wyjdą one poza krawędzie okna czasowego zdefiniowanego parametrem `Time span` w GUI LWE, impuls „zawinie się” w czasie na przeciwległy koniec siatki (wynik cykliczności FFT)[cite: 3]. Spowoduje to niszczącą, niefizyczną interferencję widoczną w `_Ext.dat`[cite: 3]. Przy symulacji silnie rozciągniętych impulsów należy odpowiednio zwiększyć `Time span` oraz liczbę punktów siatki `Ntime` w LWE[cite: 3].