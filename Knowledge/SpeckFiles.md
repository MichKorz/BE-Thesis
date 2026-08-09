# Dokumentacja formatu plików profilu impulsu i plików wyjściowych symulacji (LightwaveExplorer / FROG)

Pliki tekstowe definiujące profile impulsów optycznych w dziedzinie częstotliwości (zwyczajowo nazywane `Speck.dat`, `pulse1.dat` itp.), które są odczytywane przez symulator LightwaveExplorer (LWE) oraz oprogramowanie FROG, posiadają ścisłą, 5-kolumnową strukturę numeryczną. Reprezentuje ona wynik dyskretnej, szybkiej transformaty Fouriera (FFT) zespolonego pola elektrycznego.

Należy podkreślić, że **LWE nie generuje bezpośrednio gotowych plików w tym formacie po zakończeniu symulacji**. Proces odzyskiwania pełnego widma z fazą wymaga przetworzenia surowych plików binarnych w dziedzinie czasu, które program "wypluwa" jako swój standardowy wynik.

---

## 1. Tabela Struktury Danych wejściowych (Speck.dat / pulse.dat)

| Kolumna (Indeks) | Parametr | Wartości w Domenie FFT | Jednostka domyślna |
| --- | --- | --- | --- |
| **0** | Długość fali ($\lambda$) | Ujemne i Dodatnie | nm |
| **1** | Intensywność spektralna | Wartości dodatnie (oraz zera) znormalizowane | j.u. (max 1.0) |
| **2** | Faza spektralna | Ujemne i Dodatnie | rad |
| **3** | Rzeczywista część pola ($Re$) | Ujemne i Dodatnie | znormalizowane |
| **4** | Urojona część pola ($Im$) | Ujemne i Dodatnie | znormalizowane |

---

## 2. Pliki wyjściowe generowane przez symulację LWE

Po zakończeniu propagacji, LightwaveExplorer generuje standardowo pakiet trzech głównych plików (lub archiwum `.zip` je zawierające). Nazwa każdego pliku oparta jest na nazwie nadanej symulacji (np. `NazwaSymulacji`).

### A. Plik konfiguracyjny: `[NazwaSymulacji].txt`

* **Format:** Zwykły tekst (ASCII/UTF-8).
* **Zawartość:** Zapis pełnej konfiguracji numerycznej i fizycznej użytej w symulacji (wartości GDD, krok czasowy, wymiary siatki `Ntime` i `Nspace`, skany pętli *batch*, parametry kryształów).
* **Jednostki:** Standardowe jednostki fizyczne SI lub przedrostki zdefiniowane w silniku LWE (np. GDD w $10^{-30} \text{ s}^2$, czas w sekundach, dystans w metrach).
* **Zastosowanie:** Parametry z tego pliku są niezbędne do nadania fizycznego sensu surowym danym z plików binarnych. Służą do rekonstrukcji osi czasu (`timeVector`), osi przestrzennej (`spaceVector`) oraz osi częstotliwości (`frequencyVectorSpectrum`). Bez odczytania kroku czasowego $dt$ (zapisywanego jako `Time step`), pliki `.dat` są bezużytecznymi ciągami liczb.

### B. Pełne pole elektryczne: `[NazwaSymulacji]_Ext.dat`

* **Format:** Binarny, ciągły strumień liczb zmiennoprzecinkowych podwójnej precyzji (`double`, 64-bit). Dane są ułożone w strukturze kolumnowej (konwencja Fortran-order, `order='F'` w NumPy), bez żadnych nagłówków tekstu.
* **Zawartość:** Rzeczywiste wartości amplitudy pola elektrycznego w dziedzinie czasu $E(t)$ wyciągnięte na samym końcu zdefiniowanej propagacji przestrzennej.
* **Wymiary:** Macierz wielowymiarowa. Rozmiar i kształt zależą od parametrów wpisanych w `.txt`. Tablica zawiera osobne bloki dla polaryzacji X (`Ext_x`) oraz polaryzacji Y (`Ext_y`). Podstawowy układ to siatka czasowo-przestrzenna `(Ntime, Nspace)`. Przy aktywnych pętlach skanowania, dodawane są kolejne wymiary (np. `Nsims`).
* **Jednostki:** Pole elektryczne w woltach na metr (V/m) w zależności od zastosowanego skalowania, ułożone na dyskretnej siatce punktów czasowych i przestrzennych.
* **Zastosowanie:** Podstawa do pełnej analizy. Poprzez operację 1D FFT (lub 2D FFT) na tych danych odzyskuje się pełne zespolone widmo, fazę spektralną, dyspersję i tworzy pliki dla algorytmów FROG (jak `Speck.dat`).

### C. Zintegrowane widmo: `[NazwaSymulacji]_spectrum.dat`

* **Format:** Binarny, ciągły strumień liczb zmiennoprzecinkowych podwójnej precyzji (`double`, 64-bit), układ Fortran-order. Brak nagłówków.
* **Zawartość:** Natężenie (intensywność) widmowa scałkowana po profilu przestrzennym. LWE wylicza to natężenie automatycznie z modułu do kwadratu wewnątrz kodu C++. Tablica podzielona jest na całkowite natężenie (`spectrumTotal`) oraz rzuty na osie polaryzacji (`spectrum_x`, `spectrum_y`).
* **Wymiary:** Zredukowane względem pliku `_Ext.dat`. Oś przestrzenna jest zintegrowana (znika), a długość wektora odpowiada liczbie dodatnich częstotliwości FFT (`Nfreq = Ntime / 2 + 1`). Układ obejmuje `(Nfreq, 3, Nsims)`.
* **Zastosowanie i Ograniczenia:** Służy do błyskawicznego wykreślenia rozkładu energii i wygenerowanych harmonicznych bez konieczności ręcznego aplikowania transformaty Fouriera. **Plik ten nie posiada żadnych informacji o fazie spektralnej** (zostały utracone w operacji potęgowania modułu zespolonego) ani o czasie nadejścia impulsu. Nie nadaje się do odtwarzania parametrów fazowych.

---

## 3. Oś dziedziny w plikach wejściowych i siatka FFT

Największą pułapką w interpretacji wygenerowanego pliku `Speck.dat` jest zawartość pierwszej kolumny. Choć siatka FFT operuje w dziedzinie częstotliwości, plik zapisuje ją w **nanometrach**.

* **Asymptota i przeskok:** Siatka FFT obejmuje częstotliwości dodatnie, ujemne oraz przejście przez zerową częstotliwość ($f \rightarrow 0$). Przez relację odwrotności $\lambda = c / f$, oś dziedziny ulega rozerwaniu. W plikach objawia się to gwałtownym skokiem wartości (np. z rzędu setek tysięcy na minusie bezpośrednio do miliona na plusie).
* **Wymóg konwersji:** Do wszelkich operacji na fazie wymagana jest konwersja osi długości fal do pulsacji ($\omega$) za pomocą stałej prędkości światła:

$$\omega = \frac{2\pi c}{\lambda}$$

gdzie $c = 299.792458$ nm/fs, co zwraca pulsację w jednostce rad/fs. W celu uniknięcia błędów parsowania, zapisywane wiersze powinny być zawsze posortowane rosnąco według wartości z pierwszej kolumny.

---

## 4. Zespolone Pole Elektryczne i Intensywność

* Fala elektromagnetyczna całkowicie zdefiniowana jest w dwóch ostatnich kolumnach: rzeczywistej ($Re$) oraz urojonej ($Im$).
* **Normalizacja:** Intensywność spektralna (Kolumna 1) nie jest niezależnym parametrem, lecz kwadratem modułu wektora. Relacja ta musi być bezwzględnie zachowana:

$$I = Re^2 + Im^2$$

* Złamanie tej tożsamości poprzez sztuczną normalizację samej intensywności (z pominięciem $Re$ i $Im$) prowadzi do odrzucenia pliku przez program. Wartość $I$ powinna osiągać w szczycie amplitudę 1.0, co wymaga znormalizowania początkowego pola zespolonego przed jego podziałem na część rzeczywistą i urojoną.
* Wprowadzanie bezstratnych zmian fazowych wymaga zachowania poprawnej relacji trygonometrycznej. Sama zmiana kolumny 2 (Fazy) nie przyniesie efektu w solwerze, jeśli składowe pola nie zostaną obrócone na płaszczyźnie zespolonej o kąt $\Delta\phi$:

$$Re_{nowe} = Re \cos(\Delta\phi) - Im \sin(\Delta\phi)$$

$$Im_{nowe} = Re \sin(\Delta\phi) + Im \cos(\Delta\phi)$$

---

## 5. Odtwarzanie pliku Speck.dat z pliku _Ext.dat

Z powodu braku fazy w pliku `_spectrum.dat`, odtworzenie struktury wejściowej z wyników LWE odbywa się przez obróbkę pliku `_Ext.dat`.

Kluczowe kroki matematyczne podczas odtwarzania:

1. **Centrowanie impulsu (Redukcja opóźnienia grupowego):** Pole elektryczne w pliku binarnym rzadko znajduje się w punkcie $t=0$. Bezpośrednie nałożenie FFT na przesunięty impuls skutkuje nałożeniem na widmo bardzo stromej, oscylującej fazy liniowej $e^{-i\omega t_0}$. Faza ta ulega "owijaniu", co przy ponownym wczytywaniu do LWE niszczy dokładność interpolacji. Należy użyć kołowego przesunięcia tablicy numerycznej, aby pik amplitudy znalazł się pod indeksem `0`.
2. **Transformata Fouriera:** Wykonanie algorytmu 1D FFT na wycentrowanym w czasie polu z pominięciem częstotliwości stałej ($0\text{ Hz}$).
3. **Zespolona normalizacja:** Znalezienie absolutnego maksimum na wyjściowej tablicy zespolonej i podzielenie przez nie całej tablicy, gwarantujące spójność trygonometryczną.
4. **Ekstrakcja składowych:** Obliczenie intensywności, uciąglenie fazy (`unwrap`), rozdzielenie znormalizowanego pola na składowe $Re$ i $Im$ oraz zachowanie posortowanego ułożenia dla długości fali w nm.

---

## 6. Symetria Hermitowska i Konwencja Fazy LWE

Aby wynikowa fala świetlna po transformacie na domenę czasu była ścisłą wartością rzeczywistą, pliki wejściowe muszą zachować **symetrię hermitowską**. Dla fazy oznacza to, że jej znak po ujemnej stronie osi częstotliwości musi być odwrotny do strony dodatniej: $\Delta\phi(-\omega) = -\Delta\phi(\omega)$.

Dodatkowo, LightwaveExplorer stosuje wewnętrzną konwencję znaku ewolucji fazy, która wymusza matematyczną inwersję wyniku przed rotacją wektorową. Wzór na całkowitą zmianę fazy zadaną przez dyspersję przyjmuje formę:

1. Środek impulsu (nośna $\omega_0$) identyfikuje się wyłącznie w **dodatniej** części widma.
2. Wyznacza się odstrojenie poprzez bezwzględną odległość od nośnej:

$$\Delta\omega = \vert{}\omega\vert{} - \omega_0$$

3. Przesunięcie z rozwinięcia Taylora jest obliczane, a następnie mnożone przez inwersję konwencji i znak częstotliwości:

$$\Delta\phi = -\text{sgn}(\omega) \sum_{n=2}^{N} \frac{c_n}{n!} (\vert{}\omega\vert{} - \omega_0)^n$$

Gdzie współczynniki $c_n$ oznaczają odpowiednio GDD ($n=2$), TOD ($n=3$), FOD ($n=4$) itd.

---

## 7. Artefakty i zjawiska numeryczne w solwerze LWE

* **Brak okienkowania (Windowing):** LWE narzuca zdefiniowaną zmianę fazy na całą siatkę, ignorując rzeczywisty poziom sygnału. Wielomian dyspersji rozwija się w nieskończoność dla całego wektora.
* **Eksplozja Fazy:** Siatka FFT obejmuje mikroskopijne częstotliwości $f \rightarrow 0$ (ogromne długości fal z obrzeży widma, setki tysięcy nm). Ujęcie ich we wzorze powoduje powstawanie przesunięć fazowych rzędu miliardów radianów, nakładanych na czysty szum numeryczny zera precyzyjnego maszyny. Oprogramowanie musi radzić sobie z utratą dokładności na obrzeżach.
* **Aliasing w oknie czasowym:** Jeśli nałożona zewnętrznie dyspersja (np. potężne GDD wpisane ręcznie w `Speck.dat`) opóźni grupowe częstotliwości impulsu tak silnie, że wyjdą one poza krawędzie zdefiniowane parametrem `Time span` w GUI, impuls "zawinie się" w czasie, uderzając z drugiej strony siatki. Spowoduje to niszczycielską numeryczną interferencję, widoczną w `_Ext.dat` jako szum wysokiej częstotliwości. Zewnętrzne filtrowanie i limitowanie dodawanej fazy zapobiega temu zjawisku.
