import numpy as np
import os
import zipfile

try:
    import lwe
except ImportError:
    import LightwaveExplorer as lwe

class FitnessEvaluatorBase:
    def __init__(self, **kwargs):
        self.output_path = kwargs.get("output_path")

    def evaluate(self):
        raise NotImplementedError("Metoda evaluate musi zostać zaimplementowana.")

class MidIREnergyEvaluator(FitnessEvaluatorBase):
    def __init__(self, output_path, f_min=0.1e14, f_max=0.58e14, **kwargs):
        super().__init__(output_path=output_path, **kwargs)
        self.f_min = f_min
        self.f_max = f_max

    def evaluate(self):
        if not os.path.exists(self.output_path):
            return 0.0, None, None, None
            
        try:
            sim_results = lwe.load(self.output_path)
            frequency = sim_results.frequencyVectorSpectrum
            spectrum = sim_results.spectrumTotal
            
            phase = None
            # Próba pobrania fazy z parsera LWE
            try:
                phase = sim_results.phaseSpectrumTotal
            except AttributeError:
                pass
                
            # Fallback: Ręczne rozpakowanie ZIP i odczyt z pliku
            if phase is None or len(phase) == 0:
                try:
                    with zipfile.ZipFile(self.output_path, 'r') as z:
                        # LWE nazywa widma SpectrumTotal.dat, Spectrum_1.dat, etc.
                        possible_names = ['SpectrumTotal.dat', 'Spectrum_1.dat', 'Spectrum_2.dat']
                        for name in possible_names:
                            if name in z.namelist():
                                with z.open(name) as f:
                                    dat = np.loadtxt(f)
                                    phase = dat[:, 2] # 3 kolumna to faza
                                break
                except Exception:
                    pass
            
            # Zabezpieczenie przed niewłaściwym rozmiarem macierzy
            if phase is None or len(phase) != len(spectrum):
                phase = np.zeros_like(spectrum)

            mask = (frequency >= self.f_min) & (frequency <= self.f_max)
            freq_range = frequency[mask]
            spec_range = spectrum[mask]

            if len(freq_range) == 0:
                return 0.0, frequency, spectrum, phase
                
            fitness = float(np.trapezoid(spec_range, freq_range))
            return max(0.0, fitness), frequency, spectrum, phase
            
        except Exception as e:
            print(f"[Błąd Ewaluatora] {e}")
            return 0.0, None, None, None