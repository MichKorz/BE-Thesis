import subprocess
import os

# Flaga sterująca trybem debugowania (True = pełne logi, False = wyciszenie)
DEBUG_MODE = False

class SimulatorEngineBase:
    def __init__(self, **kwargs):
        self.binary_path = kwargs.get("binary_path")
        self.config_path = kwargs.get("config_path")

    def run(self):
        raise NotImplementedError("Metoda run musi zostać zaimplementowana.")

class LWESimulator(SimulatorEngineBase):
    def run(self):
        binary_dir = os.path.dirname(self.binary_path)
        
        if DEBUG_MODE:
            print("\n=== [LWE Symulator] Start diagnostyki ===")
            print(f"Katalog roboczy (CWD): {binary_dir}")
            print(f"Binarne uruchomienie: {self.binary_path}")
            print(f"Czy binarka istnieje? {os.path.exists(self.binary_path)}")
            print(f"Plik konfiguracyjny: {self.config_path}")
            print(f"Czy plik configu istnieje? {os.path.exists(self.config_path)}")
            
            try:
                result = subprocess.run(
                    [self.binary_path, self.config_path],
                    cwd=binary_dir,
                    check=False,
                    capture_output=True,
                    text=True
                )
                
                print(f"[LWE Symulator] Zakończono (Kod: {result.returncode})")
                
                if result.stdout and result.stdout.strip():
                    print("--- LWE STDOUT ---")
                    print(result.stdout.strip())
                
                if result.stderr and result.stderr.strip():
                    print("--- LWE STDERR ---")
                    print(result.stderr.strip())
                    
            except Exception as e:
                print(f"[Krytyczny Wyjątek Python] {e}")
            finally:
                print("=== [LWE Symulator] Koniec diagnostyki ===\n")
                
        else:
            try:
                subprocess.run(
                    [self.binary_path, self.config_path],
                    cwd=binary_dir,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
            except subprocess.CalledProcessError as e:
                print(f"[Błąd Symulatora LWE] Proces zakończył się błędem. Kod wyjścia: {e.returncode}")