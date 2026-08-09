import subprocess
import os

class SimulatorEngineBase:
    def __init__(self, **kwargs):
        self.binary_path = kwargs.get("binary_path")
        self.config_path = kwargs.get("config_path")

    def run(self):
        raise NotImplementedError("Metoda run musi zostać zaimplementowana.")

class LWESimulator(SimulatorEngineBase):
    def run(self):
        # Pobranie katalogu, w którym znajduje się binarka (i CrystalDatabase.txt)
        binary_dir = os.path.dirname(self.binary_path)
        
        try:
            subprocess.run(
                [self.binary_path, self.config_path],
                cwd=binary_dir,  # WYMUSZENIE KATALOGU ROBOCZEGO
                check=True,
                stdout=subprocess.DEVNULL,  # Teraz można bezpiecznie z powrotem wyciszyć
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError as e:
            print(f"[Błąd Symulatora LWE] Proces zakończył się błędem. Kod wyjścia: {e.returncode}")