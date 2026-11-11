import json
import os
import time
import importlib
import tempfile
import shutil
import numpy as np
from core.visualizer import ExperimentVisualizer

class ExperimentController:
    def __init__(self, json_path):
        self.json_path = json_path
        self.data = self._load_json(json_path)
        
        self.config = self.data["config"]
        self.state = self.data.setdefault("state", {
            "status": "initialized",
            "current_iteration": 0,
            "best_fitness": float('-inf'),
            "best_vector": None
        })
        self.history = self.data.setdefault("history", [])
        
        self._initialize_components()

    def _load_json(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Brak pliku konfiguracyjnego: {path}")
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)

    def _save_json(self):
        dir_name = os.path.dirname(self.json_path)
        fd, temp_path = tempfile.mkstemp(dir=dir_name, prefix="tmp_chkpt_", suffix=".json")
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as f:
                json.dump(self.data, f, indent=2)
            shutil.move(temp_path, self.json_path)
        except Exception as e:
            os.remove(temp_path)
            raise e

    def _instantiate_module(self, module_path, class_name, kwargs, optional):
        module = importlib.import_module(module_path)
        cls = getattr(module, class_name)
        init_params = kwargs.copy()
        if optional:
            init_params.update(optional)
        return cls(**init_params)

    def _initialize_components(self):
        mod_conf = self.config["modulator"]
        self.modulator = self._instantiate_module(mod_conf["module"], mod_conf["class"], mod_conf.get("kwargs", {}), mod_conf.get("optional", {}))
        
        eval_conf = self.config["evaluator"]
        self.evaluator = self._instantiate_module(eval_conf["module"], eval_conf["class"], eval_conf.get("kwargs", {}), eval_conf.get("optional", {}))
        
        sim_conf = self.config["simulator"]
        self.simulator = self._instantiate_module(sim_conf["module"], sim_conf["class"], sim_conf.get("kwargs", {}), sim_conf.get("optional", {}))
        
        opt_conf = self.config["optimizer"]
        opt_kwargs = opt_conf.get("kwargs", {}).copy()
        opt_kwargs["bounds"] = self.modulator.get_bounds()
        
        self.optimizer = self._instantiate_module(opt_conf["module"], opt_conf["class"], opt_kwargs, opt_conf.get("optional", {}))
        
        if self.history:
            self.optimizer.load_state(self.history)

    def _generate_best_speck(self, zip_path, output_speck_path):
        import LightwaveExplorer as lwe
        try:
            symulacja = lwe.load(zip_path, True)
            srodek_x = symulacja.Nspace // 2

            if len(symulacja.Ext_x.shape) == 2:
                E_t = symulacja.Ext_x[:, srodek_x]
            else:
                E_t = symulacja.Ext_x[:, srodek_x].flatten()[:symulacja.Ntime]

            indeks_szczytu = np.argmax(np.abs(E_t))
            E_t_wycentrowane = np.roll(E_t, -indeks_szczytu)

            E_f = np.fft.fft(E_t_wycentrowane)
            E_f_dodatnie = E_f[1:symulacja.Nfreq]
            czestotliwosci = symulacja.frequencyVectorSpectrum[1:]

            max_amplituda = np.max(np.abs(E_f_dodatnie))
            if max_amplituda > 0:
                E_f_norm = E_f_dodatnie / max_amplituda
            else:
                E_f_norm = E_f_dodatnie

            intensywnosc = np.abs(E_f_norm)**2
            faza = np.unwrap(np.angle(E_f_norm))
            czesc_rzeczywista = np.real(E_f_norm)
            czesc_urojona = np.imag(E_f_norm)

            dlugosc_fali_nm = (299792458.0 / czestotliwosci) * 1e9
            dane = np.column_stack((dlugosc_fali_nm, intensywnosc, faza, czesc_rzeczywista, czesc_urojona))

            sort_idx = np.argsort(dane[:, 0])
            dane_posortowane = dane[sort_idx]
            np.savetxt(output_speck_path, dane_posortowane, fmt='%1.6e', delimiter='\t')
            return True
        except Exception as e:
            print(f"[Błąd generowania Speck.dat z zip] {e}")
            return False

    def run(self, max_iterations):
        self.state["status"] = "running"
        self._save_json()
        
        visualizer = ExperimentVisualizer()
        
        original_pulse_data = np.loadtxt(self.modulator.input_path)
        
        print("Obliczanie przebiegu bazowego (Baseline)...")
        zero_vector = [0.0] * len(self.modulator.get_bounds())
        self.modulator.apply_phase(zero_vector)
        self.simulator.run()
        
        base_zip_path = os.path.join(os.path.dirname(self.evaluator.output_path), "base/results.zip")
        base_speck_path = "/mnt/ramdisk/base_Speck.dat"
        if os.path.exists(self.evaluator.output_path):
            shutil.copy2(self.evaluator.output_path, base_zip_path)
            self._generate_best_speck(base_zip_path, base_speck_path)
            print(f"Zapisano wyniki bazowe do: {base_zip_path}")
            
        base_out_pulse_data = None
        if os.path.exists(base_speck_path):
            base_out_pulse_data = np.loadtxt(base_speck_path)
        
        best_result_data = None
        best_modified_data = None
        
        if self.state["best_vector"] is not None:
            self.modulator.apply_phase(self.state["best_vector"])
            self.simulator.run()
            fit, f, s, p = self.evaluator.evaluate()
            best_result_data = {"freq": f, "spec": s, "phase": p, "fitness": fit}
            best_modified_data = np.loadtxt(self.modulator.output_path)
            
        while self.state["current_iteration"] < max_iterations:
            iteration_start_time = time.time()
            iteration = self.state["current_iteration"] + 1
            
            population = self.optimizer.ask()
            pop_size = len(population)
            
            evaluations = []
            iter_best_fitness = float('-inf')
            iter_best_mask = None
            
            for eval_id, vector in enumerate(population):
                eval_start_time = time.time()
                
                self.modulator.apply_phase(vector)
                
                self.simulator.run()
                fitness, f, s, p = self.evaluator.evaluate()
                execution_time = int((time.time() - eval_start_time) * 1000)
                
                modified_data = np.loadtxt(self.modulator.output_path)
                
                evaluations.append({
                    "id": eval_id,
                    "vector": vector,
                    "fitness": fitness,
                    "execution_time_ms": execution_time
                })
                
                if fitness >= iter_best_fitness:
                    iter_best_fitness = fitness
                    iter_best_mask = modified_data
                
                if fitness > self.state["best_fitness"]:
                    self.state["best_fitness"] = fitness
                    self.state["best_vector"] = vector
                    best_result_data = {"freq": f, "spec": s, "phase": p, "fitness": fitness}
                    best_modified_data = modified_data
                    
                    best_speck_path = "/mnt/ramdisk/best_Speck.dat"
                    # Otwiera archiwa na bieżąco, aby nie nadpisać ich w kolejnej iteracji
                    if os.path.exists(self.evaluator.output_path):
                        self._generate_best_speck(self.evaluator.output_path, best_speck_path)
                
                elapsed_time = time.time() - iteration_start_time
                progress_pct = ((eval_id + 1) / pop_size) * 100
                print(f"\rGeneracja {iteration}/{max_iterations} | Postęp: {progress_pct:.1f}% | Czas bieżący: {elapsed_time:.1f}s", end="", flush=True)

            print()

            evaluations.sort(key=lambda x: x["id"])
            fitness_values = [res["fitness"] for res in evaluations]
            
            self.optimizer.tell(population, fitness_values)
            
            self.history.append({
                "iteration": iteration,
                "iteration_time_s": time.time() - iteration_start_time,
                "evaluations": evaluations
            })
            
            self.state["current_iteration"] = iteration
            self._save_json()

            best_out_pulse_data = None
            if os.path.exists("/mnt/ramdisk/best_Speck.dat"):
                best_out_pulse_data = np.loadtxt("/mnt/ramdisk/best_Speck.dat")
            
            visualizer.update_plot(
                iteration=iteration,
                max_iterations=max_iterations,
                history=self.history,
                original_pulse_data=original_pulse_data,
                best_modified_pulse_data=best_modified_data,
                base_out_pulse_data=base_out_pulse_data,
                best_out_pulse_data=best_out_pulse_data,
                best_fitness_val=self.state["best_fitness"],
                eval_f_min=self.evaluator.f_min,
                eval_f_max=self.evaluator.f_max
            )
            
            print(f"Iteracja {iteration}/{max_iterations} zakończona. Najlepszy fitness: {self.state['best_fitness']:.4e}")

        self.state["status"] = "completed"
        self._save_json()