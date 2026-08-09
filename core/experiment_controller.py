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

    def run(self, max_iterations):
        self.state["status"] = "running"
        self._save_json()
        
        visualizer = ExperimentVisualizer()
        
        # Wczytanie oryginalnego pliku wejściowego do wizualizatora
        original_pulse_data = np.loadtxt(self.modulator.input_path)
        
        best_result_data = None
        best_modified_data = None
        if self.state["best_vector"] is not None:
            wl, orig_int, final_phase = self.modulator.apply_phase(self.state["best_vector"])
            self.simulator.run()
            fit, f, s, p = self.evaluator.evaluate()
            best_result_data = {"freq": f, "spec": s, "phase": p, "fitness": fit}
            # Odczyt zmodyfikowanej macierzy po restarcie z dysku
            best_modified_data = np.loadtxt(self.modulator.output_path)
            
        while self.state["current_iteration"] < max_iterations:
            iteration_start_time = time.time()
            iteration = self.state["current_iteration"] + 1
            
            population = self.optimizer.ask()
            evaluation_results = []
            
            iter_best_fitness = float('-inf')
            iter_best_mask = None
            iter_best_result = None
            
            for eval_id, vector in enumerate(population):
                eval_start_time = time.time()
                
                wl, orig_int, final_phase = self.modulator.apply_phase(vector)
                
                self.simulator.run()
                fitness, f, s, p = self.evaluator.evaluate()
                execution_time = int((time.time() - eval_start_time) * 1000)
                
                evaluation_results.append({
                    "id": eval_id,
                    "vector": vector,
                    "fitness": fitness,
                    "execution_time_ms": execution_time
                })
                
                if fitness >= iter_best_fitness:
                    iter_best_fitness = fitness
                    # Wczytanie z dysku po aplikacji fazy, przed wysłaniem do wizualizatora
                    iter_best_mask = np.loadtxt(self.modulator.output_path)
                    iter_best_result = {"freq": f, "spec": s, "phase": p, "fitness": fitness}
                
                if fitness > self.state["best_fitness"]:
                    self.state["best_fitness"] = fitness
                    self.state["best_vector"] = vector
                    best_result_data = iter_best_result
                    best_modified_data = iter_best_mask

            fitness_values = [res["fitness"] for res in evaluation_results]
            self.optimizer.tell(population, fitness_values)
            
            self.history.append({
                "iteration": iteration,
                "iteration_time_s": time.time() - iteration_start_time,
                "evaluations": evaluation_results
            })
            
            self.state["current_iteration"] = iteration
            self._save_json()
            
            visualizer.update_plot(
                iteration=iteration,
                max_iterations=max_iterations,
                history=self.history,
                modified_pulse_data=iter_best_mask,
                original_pulse_data=original_pulse_data,
                best_result_data=best_result_data,
                eval_f_min=self.evaluator.f_min,
                eval_f_max=self.evaluator.f_max
            )
            
            print(f"Iteracja {iteration}/{max_iterations} zakończona. Najlepszy fitness: {self.state['best_fitness']:.4e}")

        self.state["status"] = "completed"
        self._save_json()