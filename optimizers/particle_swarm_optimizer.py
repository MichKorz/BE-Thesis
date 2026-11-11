import numpy as np

class ParticleSwarmOptimizer:
    def __init__(self, bounds, num_agents=50, max_iterations=20, **kwargs):
        self.bounds = np.array(bounds)
        self.num_agents = num_agents
        self.max_iterations = max_iterations
        self.dim = len(self.bounds)
        
        self.current_iteration = 0
        
        # Parametry hiper-optymalizacyjne PSO (standardowe wartości Clerc-Kennedy)
        self.w = kwargs.get('w', 0.729)   # Waga bezwładności (inertia weight)
        self.c1 = kwargs.get('c1', 1.494) # Współczynnik kognitywny (personal best)
        self.c2 = kwargs.get('c2', 1.494) # Współczynnik społeczny (global best)
        
        self.population = None
        self.velocities = np.zeros((self.num_agents, self.dim))
        
        # Pamięć najlepszych wyników dla każdego osobnika (Personal Best)
        self.pbest_pos = np.zeros((self.num_agents, self.dim))
        self.pbest_scores = np.full(self.num_agents, float('-inf'))
        
        # Pamięć najlepszego wyniku globalnego (Global Best)
        self.gbest_pos = np.zeros(self.dim)
        self.gbest_score = float('-inf')

    def ask(self):
        """Zwraca nową populację do ewaluacji."""
        if self.current_iteration == 0:
            # Pierwsza iteracja: losowa inicjalizacja pozycji i zerowe prędkości
            lower_bounds = self.bounds[:, 0]
            upper_bounds = self.bounds[:, 1]
            self.population = np.random.uniform(
                low=lower_bounds, 
                high=upper_bounds, 
                size=(self.num_agents, self.dim)
            )
            # Kopiowanie początkowych pozycji jako pbest
            self.pbest_pos = self.population.copy()
            return self.population.tolist()
        
        # Losowe macierze dla aktualizacji roju
        r1 = np.random.random((self.num_agents, self.dim))
        r2 = np.random.random((self.num_agents, self.dim))
        
        # Aktualizacja wektora prędkości
        self.velocities = (self.w * self.velocities + 
                           self.c1 * r1 * (self.pbest_pos - self.population) + 
                           self.c2 * r2 * (self.gbest_pos - self.population))
        
        # Aktualizacja pozycji
        new_population = self.population + self.velocities
        
        # Obcinanie do granic i zerowanie prędkości na krawędziach
        lower_bounds = self.bounds[:, 0]
        upper_bounds = self.bounds[:, 1]
        
        out_of_bounds_low = new_population < lower_bounds
        out_of_bounds_high = new_population > upper_bounds
        
        # Jeśli cząstka uderza w granicę, traci pęd (prędkość = 0)
        self.velocities[out_of_bounds_low | out_of_bounds_high] = 0.0
        
        new_population = np.clip(new_population, lower_bounds, upper_bounds)
        self.population = new_population
        
        return self.population.tolist()

    def tell(self, population, fitnesses):
        """Aktualizuje stany roju na podstawie wyników ewaluacji."""
        population = np.array(population)
        
        for i in range(self.num_agents):
            fitness = fitnesses[i]
            
            # Aktualizacja Personal Best
            if fitness > self.pbest_scores[i]:
                self.pbest_scores[i] = fitness
                self.pbest_pos[i] = population[i].copy()
                
                # Aktualizacja Global Best
                if fitness > self.gbest_score:
                    self.gbest_score = fitness
                    self.gbest_pos = population[i].copy()
                    
        self.population = population
        self.current_iteration += 1

    def load_state(self, history):
        """Odbudowuje stan optymalizatora na podstawie zapisanej historii z JSON."""
        for iteration_data in history:
            population = []
            fitnesses = []
            
            for eval_data in iteration_data["evaluations"]:
                population.append(eval_data["vector"])
                fitnesses.append(eval_data["fitness"])
                
            population_arr = np.array(population)
            
            # Rekonstrukcja wektora prędkości na podstawie historycznych skoków
            if self.population is not None:
                self.velocities = population_arr - self.population
                
            self.tell(population, fitnesses)