import numpy as np

class GreyWolfOptimizer:
    def __init__(self, bounds, num_agents=50, max_iterations=20, **kwargs):
        self.bounds = np.array(bounds)
        self.num_agents = num_agents
        self.max_iterations = max_iterations
        self.dim = len(self.bounds)
        
        self.current_iteration = 0
        self.population = None
        
        # Inicjalizacja hierarchii z wartością domyślną minus nieskończoność (maksymalizacja)
        self.alpha_pos = np.zeros(self.dim)
        self.alpha_score = float('-inf')
        
        self.beta_pos = np.zeros(self.dim)
        self.beta_score = float('-inf')
        
        self.delta_pos = np.zeros(self.dim)
        self.delta_score = float('-inf')

    def ask(self):
        """Zwraca nową populację do ewaluacji."""
        if self.current_iteration == 0:
            # Pierwsza iteracja: losowa inicjalizacja agentów w granicach
            lower_bounds = self.bounds[:, 0]
            upper_bounds = self.bounds[:, 1]
            self.population = np.random.uniform(
                low=lower_bounds, 
                high=upper_bounds, 
                size=(self.num_agents, self.dim)
            )
            return self.population.tolist()
        
        # Kolejne iteracje: aktualizacja pozycji na podstawie równań GWO
        # Parametr 'a' maleje liniowo od 2 do 0
        a = 2.0 - 2.0 * (self.current_iteration / self.max_iterations)
        
        new_population = np.zeros((self.num_agents, self.dim))
        
        for i in range(self.num_agents):
            for j in range(self.dim):
                # Parametry dla wilka Alpha
                r1 = np.random.random()
                r2 = np.random.random()
                A1 = 2.0 * a * r1 - a
                C1 = 2.0 * r2
                D_alpha = abs(C1 * self.alpha_pos[j] - self.population[i, j])
                X1 = self.alpha_pos[j] - A1 * D_alpha
                
                # Parametry dla wilka Beta
                r1 = np.random.random()
                r2 = np.random.random()
                A2 = 2.0 * a * r1 - a
                C2 = 2.0 * r2
                D_beta = abs(C2 * self.beta_pos[j] - self.population[i, j])
                X2 = self.beta_pos[j] - A2 * D_beta
                
                # Parametry dla wilka Delta
                r1 = np.random.random()
                r2 = np.random.random()
                A3 = 2.0 * a * r1 - a
                C3 = 2.0 * r2
                D_delta = abs(C3 * self.delta_pos[j] - self.population[i, j])
                X3 = self.delta_pos[j] - A3 * D_delta
                
                # Nowa pozycja (średnia)
                new_population[i, j] = (X1 + X2 + X3) / 3.0
                
        # Obcinanie do granic poszukiwań (zapobieganie opuszczeniu dopuszczalnej przestrzeni)
        new_population = np.clip(new_population, self.bounds[:, 0], self.bounds[:, 1])
        self.population = new_population
        
        return self.population.tolist()

    def tell(self, population, fitnesses):
        """Aktualizuje stany hierarchii na podstawie wyników ewaluacji."""
        population = np.array(population)
        
        for i in range(self.num_agents):
            fitness = fitnesses[i]
            
            # Aktualizacja Alpha, Beta i Delta (Maksymalizacja)
            if fitness > self.alpha_score:
                self.delta_score = self.beta_score
                self.delta_pos = self.beta_pos.copy()
                self.beta_score = self.alpha_score
                self.beta_pos = self.alpha_pos.copy()
                self.alpha_score = fitness
                self.alpha_pos = population[i].copy()
            
            elif fitness > self.beta_score:
                self.delta_score = self.beta_score
                self.delta_pos = self.beta_pos.copy()
                self.beta_score = fitness
                self.beta_pos = population[i].copy()
                
            elif fitness > self.delta_score:
                self.delta_score = fitness
                self.delta_pos = population[i].copy()
                
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
                
            self.tell(population, fitnesses)