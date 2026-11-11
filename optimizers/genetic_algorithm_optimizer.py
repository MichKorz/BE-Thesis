import numpy as np

class GeneticAlgorithmOptimizer:
    def __init__(self, bounds, num_agents=50, max_iterations=20, **kwargs):
        self.bounds = np.array(bounds)
        self.num_agents = num_agents
        self.max_iterations = max_iterations
        self.dim = len(self.bounds)
        
        self.current_iteration = 0
        
        # Parametry specyficzne dla algorytmu genetycznego (ciągłego)
        self.mutation_rate = kwargs.get('mutation_rate', 0.1)
        self.crossover_rate = kwargs.get('crossover_rate', 0.8)
        self.elite_count = kwargs.get('elite_count', 2) # Liczba najlepszych osobników przechodzących bez zmian
        self.tournament_size = kwargs.get('tournament_size', 3)
        
        self.population = None
        self.fitnesses = None
        
        # Pamięć najlepszego wyniku globalnego
        self.best_pos = np.zeros(self.dim)
        self.best_score = float('-inf')

    def ask(self):
        """Zwraca nową populację do ewaluacji."""
        if self.current_iteration == 0:
            # Pierwsza iteracja: losowa inicjalizacja w granicach
            lower_bounds = self.bounds[:, 0]
            upper_bounds = self.bounds[:, 1]
            self.population = np.random.uniform(
                low=lower_bounds, 
                high=upper_bounds, 
                size=(self.num_agents, self.dim)
            )
            return self.population.tolist()
        
        new_population = []
        
        # 1. Elitaryzm: Kopiowanie najlepszych osobników z poprzedniej generacji
        sorted_indices = np.argsort(self.fitnesses)[::-1] # Sortowanie malejąco
        for i in range(self.elite_count):
            new_population.append(self.population[sorted_indices[i]].copy())
            
        # 2. Generowanie reszty populacji (krzyżowanie i mutacja)
        while len(new_population) < self.num_agents:
            # Selekcja
            parent1 = self._tournament_selection()
            parent2 = self._tournament_selection()
            
            # Krzyżowanie (Arytmetyczne dla zmiennych ciągłych)
            if np.random.random() < self.crossover_rate:
                alpha = np.random.random(self.dim)
                child1 = alpha * parent1 + (1.0 - alpha) * parent2
                child2 = (1.0 - alpha) * parent1 + alpha * parent2
            else:
                child1, child2 = parent1.copy(), parent2.copy()
                
            # Mutacja
            child1 = self._mutate(child1)
            child2 = self._mutate(child2)
            
            new_population.append(child1)
            if len(new_population) < self.num_agents:
                new_population.append(child2)
                
        new_population = np.array(new_population)
        
        # Obcinanie do dopuszczalnych granic poszukiwań
        new_population = np.clip(new_population, self.bounds[:, 0], self.bounds[:, 1])
        self.population = new_population
        
        return self.population.tolist()

    def tell(self, population, fitnesses):
        """Aktualizuje stany na podstawie wyników ewaluacji."""
        self.population = np.array(population)
        self.fitnesses = np.array(fitnesses)
        
        # Aktualizacja globalnie najlepszego wyniku (Elitaryzm)
        best_idx = np.argmax(self.fitnesses)
        if self.fitnesses[best_idx] > self.best_score:
            self.best_score = self.fitnesses[best_idx]
            self.best_pos = self.population[best_idx].copy()
            
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

    def _tournament_selection(self):
        """Selekcja turniejowa wybierająca najlepszego osobnika z wylosowanej grupy."""
        indices = np.random.choice(self.num_agents, self.tournament_size, replace=False)
        tournament_fitnesses = self.fitnesses[indices]
        best_local_idx = np.argmax(tournament_fitnesses)
        best_global_idx = indices[best_local_idx]
        return self.population[best_global_idx]

    def _mutate(self, child):
        """Mutacja Gaussowska z rozrzutem bazującym na szerokości granic."""
        for i in range(self.dim):
            if np.random.random() < self.mutation_rate:
                # Standardowe odchylenie mutacji to 10% przedziału danej zmiennej
                mutation_range = (self.bounds[i, 1] - self.bounds[i, 0]) * 0.1
                child[i] += np.random.normal(0, mutation_range)
        return child