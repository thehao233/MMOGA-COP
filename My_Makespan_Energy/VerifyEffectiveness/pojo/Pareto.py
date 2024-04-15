class Pareto:
    # def __deepcopy__(self, memodict={}):
    #     info = Pareto()
    #     info.fitness = self.fitness
    #     info.temp_fitness = self.temp_fitness
    #     info.chromosome = self.chromosome
    #     return info

    def __init__(self):
        self.chromosome = None
        self.fitness = []
        self.temp_fitness = None  #排序使用