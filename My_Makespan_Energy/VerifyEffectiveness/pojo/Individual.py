class Individual:

    # def __deepcopy__(self, memodict={}):
    #     info = Individual()
    #     info.chromosome = self.chromosome # 基因位是SMD类型    保存的为SMD个体    染色体
    #     info.fitness = self.fitness
    #     info.isFeasible = self.isFeasible  # 判断该个体是否合法
    #     info.temp_fitness = self.temp_fitness  # 临时适应度，计算拥挤距离的时候，按每个目标值来对类列表进行升序排序
    #     info.distance = self.distance
    #     info.rank = self.rank  # 支配等级 用于快速非支配排序
    #     info.S_p = self.S_p# 种群中此个体支配的个体集合
    #     info.n = self.n
    #     return info

    def __init__(self):

        self.chromosome = []      #基因位是SMD类型    保存的为SMD个体    染色体
        self.fitness = []
        self.isFeasible = True    #判断该个体是否合法
        self.temp_fitness = None  #临时适应度，计算拥挤距离的时候，按每个目标值来对类列表进行升序排序
        self.distance = 0.0
        self.rank = None # 支配等级 用于快速非支配排序
        self.S_p = []  #种群中此个体支配的个体集合
        self.n = 0  #种群中支配此个体的个数

        self.exe_on_edge = 0