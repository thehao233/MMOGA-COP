class Schedule:

    # def __deepcopy__(self, memodict={}):
    #     info = Schedule()
    #     info.taskSet = self.taskSet
    #     return info

    def __init__(self):

        self.taskSet = {}

        # Record the set of task that is executed certain execution unit selection.
        # eg. S[3]=[v1,v3,v5,v7,v9,v10]
        self.S = {1:[], 2:[], 3:[], 4:[]}
        # Index is core number, its element denotes the current time point on the core.
        self.coreTP = {1:[0], 2:[0], 3:[0]}
        # The current time point on the wireless sending channel.
        self.wsTP = [0] # 每个任务完成数据传输(到MEC)的时间
        # The current time point on the cloud.
        self.MECTP = [0] # 每个task完成数据执行(在MEC)的时间
        # The current time point on the wireless receiving channel.
        self.wrTP = [0]  # 每个task完成数据返回(到SMD)的时间
        self.T_total = None
        self.E_total = 0
        self.TimeEnergy = []