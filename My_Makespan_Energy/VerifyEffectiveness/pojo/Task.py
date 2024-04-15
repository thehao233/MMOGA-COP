class Task:

    # def __deepcopy__(self, memodict={}):
    #     info = Task()
    #     info.id = self.id
    #
    #     # info.islocal = self.islocal
    #     info.preTaskSet = self.preTaskSet
    #     info.sucTaskSet = self.sucTaskSet
    #     # info.exePosition = self.exePosition
    #     # info.actualFre = self.actualFre
    #
    #     info.c_i_j_k = self.c_i_j_k
    #     info.d_i_j_k = self.d_i_j_k
    #     info.o_i_j_k = self.o_i_j_k
    #
    #     # info.RT_i_l = self.RT_i_l
    #     # info.RT_i_ws = self.RT_i_ws
    #     # info.RT_i_c = self.RT_i_c
    #     # info.RT_i_wr = self.RT_i_wr
    #
    #     return info

    def __init__(self):
        self.id = None

        self.islocal = None    # Denote the task is executed locally or on cloud.
        self.preTaskSet = []   #The set of predecessor task (element is Task class).
        self.sucTaskSet = []   #The set of successor task (element is Task class).
        self.exePosition = None  # it denotes execution position (i.e., [1,2,3,4])of the task.
        self.actualFre = 1    # The actual frequency scaling factors. 实际频率
        # task的三元组
        self.c_i_j_k = None    # The number of CPU cycles required to perform task
        self.d_i_j_k = None    # The data size of the task.
        self.o_i_j_k = None    # The output data size of the task.

        # 准备完全可以开始的时间
        self.RT_i_l = None     # The ready time of task vi on a local core.
        self.RT_i_ws = None    # The ready time of task vi on the wireless sending channel.
        self.RT_i_c = None     # The ready time of task vi on the [10,20] server.
        self.RT_i_wr = None    # The ready time for the cloud to transmit back the results of task vi

        # 开始进行操作的时间
        self.ST_i_l = None     # The start time of task vi on a local core.
        self.ST_i_ws = None    # The start time of task vi on the wireless sending channel.
        self.ST_i_c = None     # The start time of task vi on the [10,20] server.
        self.ST_i_wr = None    # The start time for the cloud to transmit back the results of task vi

        # 完成操作的时间
        # 在本地的执行完成时间
        self.FT_i_l = None     # The finish time of task vj on a local core.
        # 发送到MEC的发送完成时间
        self.FT_i_ws = None    # The finish time of task vj on the wireless sending channel.
        # 在MEC执行的执行完成时间
        self.FT_i_c = None     # The finish time of task vj on the [10,20] server.
        # SMD接收完数据的接收完成时间
        self.FT_i_wr = None    # The finish time of task vj on the wireless receiving channel.
        self.energy = 0