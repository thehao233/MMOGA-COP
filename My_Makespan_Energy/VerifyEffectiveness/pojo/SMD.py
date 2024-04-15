import random

from My_Makespan_Energy.VerifyEffectiveness.pojo.Workflow import Workflow


class SMD:

    # def __deepcopy__(self, memodict={}):
    #     info = SMD()
    #     info.coordinate = self.coordinate
    #     info.workflow = self.workflow
    #     info.channel = self.channel
    #     info.coreLevel = self.coreLevel
    #     info.numOfCoreLevel = self.numOfCoreLevel
    #     info.r = self.r
    #     info.g_i_j = self.g_i_j
    #     info.R_i_j = self.R_i_j
    #     info.I_i_j = self.I_i_j
    #     info.coreCC = self.coreCC
    #     info.pcc_i_j = self.pcc_i_j
    #     info.pws_i_j = self.pws_i_j
    #     info.pwr_i_j = self.pwr_i_j
    #     return info

    def __init__(self):
        self.coordinate = []    # The position coordination of the SMD
        self.workflow = Workflow()      #The workflow of the SMD 一个SMD包含一个工作流
        self.channel = None     # Gaining channel index  所获得的和SeNB之间的通道

        self.coreLevel = [0.2,0.5,0.8,1]
        self.numOfCoreLevel = self.coreLevel.__len__()
        self.r = 2

        # SMD与SeNB之间的信道增益
        self.g_i_j = None       # The channel gain between the SMD and SeNB Sm
        # SMD所能实现的上行传输速率
        self.R_i_j = None       # The data transmission rate of the SMD
        # 与SMD有关的干扰系数，指示信道共享的严重程度
        self.I_i_j = None       # The interference at the SMD

        #The SMD is modeled as a 3-tuple
        # 分别包含三个核心和MEC的计算能力
        self.coreCC = {1:None, 2:None, 3:None, 4:4}        # The computing capacity of three core.
        # 分别表示三个核心在最大功率下的功耗
        self.pcc_i_j = {1:4, 2:2, 3:1}  # The power consumption of the three cores under the maximum operating frequency.
        # SMD发送数据的功率(w)
        self.pws_i_j = 0.5  # The send data power (w) of the SMD
        # SMD接收数据的功率(w)
        self.pwr_i_j = 0.1  # The receive data power (w) of the SMD

        self.taskSetTxt = []  # 记录SMD中任务集的文本