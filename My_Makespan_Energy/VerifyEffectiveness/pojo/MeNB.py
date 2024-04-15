import math, os, xlrd

from My_Makespan_Energy.VerifyEffectiveness.pojo.SMD import SMD
from My_Makespan_Energy.VerifyEffectiveness.pojo.LeNB import LeNB
from My_Makespan_Energy.VerifyEffectiveness.pojo.Task import Task
from My_Makespan_Energy.VerifyEffectiveness.pojo.Workflow import Workflow

import numpy as np



class MeNB:

    # def __deepcopy__(self, memodict={}):
    #     info = MeNB(self.taskNumberRange)
    #     info.taskNumberRange = self.taskNumberRange
    #     info.f_MEC = self.f_MEC
    #     info.MEC_radius = self.MEC_radius
    #     info.SeNB_radius = self.SeNB_radius
    #     info.SeNBSet = self.SeNBSet
    #     info.Bandwidth = self.Bandwidth
    #     info.N = self.N
    #     info.w = self.w
    #     info.noisePower = self.noisePower
    #     info.kk = self.kk
    #     info.totalSMDNumber = self.totalSMDNumber
    #     info.codeLength = self.codeLength
    #     info.M = self.M
    #     info.H = self.H
    #     return info


    def __init__(self,taskNumberRange):

        self.taskNumberRange = taskNumberRange
        # self.taskNumberRange = '[25,35]'
        # self.taskNumberRange = '[10,40]'
        # self.taskNumberRange = '[20,30]'
        # self.taskNumberRange = '[15,25]'
        # self.taskNumberRange = '[30,40]'
        # self.taskNumberRange = '[10,20]'
        #---------------------------------------Problem Notation----------------------------------------

        self.f_MEC = 1  # Computation capacity of [10,20] server   MEC server的计算能力
        self.MEC_radius = 100 # MEC半径
        self.SeNB_radius = 50 # SeNB半径
        self.SeNBSet = [] # SeNB集合

        self.Bandwidth = 20             #Bandwidth
        self.N = 10             #Number of channel
        self.w = (self.Bandwidth / self.N) * pow(10, 6)  #The bandwidth of channel 每个通道的带宽
        self.noisePower = pow(10, -176/10)*pow(10, -3)  #The background noise power (50dBm = 100w)
        self.kk = pow(10, -11)  #It is a coefficient depending on the chip architecture ？
        self.totalSMDNumber = 0  # The total number of SMD in the network system
        self.codeLength = 0

        self.M = 0 # 所包含的SeNB结点数
        self.H = 3  # The number of the core in a SMD.  每个SMD中的核心数

        # 获取SeNBSet集合信息 包含：
        # 1.每个SeNB 2.每个SeNB中SMD信息
        # 3.每个SMD的基本信息(如：每个core的计算能力、功耗等)和包含的工作流信息
        # 4.每个工作流的基本信息(开始任务索引、EL和EO、任务集合等)和 Schedule信息
        # 为task三元组(所需CPU周期数、输入数据大小、输出数据大小)设置值
        # 5.Schedule信息 ？
        self.readMECNetwork()

        """
            计算上行传输速率 R i,j
        """
        self.M = self.SeNBSet.__len__() # self.M ：SeNB结点数
        self.calculateInterference() # 为每个SMD都设置信道干扰参数
        self.calculateDataTransmissionRate() # 计算上行传输速率 R i,j



        print("The total SMD number: ", self.totalSMDNumber)
        print("The code length: ", self.codeLength)
        # self.PF_ref = self.get_PF_ref()
        # self.IGDValue = None
        # self.IGD_list = []  # 保存100代的IGD值


    # 读取MEC工作网络
    def readMECNetwork(self):
        file_SMD_task_cpu =  open(self.getCurrentPath()+'\\'+self.taskNumberRange+'\SMD_Task_CPU_Cycles_Number.txt', 'r')
        file_SMD_task_data = open(self.getCurrentPath()+'\\'+self.taskNumberRange+'\SMD_Task_Data_Size.txt', 'r')
        file_SMD_output_task_data = open(self.getCurrentPath()+'\\'+self.taskNumberRange+'\SMD_Task_Output_Data_Size.txt', 'r')

        SeNB_count = -1
        with open(self.getCurrentPath()+'\\'+self.taskNumberRange+'\MEC_Network.txt', 'r') as MEC_Network:
            for line in MEC_Network:
                if(line == '---file end---\n'):
                    break
                elif(line == 'SeNB:\n'):
                    SeNB_count += 1
                    senb = LeNB()     #create SeNB cell
                    if(MEC_Network.readline() == 'Coordinate:\n'): # 放入SeNB的坐标
                        SeNB_crd = MEC_Network.readline()
                        SeNB_crd = SeNB_crd.splitlines()
                        SeNB_crd = SeNB_crd[0].split('  ')
                        senb.coordinate.append(float(SeNB_crd[0]))
                        senb.coordinate.append(float(SeNB_crd[1]))

                        if(MEC_Network.readline() == 'SMD number:\n'): # 设置该SeNB中的任务数量
                            senb.SMDNumber = int(MEC_Network.readline())

                        for line1 in MEC_Network:
                            if (line1 == '---SeNB end---\n'):
                                break
                            elif(line1 == 'SMD:\n'):
                                self.totalSMDNumber += 1
                                smd = SMD()
                                if (MEC_Network.readline() == 'Coordinate:\n'): # 设置SMD坐标
                                    SMD_crd = MEC_Network.readline()
                                    SMD_crd = SMD_crd.splitlines()
                                    SMD_crd = SMD_crd[0].split('  ')
                                    smd.coordinate.append(float(SMD_crd[0]))
                                    smd.coordinate.append(float(SMD_crd[1]))

                                if (MEC_Network.readline() == 'Computation capacity:\n'): # 设置SMD3个core的计算能力
                                    SMD_cc = MEC_Network.readline()
                                    SMD_cc = SMD_cc.splitlines()
                                    SMD_cc = SMD_cc[0].split('  ') # core的计算能力
                                    smd.coreCC[1] = float(SMD_cc[0])
                                    smd.coreCC[2] = float(SMD_cc[1])
                                    smd.coreCC[3] = float(SMD_cc[2])

                                if (MEC_Network.readline() == 'The number of task:\n'):  #在SeNB（SeNB_count）下得到一个工作流
                                    taskNumber = int(MEC_Network.readline()) # 该SMD中的task数量
                                    SeNB_directory = "SeNB-"+str(SeNB_count)+"\\t"+str(taskNumber)+".txt"
                                    wf_directory = self.getCurrentPath()+"\workflowSet\\"+SeNB_directory
                                    smd.workflow = self.getWorkflow(wf_directory)
                                    smd.workflow.taskNumber = taskNumber
                                    self.codeLength += taskNumber
                                    for task in smd.workflow.taskSet: # 为每个任务的三元组设置信息
                                        task.c_i_j_k = float(file_SMD_task_cpu.readline())    # task需要的cpu周期数
                                        task.d_i_j_k = float(file_SMD_task_data.readline()) * 1024 # task输入数据大小
                                        task.o_i_j_k = float(file_SMD_output_task_data.readline()) * 1024 # task输出数据大小

                                if (MEC_Network.readline() == 'Channel:\n'): # 设置SMD所占据的通道索引
                                    channel = MEC_Network.readline()
                                    smd.channel = int(channel)

                                senb.SMDSet.append(smd) # 向该SeNB中添加该SMD
                    self.SeNBSet.append(senb) # 将该SeNB添加到MeNB中

        file_SMD_task_data.close()
        file_SMD_task_cpu.close()

    # 为每个SMD设置信道增益 I i,j
    def calculateInterference(self):
        for i in range(self.M):
            for j in range(self.SeNBSet[i].SMDNumber):
                I_i_j = 0
                for l in range(self.M):
                    if(self.SeNBSet[l] != self.SeNBSet[i]):
                        for k in range(self.SeNBSet[l].SMDNumber):
                            # 判断不同SeNB中的两个SMD所占用的channel是否相等
                            if(self.SeNBSet[l].SMDSet[k].channel == self.SeNBSet[i].SMDSet[j].channel):  # U_m_j and U_l_i have the same channel
                                g_i_l_k = self.getChannelGain(self.SeNBSet[l].SMDSet[k].coordinate,
                                                              self.SeNBSet[i].coordinate)
                                # 这里的λ设置的为1
                                I_i_j += self.SeNBSet[l].SMDSet[k].pws_i_j * g_i_l_k
                self.SeNBSet[i].SMDSet[j].I_i_j = I_i_j

    # 计算上行数据传输速率 R i,j
    def calculateDataTransmissionRate(self):
        self.calculateChannelGain()
        for i in range(self.M):
            for j in range(self.SeNBSet[i].SMDNumber):
                log_v = 1 + (self.SeNBSet[i].SMDSet[j].pws_i_j*self.SeNBSet[i].SMDSet[j].g_i_j) / (self.noisePower + self.SeNBSet[i].SMDSet[j].I_i_j)
                self.SeNBSet[i].SMDSet[j].R_i_j = self.w * math.log(log_v, 2)

    # 获得参考前沿
    def get_PF_ref(self):
        readReferPFPath = self.getCurrentPath() + "\ExperimentResult\\" + self.taskNumberRange + "\\referPF.xls"
        data = xlrd.open_workbook(readReferPFPath)
        table = data.sheet_by_name('total')
        return table._cell_values

    def getCurrentPath(self):
        return os.path.dirname(os.path.realpath(__file__))

    # channel gain= D^(-pl),
    # where D is the distance between U_m_j and S_m, pl=4 is the path loss factor
    # 传递过来的条件：处在不同SeNB的不同SMD具有相同channel
    # 传递来的参数为：一个SeNB中SMD的坐标、另一个SeNB的坐标
    # g i,(l,k)    p5页
    def getChannelGain(self, U_l_k, S_i):
        distance = self.getDistance(U_l_k, S_i)
        channelGain = pow(distance, -4)
        return channelGain

    # calculate G_m_j between SMD U_m_j and SeNB S_m
    def calculateChannelGain(self):
        for i in range(self.M):
            for j in range(self.SeNBSet[i].SMDNumber):
                self.SeNBSet[i].SMDSet[j].g_i_j = self.getChannelGain(self.SeNBSet[i].SMDSet[j].coordinate,
                                                                      self.SeNBSet[i].coordinate)

    # 获得欧氏距离
    def getDistance(self, point1, point2):
        return np.sqrt(np.sum(np.square([point1[i] - point2[i] for i in range(2)])))

    def getWorkflow(self, filename):
        wf = Workflow()
        with open(filename, 'r') as readFile:
            for line in readFile:
                task = Task()
                s = line.splitlines()
                s = s[0].split(':')
                predecessor = s[0]  # 前继任务
                id = s[1]           # 该任务id
                successor = s[2]    # 该任务的后继任务
                if (predecessor != ''):
                    predecessor = predecessor.split(',')
                    for pt in predecessor:
                        task.preTaskSet.append(int(pt))
                else:
                    wf.entryTask = int(id) # 设置任务流的开始任务为该task(前序任务为空的任务)
                task.id = int(id)

                if (successor != ''):
                    successor = successor.split(',')
                    for st in successor:
                        task.sucTaskSet.append(int(st))
                else:
                    wf.exitTask = int(id)

                wf.taskSet.append(task)
        return wf