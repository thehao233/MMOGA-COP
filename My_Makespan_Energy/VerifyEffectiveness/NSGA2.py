import random, copy, math, turtle, os, xlrd

# random.seed(20)
import pickle
import matplotlib.pyplot as plt
import numpy as np
import time
from Tool import myRandom
# 在NSGA2的基础上
# min energy,time 在初始化种群阶段
# DVFS

from My_Makespan_Energy.VerifyEffectiveness.pojo.Individual import Individual
from My_Makespan_Energy.VerifyEffectiveness.pojo.MeNB import MeNB
from My_Makespan_Energy.VerifyEffectiveness.pojo.Pareto import Pareto
from My_Makespan_Energy.VerifyEffectiveness.pojo.SMD import SMD
from My_Makespan_Energy.VerifyEffectiveness.pojo.LeNB import LeNB
from My_Makespan_Energy.VerifyEffectiveness.pojo.Task import Task
from My_Makespan_Energy.VerifyEffectiveness.pojo.Workflow import Workflow

class NSGA2:

    def run(self):

        # 随机初始化种群
        self.initializePopulation()
        self.ind0 = pickle.dumps(self.ind0)  # 序列化
        self.fast_non_dominated_sort(self.P_population)

        for i in range(1, self.F_rank.__len__()):
            self.crowding_distance_assignment(self.F_rank[i])
        self.update_EP_FromSet(self.EP, self.F_rank[1])  # Updating self.EP using self.F_rank[1]

        t = 1
        while (t <= self.maxGen):
            self.Q_population = self.make_new_population(self.P_population)  # 使用父代产生子代
            self.R_population = self.combine_Pt_and_Qt(self.P_population, self.Q_population)
            self.fast_non_dominated_sort(self.R_population)
            self.update_EP_FromSet(self.EP, self.F_rank[1])
            self.P_population = []
            i = 1
            while ((self.P_population.__len__() + self.F_rank[i].__len__()) <= self.popSize):
                self.P_population += self.F_rank[i]
                i += 1
            cha_N = self.popSize - self.P_population.__len__()
            if (cha_N > 0):  # 补足pSize个解
                self.crowding_distance_assignment(self.F_rank[i])  # 新加入的 !!
                self.F_rank[i] = sorted(self.F_rank[i], key=lambda Individual: Individual.distance,
                                        reverse=True)  # 按拥挤距离的降序排序
                self.P_population += self.F_rank[i][:cha_N]

            print('Generation ', t, ': ')

            t += 1

            import gc
            gc.collect()

        for ep in self.EP:
            ep.temp_fitness = ep.fitness[0]

        test_fast = sorted(self.EP, key=lambda Pareto: Pareto.temp_fitness)
        EP_list = [ind.fitness for ind in test_fast]
        return EP_list  # 返回最终的非支配解集

    def __init__(self, popSize, maxGen, pc, pm_SMD, pm_bit, taskNumberRange):

        self.menb = MeNB(taskNumberRange)  ## 新添加的，有错就改正
        self.smd = SMD()
        self.ind0 = Individual()

        self.popSize = popSize
        self.maxGen = maxGen
        self.pc = pc
        self.pm_SMD = pm_SMD
        self.pm_bit = pm_bit
        self.taskNumberRange = taskNumberRange  # 用于选取实例
        self.objectNumber = 2
        self.P_population = []  # 父代种群
        self.Q_population = []  # 子代种群
        self.R_population = []  # 临时种群

        self.Pop1 = []  # 层次遗传
        self.Pop2 = []

        self.F_rank = []  # 将种群非支配排序分层, 用种群中的个体的下标来表示，一个元素表示第一层,下标从1开始
        self.EP = []  # 保存当前代的历史非支配解
        self.PF_ref = []
        # self.PF_ref = self.menb.get_PF_ref()
        self.IGDValue = None
        self.IGD_list = []  # 保存100代的IGD值

    """
        **********************************************run**********************************************
    """

    def DVFS_Schedule_Algorithm(self, population):
        for ind in population:
            energy = []
            for gene in ind.chromosome:
                workflow = gene.workflow
                schedule = gene.workflow.schedule
                for coreId in schedule.S:
                    if coreId < 4 and schedule.S[coreId] != []:
                        for i in range(len(schedule.S[coreId])):
                            taskId = schedule.S[coreId][i]
                            vi = workflow.taskSet[taskId]
                            flag = 0
                            m = 1
                            while flag == 0 and m < self.smd.numOfCoreLevel:
                                FT_i_new = self.DVFS_calculateNewFinishTime(m, vi)
                                if i != (len(schedule.S[coreId]) - 1):  # There is next task v_j on the same core
                                    vj = workflow.taskSet[schedule.S[coreId][i + 1]]
                                    lim1 = vj.ST_i_l
                                else:
                                    lim1 = workflow.schedule.T_total
                                if vi.id != workflow.exitTask:
                                    lim2 = self.get_min_Succ_ST(workflow, vi)
                                else:
                                    lim2 = workflow.schedule.T_total
                                if FT_i_new <= lim1 and FT_i_new <= lim2:
                                    flag = 1
                                    vi.actualFre = m
                                    vi.FT_i_l = FT_i_new
                                    schedule.E_total -= vi.energy
                                    vi.energy = self.smd.coreLevel[vi.actualFre - 1] * vi.energy
                                    schedule.E_total += vi.energy
                                    schedule.TimeEnergy[1] = schedule.E_total
                                    ind.fitness = schedule.TimeEnergy
                                    break
                                m += 1
            #     energy.append(schedule.TimeEnergy[1])
            # ind.fitness[1] = np.mean(energy)

    def DVFS_calculateNewFinishTime(self, m, task):
        return task.ST_i_l + (task.FT_i_l - task.ST_i_l) / self.smd.coreLevel[m - 1]

    def get_min_Succ_ST(self, workflow, task):
        minST = []
        for succ_taskId in task.sucTaskSet:
            succ = workflow.taskSet[succ_taskId]
            if succ.islocal == True:
                minST.append(succ.ST_i_l)
            else:
                minST.append(succ.ST_i_ws)
        return min(minST)

    # 获得参考前沿
    def get_PF_ref(self):
        readReferPFPath = self.getCurrentPath() + "\ExperimentResult\\" + self.taskNumberRange + "\\referPF.xls"
        data = xlrd.open_workbook(readReferPFPath)
        table = data.sheet_by_name('total')
        return table._cell_values

    def getIGDValue(self, PF_ref, PF_know):
        sum = []
        for v in PF_ref:
            distance = self.d_v_PFSet(v, PF_know)
            sum.append(distance)
        return np.average(sum)

    # 与getIGDValue一同使用
    def d_v_PFSet(self, v, PFSet):  # 求v和PFSet中最近的距离
        dList = []
        for pf in PFSet:
            distance = self.getDistance(v, pf)
            dList.append(distance)
        return min(dList)

    # 获得欧氏距离
    def getDistance(self, point1, point2):
        return np.sqrt(np.sum(np.square([point1[i] - point2[i] for i in range(2)])))

    def update_EP_FromSet(self, EP, F_rank):  # 用当前代的非支配排序后的第一层的非支配解来更新EP
        if (EP == []):
            for ind in F_rank:
                ind_ = Individual()
                ind_.fitness = copy.deepcopy(ind.fitness)
                ind_.chromosome = copy.deepcopy(ind.chromosome)
                EP.append(ind_)
        else:
            for ind in F_rank:
                if (self.isExist(ind, EP) == False):  # 先判断ind是否在EP中，若在，则返回True。
                    if (self.isEP_Dominated_ind(ind, EP) == False):  # 然后再判断EP是否支配ind
                        i = 0
                        while (i < EP.__len__()):  # 判断ind是否支配EP中的非支配解，若支配，则删除它所支配的解
                            if (self.isDominated(ind.fitness, EP[i].fitness) == True):
                                EP.remove(EP[i])
                                i -= 1
                            i += 1
                        ind_ = Individual()
                        ind_.fitness = copy.deepcopy(ind.fitness)
                        ind_.chromosome = copy.deepcopy(ind.chromosome)
                        EP.append(ind_)

    def isExist(self, ind, EP):  # 判断个体ind的适应度是否与EP中某个个体的适应度相对，若相等，则返回True
        for ep in EP:
            if ind.fitness == ep.fitness:  # 判断两个列表对应元素的值是否相等
                return True
        return False

    def isEP_Dominated_ind(self, ind, EP):  # 判断EP中的某个个体是否支配ind，若支配，则返回True
        for ep in EP:
            if self.isDominated(ep.fitness, ind.fitness):
                return True
        return False

    def initializePopulation(self):
        for i in range(self.popSize):  # Generate randomly
            ind = Individual()
            for senb in self.menb.SeNBSet:
                for smd in senb.SMDSet:

                    # 保留smd集原始备份
                    if i == 0:
                        self.ind0.chromosome.append(copy.deepcopy(smd))

                    temp_smd = copy.deepcopy(smd)

                    # 为每个SMD中的每个task随机分配执行位置
                    for j in range(temp_smd.workflow.taskNumber):
                        temp_smd.workflow.position.append(random.randint(1, self.menb.H + 1))

                    # 为任务流中的task分配执行顺序
                    temp_smd.workflow.sequence = self.initializeWorkflowSequence(temp_smd.workflow)
                    ind.chromosome.append(temp_smd)  # 将分配好EL和EO的SMD作为Individual中的染色体

            self.calculateFitness(ind)
            self.P_population.append(ind)

    def fast_non_dominated_sort(self, population):
        for p in population:
            p.S_p = []
            p.rank = None
            p.n = 0

        self.F_rank = []
        F1 = []  # 第一个非支配解集前端
        self.F_rank.append(None)
        for p in population:
            for q in population:
                # 每次循环只设置与p有关的信息(支配的集合、被支配的个体的数量)
                if self.isDominated(p.fitness, q.fitness):
                    p.S_p.append(q)
                elif self.isDominated(q.fitness, p.fitness):
                    p.n += 1
            if (p.n == 0):  # 若p被支配数为0，则将p的支配等级设置为1
                p.rank = 1
                F1.append(p)
        self.F_rank.append(F1)

        i = 1
        while (self.F_rank[i] != []):
            Q = []
            for p in self.F_rank[i]:
                for q in p.S_p:
                    q.n -= 1
                    if (q.n == 0):
                        q.rank = i + 1
                        Q.append(q)

            if (Q != []):
                i += 1
                self.F_rank.append(Q)
            else:
                break

    def isDominated(self, fitness_1, fitness_2):  # 前者是否支配后者
        flag = -1
        for i in range(self.objectNumber):
            if fitness_1[i] < fitness_2[i]:
                flag = 0
            if fitness_1[i] > fitness_2[i]:
                return False
        if flag == 0:
            return True
        else:
            return False

    def crowding_distance_assignment(self, population):  # 计算拥挤距离
        for ind in population:  # initialize distance
            ind.distance = 0.0

        for m in range(self.objectNumber):
            for ind in population:
                ind.temp_fitness = ind.fitness[m]

            population = sorted(population, key=lambda
                Individual: Individual.temp_fitness)  # sort population by Individual.temp_fitness

            if (population.__len__() == 1):
                population[0].distance = float('inf')
            elif (population.__len__() == 2):
                population[0].distance = float('inf')
                population[1].distance = float('inf')
            else:
                population[0].distance = float('inf')
                population[-1].distance = float('inf')
                f_min = population[0].fitness[m]
                f_max = population[-1].fitness[m]
                for i in range(1, population.__len__() - 1):
                    population[i].distance += (population[i + 1].fitness[m] - population[i - 1].fitness[m]) / (
                            f_max - f_min + 1)  # 加一为了防止分母为零

    # 遗传操作产生新种群
    def make_new_population(self, population):
        new_population = []
        selectionResult = self.tournamentSelectionOperator(population)
        for ind in selectionResult:
            new_ind = pickle.loads(self.ind0)

            # for gene in ind.chromosome:

            for i in range(self.menb.totalSMDNumber):
                position = copy.deepcopy(ind.chromosome[i].workflow.position)
                sequence = copy.deepcopy(ind.chromosome[i].workflow.sequence)

                new_ind.chromosome[i].workflow.position = position
                new_ind.chromosome[i].workflow.sequence = sequence

            new_population.append(new_ind)
        self.crossoverOperator(new_population)
        self.mutantOperator(new_population)
        self.calculatePopulationFitness(new_population)
        return new_population

    # 用锦标赛选择法从父代中选择待进化子代
    # def tournamentSelectionOperator(self, population):
    #     selectionResult = []
    #     N = self.popSize
    #     while (N != 0):
    #         ind_1 = 0
    #         ind_2 = 0
    #         while (ind_1 == ind_2):
    #             ind_1 = random.randint(0, self.popSize - 1)
    #             ind_2 = random.randint(0, self.popSize - 1)
    #         if (population[ind_1].rank < population[ind_2].rank) or (  # 优先选取rank小的，distance大的
    #                 (population[ind_1].rank == population[ind_2].rank) and (
    #                 population[ind_1].distance > population[ind_2].distance)):
    #             selectionResult.append(population[ind_1])
    #         else:
    #             selectionResult.append(population[ind_2])
    #         N -= 1
    #     return selectionResult

    def tournamentSelectionOperator(self, population):
        selectionResult = []
        N = self.popSize
        while (N != 0):
            ind_1 = 0
            ind_2 = 0
            while (ind_1 == ind_2):
                ind_1 = random.randint(0, self.popSize - 1)
                ind_2 = random.randint(0, self.popSize - 1)


            if population[ind_1].fitness[0] < population[ind_2].fitness[0]:
                selectionResult.append(population[ind_1])

            else:
                selectionResult.append(population[ind_2])

            # if (population[ind_1].rank < population[ind_2].rank) or ( # 优先选取rank小的，distance大的
            #         (population[ind_1].rank == population[ind_2].rank) and (
            #         population[ind_1].distance > population[ind_2].distance)):
            #     selectionResult.append(population[ind_1])
            # else:
            #     selectionResult.append(population[ind_2])
            N -= 1
        return selectionResult


    # 交叉
    def crossoverOperator(self, population):
        random.shuffle(population)
        for i in range(0, len(population), 2):  # 每隔两个个体取一次
            ind_1 = population[i]
            ind_2 = population[i + 1]
            rnd = myRandom.get_0to1_RandomNumber()
            if (rnd < self.pc):
                for i in range(self.menb.totalSMDNumber):
                    gene_1 = ind_1.chromosome[i]
                    gene_2 = ind_2.chromosome[i]

                    # 获取随机交叉点
                    cpt = random.randint(0, len(gene_1.workflow.position) - 1)

                    # EL交叉
                    p1 = gene_2.workflow.position[0:cpt]
                    p2 = gene_1.workflow.position[0:cpt]
                    gene_1.workflow.position = p1 + gene_1.workflow.position[cpt:]
                    gene_2.workflow.position = p2 + gene_2.workflow.position[cpt:]

                    # EO交叉
                    sequence1 = gene_2.workflow.sequence[0:cpt] + gene_1.workflow.sequence
                    sequence2 = gene_1.workflow.sequence[0:cpt] + gene_2.workflow.sequence
                    gene_1.workflow.sequence = sorted(list(set(sequence1)), key=sequence1.index)
                    gene_2.workflow.sequence = sorted(list(set(sequence2)), key=sequence2.index)

    # 变异
    def mutantOperator(self, population):
        for ind in population:
            for gene in ind.chromosome:
                rnd_SMD = myRandom.get_0to1_RandomNumber()
                if (rnd_SMD < 1.0 / self.menb.totalSMDNumber):  # 针对每一个基因（SMD）判断是否变异
                    for i in range(len(gene.workflow.position)):  # EL变异
                        rnd_bit = myRandom.get_0to1_RandomNumber()
                        if (rnd_bit < 1.0 / (gene.workflow.taskNumber)):
                            pos = gene.workflow.position[i]
                            rand = [1, 2, 3, 4]
                            rand.remove(pos)
                            # 从除当前EL外的另外三个EL中任选一个
                            gene.workflow.position[i] = random.choice(rand)

                    wrap_list = []  # save two adjacent task id without precedence relationship
                    for j in range(1, len(gene.workflow.sequence) - 1):  # EO变异
                        first_taskId = gene.workflow.sequence[j]
                        second_taskId = gene.workflow.sequence[j + 1]
                        if first_taskId not in gene.workflow.taskSet[second_taskId].preTaskSet:
                            wrap_list.append((j, j + 1))
                    if wrap_list != []:
                        index = np.random.randint(0, len(wrap_list))
                        first_taskIndex = wrap_list[index][0]
                        second_taskIndex = wrap_list[index][1]
                        gene.workflow.sequence[first_taskIndex], gene.workflow.sequence[second_taskIndex] = \
                            gene.workflow.sequence[second_taskIndex], gene.workflow.sequence[first_taskIndex]
                    else:
                        gene.workflow.sequence = self.initializeWorkflowSequence(gene.workflow)

    # 组合两个种群
    def combine_Pt_and_Qt(self, P_population, Q_population):
        # population = []
        # for p in P_population:
        #     population.append(p)
        # for p in Q_population:
        #     population.append(p)
        return P_population + Q_population

    # 把pop2加到pop1中
    def combine_pop1_and_pop2(self, pop1, pop2):
        # for p in pop2:
        #     pop1.append(p)
        pop1 = pop1 + pop2

    # 计算种群适应度值
    def calculatePopulationFitness(self, population):
        for ind in population:
            self.calculateFitness(ind)

    # 平均能耗和平均执行时间
    def calculateFitness(self, ind):
        W = []
        Q = []
        for workflowIndex in range(len(ind.chromosome)):

            smd = ind.chromosome[workflowIndex]
            workflow = smd.workflow

            for i in range(len(workflow.sequence)):
                taskId = workflow.sequence[i]
                pos = workflow.position[i]
                task = workflow.taskSet[taskId]
                task.exePosition = pos  # 设置task的执行位置
                T = task.c_i_j_k / smd.coreCC[pos]

                if pos == 4:  # The task is executed on [10,20] server. 在MEC执行
                    task.islocal = False
                    if task.id == workflow.entryTask:  # 若该task为工作流的起始任务
                        task.RT_i_l = task.ST_i_l = task.FT_i_l = 0
                        task.RT_i_ws = task.ST_i_ws = 0.0  # 准备好可以进行传输的时间 和 开始进行传输的时间
                        task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j  # task完成传输的时间
                        task.RT_i_c = task.ST_i_c = task.FT_i_ws  # 任务准备好进行计算的时间 和 开始进行计算的时间
                        task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]  # task完成执行的时间
                        task.RT_i_wr = task.ST_i_wr = task.FT_i_c  # task准备好返回数据的时间 和 task开始进行返回数据的时间
                        task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j  # task完成数据返回的时间
                        workflow.schedule.wsTP.append(task.FT_i_ws)
                        workflow.schedule.MECTP.append(task.FT_i_c)
                        workflow.schedule.wrTP.append(task.FT_i_wr)
                    else:  # 该task为中间task
                        task.RT_i_ws = self.get_RT_i_ws(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                        task.ST_i_l = float("inf")
                        task.FT_i_l = float("inf")
                        if workflow.schedule.wsTP[-1] < task.RT_i_ws:
                            task.ST_i_ws = task.RT_i_ws
                            task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                        else:
                            task.ST_i_ws = workflow.schedule.wsTP[-1]
                            task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                        workflow.schedule.wsTP.append(task.FT_i_ws)

                        task.RT_i_c = self.get_RT_i_c(task, workflow)
                        if workflow.schedule.MECTP[-1] < task.RT_i_c:
                            task.ST_i_c = task.RT_i_c
                            task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                        else:
                            task.ST_i_c = workflow.schedule.MECTP[-1]
                            task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                        T = task.c_i_j_k / smd.coreCC[pos]
                        workflow.schedule.MECTP.append(task.FT_i_c)

                        task.RT_i_wr = task.FT_i_c
                        if workflow.schedule.wrTP[-1] < task.RT_i_wr:
                            task.ST_i_wr = task.RT_i_wr
                            task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                        else:
                            task.ST_i_wr = workflow.schedule.wrTP[-1]
                            task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                        workflow.schedule.wrTP.append(task.FT_i_wr)
                    task.energy += smd.pws_i_j * (task.FT_i_ws - task.ST_i_ws)  # 提交能耗
                    task.energy += smd.pwr_i_j * (task.FT_i_wr - task.ST_i_wr)  # 回传能耗
                    workflow.schedule.E_total += task.energy
                    W.append([task, workflowIndex, i])
                    break
                else:  # The task is executed on a local core.
                    task.islocal = True
                    task.RT_i_ws = task.RT_i_c = task.RT_i_wr = 0.0
                    task.ST_i_ws = task.ST_i_c = task.ST_i_wr = 0.0
                    task.FT_i_ws = task.FT_i_c = task.FT_i_wr = 0.0

                    T = task.c_i_j_k / smd.coreCC[pos]
                    if task.id == workflow.entryTask:
                        task.RT_i_l = task.ST_i_l = 0
                        task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                    else:  # 该task为中间task
                        task.RT_i_l = self.get_RT_i_l(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                        if task.RT_i_l > workflow.schedule.coreTP[pos][-1]:
                            task.ST_i_l = task.RT_i_l
                        else:
                            task.ST_i_l = workflow.schedule.coreTP[pos][-1]
                        task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                    workflow.schedule.coreTP[pos].append(task.FT_i_l)
                    task.energy = smd.pcc_i_j[pos] * (task.FT_i_l - task.ST_i_l)
                    workflow.schedule.E_total += task.energy
                workflow.schedule.S[pos].append(task.id)

        while len(W) != 0:
            temp = W[0]
            minCT_send = temp[0].FT_i_ws
            minCT_sendIndex = 0
            for i in range(len(W)):
                temp = W[i]
                temp_task = temp[0]
                if temp_task.FT_i_ws < minCT_send:
                    minCT_sendIndex = i
                    minCT_send = temp_task.FT_i_ws
            m, workflowIndex, sequence = W[minCT_sendIndex]
            W.remove(W[minCT_sendIndex])

            Q_ = [t for t in Q if t.FT_i_c > m.FT_i_ws]
            Q = Q_

            if(len(Q) == self.N_VM):

                earlistTask = self.earlistTaskInQ(Q)
                m.RT_i_c = earlistTask.FT_i_c
                Q.remove(earlistTask)
                Q.append(m)

                smd = ind.chromosome[workflowIndex]
                workflow = smd.workflow
                workflow.schedule.MECTP.append(earlistTask.FT_i_c)

                for i in range(sequence, len(workflow.sequence)):
                    taskId = workflow.sequence[i]
                    pos = workflow.position[i]
                    task = workflow.taskSet[taskId]
                    task.exePosition = pos  # 设置task的执行位置

                    if pos == 4:  # The task is executed on [10,20] server. 在MEC执行
                        task.islocal = False
                        if task.id == workflow.entryTask:  # 若该task为工作流的起始任务
                            task.RT_i_l = task.ST_i_l = task.FT_i_l = 0
                            task.RT_i_ws = task.ST_i_ws = 0.0  # 准备好可以进行传输的时间 和 开始进行传输的时间
                            task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j  # task完成传输的时间
                            task.RT_i_c = task.ST_i_c = task.FT_i_ws  # 任务准备好进行计算的时间 和 开始进行计算的时间
                            task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]  # task完成执行的时间
                            task.RT_i_wr = task.ST_i_wr = task.FT_i_c  # task准备好返回数据的时间 和 task开始进行返回数据的时间
                            task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j  # task完成数据返回的时间
                            workflow.schedule.wsTP.append(task.FT_i_ws)
                            workflow.schedule.MECTP.append(task.FT_i_c)
                            workflow.schedule.wrTP.append(task.FT_i_wr)
                        else:  # 该task为中间task
                            task.RT_i_ws = self.get_RT_i_ws(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                            task.ST_i_l = float("inf")
                            task.FT_i_l = float("inf")
                            if workflow.schedule.wsTP[-1] < task.RT_i_ws:
                                task.ST_i_ws = task.RT_i_ws
                                task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                            else:
                                task.ST_i_ws = workflow.schedule.wsTP[-1]
                                task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                            workflow.schedule.wsTP.append(task.FT_i_ws)

                            task.RT_i_c = self.get_RT_i_c(task, workflow)
                            if workflow.schedule.MECTP[-1] < task.RT_i_c:
                                task.ST_i_c = task.RT_i_c
                                task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                            else:
                                task.ST_i_c = workflow.schedule.MECTP[-1]
                                task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                            workflow.schedule.MECTP.append(task.FT_i_c)

                            task.RT_i_wr = task.FT_i_c
                            if workflow.schedule.wrTP[-1] < task.RT_i_wr:
                                task.ST_i_wr = task.RT_i_wr
                                task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                            else:
                                task.ST_i_wr = workflow.schedule.wrTP[-1]
                                task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                            workflow.schedule.wrTP.append(task.FT_i_wr)
                        task.energy += smd.pws_i_j * (task.FT_i_ws - task.ST_i_ws)  # 提交能耗
                        task.energy += smd.pwr_i_j * (task.FT_i_wr - task.ST_i_wr)  # 回传能耗
                        workflow.schedule.E_total += task.energy
                        if i != sequence :
                            W.append([task, workflowIndex, i])
                            break
                    else:  # The task is executed on a local core.
                        task.islocal = True
                        task.RT_i_ws = task.RT_i_c = task.RT_i_wr = 0.0
                        task.ST_i_ws = task.ST_i_c = task.ST_i_wr = 0.0
                        task.FT_i_ws = task.FT_i_c = task.FT_i_wr = 0.0
                        if task.id == workflow.entryTask:
                            task.RT_i_l = task.ST_i_l = 0
                            task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                        else:  # 该task为中间task
                            task.RT_i_l = self.get_RT_i_l(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                            if task.RT_i_l > workflow.schedule.coreTP[pos][-1]:
                                task.ST_i_l = task.RT_i_l
                            else:
                                task.ST_i_l = workflow.schedule.coreTP[pos][-1]
                            task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                        workflow.schedule.coreTP[pos].append(task.FT_i_l)
                        task.energy = smd.pcc_i_j[pos] * (task.FT_i_l - task.ST_i_l)
                        workflow.schedule.E_total += task.energy
                    workflow.schedule.S[pos].append(task.id)

            else:
                Q.append(m)
                smd = ind.chromosome[workflowIndex]
                workflow = smd.workflow

                for i in range(sequence+1, len(workflow.sequence)):
                    taskId = workflow.sequence[i]
                    pos = workflow.position[i]
                    task = workflow.taskSet[taskId]
                    task.exePosition = pos  # 设置task的执行位置

                    if pos == 4:  # The task is executed on [10,20] server. 在MEC执行
                        task.islocal = False
                        if task.id == workflow.entryTask:  # 若该task为工作流的起始任务
                            task.RT_i_l = task.ST_i_l = task.FT_i_l = 0
                            task.RT_i_ws = task.ST_i_ws = 0.0  # 准备好可以进行传输的时间 和 开始进行传输的时间
                            task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j  # task完成传输的时间
                            task.RT_i_c = task.ST_i_c = task.FT_i_ws  # 任务准备好进行计算的时间 和 开始进行计算的时间
                            task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]  # task完成执行的时间
                            task.RT_i_wr = task.ST_i_wr = task.FT_i_c  # task准备好返回数据的时间 和 task开始进行返回数据的时间
                            task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j  # task完成数据返回的时间
                            workflow.schedule.wsTP.append(task.FT_i_ws)
                            workflow.schedule.MECTP.append(task.FT_i_c)
                            workflow.schedule.wrTP.append(task.FT_i_wr)
                        else:  # 该task为中间task
                            task.RT_i_ws = self.get_RT_i_ws(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                            task.ST_i_l = float("inf")
                            task.FT_i_l = float("inf")
                            if workflow.schedule.wsTP[-1] < task.RT_i_ws:
                                task.ST_i_ws = task.RT_i_ws
                                task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                            else:
                                task.ST_i_ws = workflow.schedule.wsTP[-1]
                                task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                            workflow.schedule.wsTP.append(task.FT_i_ws)

                            task.RT_i_c = self.get_RT_i_c(task, workflow)
                            if workflow.schedule.MECTP[-1] < task.RT_i_c:
                                task.ST_i_c = task.RT_i_c
                                task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                            else:
                                task.ST_i_c = workflow.schedule.MECTP[-1]
                                task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                            workflow.schedule.MECTP.append(task.FT_i_c)

                            task.RT_i_wr = task.FT_i_c
                            if workflow.schedule.wrTP[-1] < task.RT_i_wr:
                                task.ST_i_wr = task.RT_i_wr
                                task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                            else:
                                task.ST_i_wr = workflow.schedule.wrTP[-1]
                                task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                            workflow.schedule.wrTP.append(task.FT_i_wr)
                        task.energy += smd.pws_i_j * (task.FT_i_ws - task.ST_i_ws)  # 提交能耗
                        task.energy += smd.pwr_i_j * (task.FT_i_wr - task.ST_i_wr)  # 回传能耗
                        workflow.schedule.E_total += task.energy
                        # if i != sequence + 1:
                        W.append([task, workflowIndex, i])
                        break

                    else:  # The task is executed on a local core.
                        task.islocal = True
                        task.RT_i_ws = task.RT_i_c = task.RT_i_wr = 0.0
                        task.ST_i_ws = task.ST_i_c = task.ST_i_wr = 0.0
                        task.FT_i_ws = task.FT_i_c = task.FT_i_wr = 0.0
                        if task.id == workflow.entryTask:
                            task.RT_i_l = task.ST_i_l = 0
                            task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                        else:  # 该task为中间task
                            task.RT_i_l = self.get_RT_i_l(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                            if task.RT_i_l > workflow.schedule.coreTP[pos][-1]:
                                task.ST_i_l = task.RT_i_l
                            else:
                                task.ST_i_l = workflow.schedule.coreTP[pos][-1]
                            task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                        workflow.schedule.coreTP[pos].append(task.FT_i_l)
                        task.energy = smd.pcc_i_j[pos] * (task.FT_i_l - task.ST_i_l)
                        workflow.schedule.E_total += task.energy
                    workflow.schedule.S[pos].append(task.id)

            if workflow.taskSet[workflow.exitTask].islocal == True:
                # 若终止task在本地执行，则整个任务流的完成时间为该终止task的完成执行时间
                workflow.schedule.T_total = workflow.taskSet[workflow.exitTask].FT_i_l
            else:
                # 否则整个任务流的完成时间为该终止task的完成数据返回的时间
                workflow.schedule.T_total = workflow.taskSet[workflow.exitTask].FT_i_wr
            workflow.schedule.TimeEnergy.append(workflow.schedule.T_total)
            workflow.schedule.TimeEnergy.append(workflow.schedule.E_total)

        ind.fitness = []
        time = []  # 记录每个workflow的执行时间
        energy = []  # 记录每个workflow的执行能耗
        for gene in ind.chromosome:
            smd = gene
            self.calculateWorkflowTimeEnergy(smd, smd.workflow)
            if smd.workflow.schedule.T_total == None or smd.workflow.schedule.E_total == None:
                self.calculateWorkflowTimeEnergy(smd, smd.workflow)
            time.append(smd.workflow.schedule.T_total)
            energy.append(smd.workflow.schedule.E_total)

        ind.fitness.append(np.average(time))
        ind.fitness.append(np.average(energy))

    def earlistTaskInQ(self, Q):
        earlistTask = Q[0]
        for t in Q:
            if t.FT_i_c < earlistTask.FT_i_c:
                earlistTask = t
        return earlistTask


    # 计算一条工作流的能耗和时延
    def calculateWorkflowTimeEnergy(self, smd, workflow):
        workflow.schedule.TimeEnergy = []
        workflow.schedule.T_total = None
        workflow.schedule.E_total = 0

        for i in range(len(workflow.sequence)):
            taskId = workflow.sequence[i]
            pos = workflow.position[i]
            task = workflow.taskSet[taskId]
            task.exePosition = pos  # 设置task的执行位置
            if pos == 4:  # The task is executed on [10,20] server. 在MEC执行
                task.islocal = False
                if task.id == workflow.entryTask:  # 若该task为工作流的起始任务
                    task.RT_i_l = task.ST_i_l = task.FT_i_l = 0
                    task.RT_i_ws = task.ST_i_ws = 0.0  # 准备好可以进行传输的时间 和 开始进行传输的时间
                    task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j  # task完成传输的时间
                    task.RT_i_c = task.ST_i_c = task.FT_i_ws  # 任务准备好进行计算的时间 和 开始进行计算的时间
                    task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]  # task完成执行的时间
                    task.RT_i_wr = task.ST_i_wr = task.FT_i_c  # task准备好返回数据的时间 和 task开始进行返回数据的时间
                    task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j  # task完成数据返回的时间
                    workflow.schedule.wsTP.append(task.FT_i_ws)
                    workflow.schedule.MECTP.append(task.FT_i_c)
                    workflow.schedule.wrTP.append(task.FT_i_wr)
                else:  # 该task为中间task
                    task.RT_i_ws = self.get_RT_i_ws(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                    task.ST_i_l = float("inf")
                    task.FT_i_l = float("inf")
                    if workflow.schedule.wsTP[-1] < task.RT_i_ws:
                        task.ST_i_ws = task.RT_i_ws
                        task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                    else:
                        task.ST_i_ws = workflow.schedule.wsTP[-1]
                        task.FT_i_ws = task.ST_i_ws + task.d_i_j_k / smd.R_i_j
                    workflow.schedule.wsTP.append(task.FT_i_ws)

                    task.RT_i_c = self.get_RT_i_c(task, workflow)
                    if workflow.schedule.MECTP[-1] < task.RT_i_c:
                        task.ST_i_c = task.RT_i_c
                        task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                    else:
                        task.ST_i_c = workflow.schedule.MECTP[-1]
                        task.FT_i_c = task.ST_i_c + task.c_i_j_k / smd.coreCC[pos]
                    workflow.schedule.MECTP.append(task.FT_i_c)

                    task.RT_i_wr = task.FT_i_c
                    if workflow.schedule.wrTP[-1] < task.RT_i_wr:
                        task.ST_i_wr = task.RT_i_wr
                        task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                    else:
                        task.ST_i_wr = workflow.schedule.wrTP[-1]
                        task.FT_i_wr = task.ST_i_wr + task.o_i_j_k / smd.R_i_j
                    workflow.schedule.wrTP.append(task.FT_i_wr)
                task.energy += smd.pws_i_j * (task.FT_i_ws - task.ST_i_ws)  # 提交能耗
                task.energy += smd.pwr_i_j * (task.FT_i_wr - task.ST_i_wr)  # 回传能耗
                workflow.schedule.E_total += task.energy
            else:  # The task is executed on a local core.
                task.islocal = True
                task.RT_i_ws = task.RT_i_c = task.RT_i_wr = 0.0
                task.ST_i_ws = task.ST_i_c = task.ST_i_wr = 0.0
                task.FT_i_ws = task.FT_i_c = task.FT_i_wr = 0.0
                if task.id == workflow.entryTask:
                    task.RT_i_l = task.ST_i_l = 0
                    task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                else:  # 该task为中间task
                    task.RT_i_l = self.get_RT_i_l(task, workflow)  # 设置该task的RT为其前序任务的完成时间
                    if task.RT_i_l > workflow.schedule.coreTP[pos][-1]:
                        task.ST_i_l = task.RT_i_l
                    else:
                        task.ST_i_l = workflow.schedule.coreTP[pos][-1]
                    task.FT_i_l = task.ST_i_l + task.c_i_j_k / smd.coreCC[pos]
                workflow.schedule.coreTP[pos].append(task.FT_i_l)
                task.energy = smd.pcc_i_j[pos] * (task.FT_i_l - task.ST_i_l)
                workflow.schedule.E_total += task.energy
            workflow.schedule.S[pos].append(task.id)

        if workflow.taskSet[workflow.exitTask].islocal == True:
            # 若终止task在本地执行，则整个任务流的完成时间为该终止task的完成执行时间
            workflow.schedule.T_total = workflow.taskSet[workflow.exitTask].FT_i_l
        else:
            # 否则整个任务流的完成时间为该终止task的完成数据返回的时间
            workflow.schedule.T_total = workflow.taskSet[workflow.exitTask].FT_i_wr
        workflow.schedule.TimeEnergy.append(workflow.schedule.T_total)
        workflow.schedule.TimeEnergy.append(workflow.schedule.E_total)

    # 设置该task的RT为其前序任务的完成时间
    def get_RT_i_ws(self, task, workflow):
        if task.id == workflow.entryTask:
            return 0.0
        else:
            pre = []  # 该task的前序task的完成时间
            for pre_taskId in task.preTaskSet:
                if workflow.taskSet[pre_taskId].islocal == True:
                    # 如果前序任务在local执行，则添加前序任务的执行完成时间
                    pre.append(workflow.taskSet[pre_taskId].FT_i_l)
                else:
                    # 如果前序任务在MEC执行，则添加前序任务的回传完成时间
                    pre.append(workflow.taskSet[pre_taskId].FT_i_ws)
            return max(pre)

    def get_RT_i_c(self, task, workflow):
        pre_max = []
        for pre_taskId in task.preTaskSet:
            pre_max.append(workflow.taskSet[pre_taskId].FT_i_c)
        return max(task.FT_i_ws, max(pre_max))

    def get_RT_i_l(self, task, workflow):
        if task.id == workflow.entryTask:
            return 0.0
        else:
            pre_max = []
            for pre_taskId in task.preTaskSet:
                if workflow.taskSet[pre_taskId].islocal == True:
                    pre_max.append(workflow.taskSet[pre_taskId].FT_i_l)
                else:
                    pre_max.append(workflow.taskSet[pre_taskId].FT_i_wr)
            return max(pre_max)

    def reInitialize_WorkflowTaskSet_Schedule(self, smd):
        for task in smd.workflow.taskSet:
            self.reInitializeTaskSet(task)
        self.reInitializeSchedule(smd.workflow.schedule)

    def reInitializeTaskSet(self, task):
        task.islocal = None
        task.exePosition = None
        task.RT_i_l = task.ST_i_l = task.FT_i_l = None
        task.RT_i_ws = task.RT_i_c = task.RT_i_wr = None
        task.ST_i_ws = task.ST_i_c = task.ST_i_wr = None
        task.FT_i_ws = task.FT_i_c = task.FT_i_wr = None
        task.energy = 0

    def reInitializeSchedule(self, schedule):
        schedule.S = {1: [], 2: [], 3: [], 4: []}
        schedule.coreTP = {1: [0], 2: [0], 3: [0]}
        schedule.wsTP = [0]
        schedule.MECTP = [0]
        schedule.wrTP = [0]
        schedule.T_total = None
        schedule.E_total = 0
        schedule.TimeEnergy = []

    def initializeWorkflowSequence(self, workflow):
        # 这个方法应该是应用了 RSEOI
        S = []  # 前序任务已经被排序的任务
        R = []  # 已排序任务
        T = []  # 待排序的任务集合
        R.append(workflow.entryTask)
        for task in workflow.taskSet:
            T.append(task.id)
        T.remove(workflow.entryTask)

        while T != []:
            for t in T:
                # 判断t的前驱节点集是否包含在R中
                if set(workflow.taskSet[t].preTaskSet).issubset(set(R)):
                    if t not in S:
                        S.append(t)
            ti = random.choice(S)  # 随机从S中选择一个元素
            S.remove(ti)
            T.remove(ti)
            R.append(ti)
        return R

    def printWorkflowSchedule(self, workflow):
        print("S--", workflow.schedule.S)
        for coreId in workflow.schedule.S:
            if coreId < 4:
                print("core", coreId, ": ", end="")
                for taskId in workflow.schedule.S[coreId]:
                    task = workflow.taskSet[taskId]
                    print(str(taskId) + "=(" + str(round(task.ST_i_l, 2)) + "," + str(round(task.FT_i_l, 2)) + ") ",
                          end="")
                print("\n")

        for coreId in workflow.schedule.S:
            if coreId == 4:
                print("WS:      ", end="")
                for taskId in workflow.schedule.S[coreId]:
                    task = workflow.taskSet[taskId]
                    print(str(taskId) + "=(" + str(round(task.ST_i_ws, 2)) + "," + str(round(task.FT_i_ws, 2)) + ") ",
                          end="")
                break
        print("\n")

        for coreId in workflow.schedule.S:
            if coreId == 4:
                print("Cloud:   ", end="")
                for taskId in workflow.schedule.S[coreId]:
                    task = workflow.taskSet[taskId]
                    print(str(taskId) + "=(" + str(round(task.ST_i_c, 2)) + "," + str(round(task.FT_i_c, 2)) + ") ",
                          end="")
                break
        print("\n")

        for coreId in workflow.schedule.S:
            if coreId == 4:
                print("WR:      ", end="")
                for taskId in workflow.schedule.S[coreId]:
                    task = workflow.taskSet[taskId]
                    print(str(taskId) + "=(" + str(round(task.ST_i_wr, 2)) + "," + str(round(task.FT_i_wr, 2)) + ") ",
                          end="")
                break
        print("\n")
        print("(Time, Energy)=", workflow.schedule.TimeEnergy)
        print("\n\n")

    # channel gain= D^(-pl),
    # where D is the distance between U_m_j and S_m, pl=4 is the path loss factor
    # 传递过来的条件：处在不同SeNB的不同SMD具有相同channel
    # 传递来的参数为：一个SeNB中SMD的坐标、另一个SeNB的坐标
    # g i,(l,k)    p5页
    def getChannelGain(self, U_l_k, S_i):
        distance = self.getDistance(U_l_k, S_i)
        channelGain = pow(distance, -4)
        return channelGain

    def getWorkflow(self, filename):
        wf = Workflow()
        with open(filename, 'r') as readFile:
            for line in readFile:
                task = Task()
                s = line.splitlines()
                s = s[0].split(':')
                predecessor = s[0]  # 前继任务
                id = s[1]  # 该任务id
                successor = s[2]  # 该任务的后继任务
                if (predecessor != ''):
                    predecessor = predecessor.split(',')
                    for pt in predecessor:
                        task.preTaskSet.append(int(pt))
                else:
                    wf.entryTask = int(id)  # 设置任务流的开始任务为该task(前序任务为空的任务)
                task.id = int(id)

                if (successor != ''):
                    successor = successor.split(',')
                    for st in successor:
                        task.sucTaskSet.append(int(st))
                else:
                    wf.exitTask = int(id)

                wf.taskSet.append(task)
        return wf

    def printPopulationFitness(self, population):
        for ind in population:
            print('Index:  ', population.index(ind), "--", ind.fitness)

    def getCurrentPath(self):
        return os.path.dirname(os.path.realpath(__file__))

    def getProjectPath(self):
        cur_path = os.path.dirname(os.path.realpath(__file__))
        return os.path.join(os.path.dirname(cur_path))

    def calculateExeOnEdgeNum(self, ind):
        for smd in ind.chromosome:
            position = smd.workflow.position
            for pos in position:
                if pos == 4:
                    ind.exe_on_edge += 1
        pass


if __name__ == "__main__":
    popSize = 100
    maxGen = 100
    pc = 0.8  # 交叉概率
    pm_SMD = 0.3  #
    pm_bit = 0.01

    # taskNumberRange= '[10,20]'
    # taskNumberRange = '[20,30]'
    taskNumberRange = '[30,40]'
    # taskNumberRange = '[40,50]'
    # taskNumberRange = '[50,60]'
    # taskNumberRange = '[60,70]'
    # taskNumberRange = '[70,80]'


    nsga2 = NSGA2(popSize, maxGen, pc, pm_SMD, pm_bit, taskNumberRange)
    nsga2.run()
