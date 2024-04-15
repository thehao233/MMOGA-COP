from My_Makespan_Energy.VerifyEffectiveness.pojo.Schedule import Schedule


class Workflow:

    # def __deepcopy__(self, memodict={}):
    #     info = Workflow()
    #     info.entryTask = self.entryTask
    #     info.exitTask = self.exitTask
    #     # info.position = self.position
    #     # info.sequence = self.sequence
    #     info.taskNumber = self.taskNumber
    #     info.taskSet = self.taskSet
    #     info.schedule = self.schedule
    #     return info


    def __init__(self):
        self.entryTask = None      #开始任务
        self.exitTask = None       #结束任务
        self.position = []         #执行位置
        self.sequence = []         #执行顺序
        self.taskNumber = None
        self.taskSet = []          #列表的索引值就是任务的id值
        self.schedule = Schedule()