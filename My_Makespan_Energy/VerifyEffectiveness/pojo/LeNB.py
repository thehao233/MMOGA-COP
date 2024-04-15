class LeNB:
    # def __deepcopy__(self, memodict={}):
    #     info = SeNB()
    #     info.coordinate = self.coordinate
    #     info.SMDNumber = self.SMDNumber
    #     info.SMDSet = self.SMDSet
    #     return info

    def __init__(self):
        self.coordinate = []  #The position coordination of the SeNB     SeNB的坐标
        self.SMDNumber = 0
        self.SMDSet = []    #The set of SMD the SeNB covers


