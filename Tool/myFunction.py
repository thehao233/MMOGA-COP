import numpy as np
import random

# 参考网址：https://blog.csdn.net/hujiaxuan1995/article/details/80629526
# size为rank个数，z为数据倾斜程度, 取值为0表示数据无倾斜，取值越大倾斜程度越高
def Zipf(size, z):
    p_list = [] # 累计概率列表
    div = 0 # 总概率
    for i in range(1, size+1):
        div += 1.0/pow(i, z)
    sum = 0
    for i in range(1, size+1):
        # The i-th probability in position i
        p = (1.0/pow(i, z))/div
        sum += p
        p_list.append(sum)

    rand = random.random() # 使用轮盘赌来选择任务
    for i in range(len(p_list)):
        if rand < p_list[i]:
            return i

