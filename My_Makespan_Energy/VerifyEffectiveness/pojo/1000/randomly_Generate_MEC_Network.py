# coding=utf-8
import turtle
import time, math, random
from Tool import myRandom
import numpy as np


def randomlyGenerateMECNetwork(MEC_radius, SeNB_radius, SeNB_number, SMD_range, channelNumber, taskNumber_range):
    SMD_N = 1
    rotateDegree = (2*math.pi)/SeNB_number
    sum_rotateDegree = rotateDegree
    SeNBCoordinateSet = []
    while(SMD_N<=SeNB_number):
        SeNBCoordinate = []
        x = MEC_radius * math.cos(sum_rotateDegree)
        y = MEC_radius * math.sin(sum_rotateDegree)
        SeNBCoordinate.append(x)
        SeNBCoordinate.append(y)
        SeNBCoordinateSet.append(SeNBCoordinate)
        sum_rotateDegree += rotateDegree
        SMD_N += 1

    turtle.setup(width=1000, height=1000, startx=10, starty=10)
    turtle.pencolor('red')
    turtle.dot(18)
    turtle.penup()
    turtle.goto(0, -MEC_radius)
    turtle.pendown()
    turtle.pensize(3)
    turtle.pencolor('black')
    turtle.speed(50)
    turtle.circle(MEC_radius)

    turtle.pencolor('blue')

    MEC_WF = open('MEC_Network.txt', 'w')
    TaskCPUCycleNumber_WF = open('SMD_Task_CPU_Cycles_Number.txt', 'w')
    TaskDataSize_WF = open('SMD_Task_Data_Size.txt', 'w')
    outputTaskDataSize_WF = open('SMD_Task_Output_Data_Size.txt', 'w')

    for crd in SeNBCoordinateSet:
        turtle.pencolor('blue')
        x = crd[0]
        y = crd[1]
        if (abs(x)<0.000001):  # 表示crd的横坐标为0
            y = y/2
        elif (abs(y)<0.000001):
            x = x/2
        else:
            x = x/2
            y = y/2
        crd[0] = x
        crd[1] = y

        turtle.penup()
        turtle.goto(x, y)
        turtle.pendown()
        turtle.dot(15)

        turtle.penup()
        turtle.goto(x, y-SeNB_radius)
        turtle.pendown()
        turtle.circle(SeNB_radius)


        MEC_WF.write('SeNB:\n')
        MEC_WF.write('Coordinate:\n')
        MEC_WF.write(str(crd[0])+'   '+str(crd[1])+'\n')


        turtle.pencolor('green')
        SMD_N = random.randint(SMD_range[0], SMD_range[1])
        MEC_WF.write('SMD number:\n')
        MEC_WF.write(str(SMD_N)+'\n')
        x1 = crd[0]
        x2 = crd[1]
        radius = SeNB_radius
        a = 2 * math.pi * np.array([random.random() for _ in range(SMD_N)])
        r = np.array([random.random() for _ in range(SMD_N)])
        x = radius * np.sqrt(r) * np.cos(a) + x1
        y = radius * np.sqrt(r) * np.sin(a) + x2

        channelList = [i for i in range(1, channelNumber + 1)]   #保证同一个SeNB中的每个SMD的信道不同
        taskNumberList = []
        while len(taskNumberList) < SMD_N:
            for i in range(taskNumber_range[0], taskNumber_range[1] + 1):
                taskNumberList.append(i)
        # taskNumberList = [i for i in range(taskNumber_range[0], taskNumber_range[1]+1)]   #保证同一个SeNB中的每个SMD的任务数不同
        # for i in range(taskNumber_range[0], taskNumber_range[1]+1):
        #     taskNumberList.append(i)
        for j in range(SMD_N):
            MEC_WF.write('SMD:\n')
            MEC_WF.write('Coordinate:\n')
            MEC_WF.write(str(x[j]) + '   ' + str(y[j]) + '\n')
            turtle.penup()
            turtle.goto(x[j], y[j])
            turtle.pendown()
            turtle.dot(10)

            MEC_WF.write('Computation capacity:\n')
            cc = random.uniform(0.5, 1)
            MEC_WF.write(str(cc) + '   '+str(cc-0.1)+'   '+str(cc-0.25)+'\n')

            MEC_WF.write('The number of task:\n')
            taskNumber = random.choice(taskNumberList)
            taskNumberList.remove(taskNumber)
            MEC_WF.write(str(taskNumber) + '\n')

            channel = random.choice(channelList)
            channelList.remove(channel)
            MEC_WF.write('Channel:\n')
            MEC_WF.write(str(channel) + '\n')

            MEC_WF.write('\n\n\n')
            for i in range(taskNumber):
                TaskCPUCycleNumber_WF.write(str(random.uniform(0.1, 0.5))+'\n')
                TaskDataSize_WF.write(str(random.uniform(5000, 6000))+'\n')
                outputTaskDataSize_WF.write(str(random.uniform(500, 1000))+'\n')

        MEC_WF.write('---SeNB end---\n\n\n')
    MEC_WF.write('---File end---\n')
    MEC_WF.close()
    TaskDataSize_WF.close()
    TaskCPUCycleNumber_WF.close()
    # 点击窗口关闭
    window = turtle.Screen()
    window.exitonclick()


MEC_radius = 100  # The radius of [10,20] network
SeNB_radius = 50  # The radius of SeNB
SeNB_number = 5
SMD_range = [200,200]
channelNumber = 200
taskNumber_range = [20, 40]
SeNBCoordinateSet = randomlyGenerateMECNetwork(MEC_radius, SeNB_radius, SeNB_number, SMD_range, channelNumber, taskNumber_range)






