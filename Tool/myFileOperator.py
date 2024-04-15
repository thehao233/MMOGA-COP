# coding=gbk
import os

def getCurrentPath():
    return os.path.dirname(os.path.realpath(__file__))


def getProjectPath():
    cur_path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(os.path.dirname(cur_path))


def copy_rename(old_name,new_name):

    # 3.备份文件写入数据(数据和原文件一样)
    # 3.1 打开原文件和备份文件
    old_f = open(old_name, 'rb')  # 以二进制读的方式，打开原文件
    new_f = open(new_name, 'wb')  # 以二进制写的方式，打开原文件

    # 3.2 原文件读入，备份文件写入
    # 注意：如果不确定目标文件的大小，循环读取写入，当读取出来的数据没有了才终止
    while True:
        con = old_f.read(1024)  # 一次读入1024个字节
        if len(con) == 0:
            break
        new_f.write(con)



# path = 'D:/每周组会/资料/文献含代码/MCOP参考代码/COP-Resource-Allocation -4-26 result/My_Makespan_Energy/VerifyEffectiveness/ExperimentResult/'
# taskNumberRangeList = ['[20,30]', '[30,40]', '[40,50]', '[50,60]', '[60,70]', '[70,80]']
# fileType = '.xls'
#
#
# for taskNumberRange in taskNumberRangeList:
#     old_name = path + taskNumberRange + '/' + 'MOWOA' + fileType
#     new_name = path + taskNumberRange + '/' + 'MOPSO' + fileType
#     copy_rename(old_name,new_name)