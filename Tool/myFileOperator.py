# coding=gbk
import os

def getCurrentPath():
    return os.path.dirname(os.path.realpath(__file__))


def getProjectPath():
    cur_path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(os.path.dirname(cur_path))


def copy_rename(old_name,new_name):

    # 3.�����ļ�д������(���ݺ�ԭ�ļ�һ��)
    # 3.1 ��ԭ�ļ��ͱ����ļ�
    old_f = open(old_name, 'rb')  # �Զ����ƶ��ķ�ʽ����ԭ�ļ�
    new_f = open(new_name, 'wb')  # �Զ�����д�ķ�ʽ����ԭ�ļ�

    # 3.2 ԭ�ļ����룬�����ļ�д��
    # ע�⣺�����ȷ��Ŀ���ļ��Ĵ�С��ѭ����ȡд�룬����ȡ����������û���˲���ֹ
    while True:
        con = old_f.read(1024)  # һ�ζ���1024���ֽ�
        if len(con) == 0:
            break
        new_f.write(con)



# path = 'D:/ÿ�����/����/���׺�����/MCOP�ο�����/COP-Resource-Allocation -4-26 result/My_Makespan_Energy/VerifyEffectiveness/ExperimentResult/'
# taskNumberRangeList = ['[20,30]', '[30,40]', '[40,50]', '[50,60]', '[60,70]', '[70,80]']
# fileType = '.xls'
#
#
# for taskNumberRange in taskNumberRangeList:
#     old_name = path + taskNumberRange + '/' + 'MOWOA' + fileType
#     new_name = path + taskNumberRange + '/' + 'MOPSO' + fileType
#     copy_rename(old_name,new_name)