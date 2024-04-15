# coding=utf-8

import random
# def generate_dag(n):
#     # 初始化图，edges存储有向边
#     edges = {i: [] for i in range(n)}
#
#     # 添加边，确保图为DAG
#     for i in range(n):
#         for j in range(i + 1, n):
#             # 以一定概率添加边，这里假设为50%
#             if random.random() > 0.5:
#                 edges[i].append(j)
#
#     # 构建输出
#     output = []
#     for i in range(n):
#         # 找先序节点集合
#         preds = [str(k) for k, succs in edges.items() if i in succs]
#         # 后继节点集合已知
#         succs = [str(j) for j in edges[i]]
#         # 格式化输出
#         output.append(f"{','.join(preds)}:{i}:{','.join(succs)}")
#
#     return output


# 示例：生成一个包含5个节点的DAG
# n = 10
# dag_output = generate_dag(n)
# for line in dag_output:
#     print(line)



# import random
#
#
# def generate_connected_dag(n):
#     # 确保至少有一条从起点到终点的路径
#     edges = [(i, i + 1) for i in range(n - 1)]
#
#     # 随机添加其他边，保证无环
#     for i in range(n):
#         for j in range(i + 2, n):  # 从 i+2 开始，避免重复添加已经存在的边和保证无环
#             if random.random() > 0.5:  # 添加边的概率
#                 edges.append((i, j))
#
#     return edges
#
#
# def format_dag_output(n, edges):
#     # 初始化每个节点的先序节点集合和后继节点集合
#     predecessors = {i: [] for i in range(n)}
#     successors = {i: [] for i in range(n)}
#
#     for edge in edges:
#         i, j = edge
#         successors[i].append(j)
#         predecessors[j].append(i)
#
#     # 构建输出格式
#     output = []
#     for i in range(n):
#         # 转换为字符串格式
#         pred_str = ",".join(map(str, sorted(predecessors[i])))
#         succ_str = ",".join(map(str, sorted(successors[i])))
#         output.append(f"{pred_str}:{i}:{succ_str}")
#
#     return output
#
#
# def generate_and_print_dag(n):
#     edges = generate_connected_dag(n)
#     output = format_dag_output(n, edges)
#     for line in output:
#         print(line)

# # 示例
# n = 100
# generate_and_print_dag(n)


# import random
#
#
# def add_edge(predecessors, successors, i, j):
#     """尝试添加边，同时确保每个节点的先序节点数量不超过5"""
#     if len(predecessors[j]) < 5:
#         predecessors[j].append(i)
#         successors[i].append(j)
#         return True
#     return False
#
#
# def generate_dag_with_constraints(n):
#     predecessors = {i: [] for i in range(n)}
#     successors = {i: [] for i in range(n)}
#
#     # 确保从起点到终点的连通性
#     for i in range(n - 1):
#         add_edge(predecessors, successors, i, i + 1)
#
#     # 随机添加其他边，同时确保每个节点的先序节点数不超过5
#     for _ in range( n):  # 尝试次数，可以根据需要调整
#         i, j = random.randint(0, n - 2), random.randint(1, n - 1)
#         if i < j:
#             add_edge(predecessors, successors, i, j)
#
#     return predecessors, successors
#
#
# def format_dag_output(predecessors, successors):
#     n = len(predecessors)
#     output = []
#     for i in range(n):
#         pred_str = ",".join(map(str, sorted(predecessors[i])))
#         succ_str = ",".join(map(str, sorted(successors[i])))
#         output.append(f"{pred_str}:{i}:{succ_str}")
#     return output
#
#
# def generate_and_print_dag(n):
#     predecessors, successors = generate_dag_with_constraints(n)
#     output = format_dag_output(predecessors, successors)
#     for line in output:
#         print(line)
#
#
# # 示例
# n = 100
# generate_and_print_dag(n)


import random


def generate_dag(n):
    # 初始化每个节点的先序节点和后继节点列表
    predecessors = {i: [] for i in range(n)}
    successors = {i: [] for i in range(n)}

    # 确保图的连通性：构建一条从起点到终点的路径
    for i in range(n - 1):
        successors[i].append(i + 1)
        predecessors[i + 1].append(i)

    # 随机添加其他边，同时遵守给定的规则
    for i in range(n):
        # 限制每个节点的后继节点的数量，确保与后继节点的序号相差不超过10
        # 并且保证每个节点的先序节点不超过5个
        possible_successors = range(i + 2, min(i + 11, n))  # 保证序号相差不超过10
        for j in possible_successors:
            if len(successors[i]) < 3 and len(predecessors[j]) < 4:  # 检查先序和后继节点的数量限制
                if random.random() > 0.95:  # 随机决定是否添加边
                    successors[i].append(j)
                    predecessors[j].append(i)

    # 格式化输出
    output = []
    for i in range(n):
        pred_str = ",".join(map(str, sorted(predecessors[i])))
        succ_str = ",".join(map(str, sorted(successors[i])))
        output.append(f"{pred_str}:{i}:{succ_str}")

    return output


def print_dag(n):
    dag_output = generate_dag(n)
    for line in dag_output:
        print(line)

def save_dag(n):
    dag_output = generate_dag(n)
    filename = './SeNB-5//t' + str(n) + '.txt'
    MEC_WF = open(filename, 'w')
    for line in dag_output:
        # print(line)
        MEC_WF.write(line + '\n')

# 示例
for n in range(10,301):
    save_dag(n)

n = 101  # 节点数

