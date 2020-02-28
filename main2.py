# 参考文档
# networkx: https://networkx.github.io/documentation/stable/tutorial.html
# NDlib: https://ndlib.readthedocs.io/en/latest/tutorial.html

import random
import bokeh
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import ndlib.models.ModelConfig as mc
import ndlib.models.CompositeModel as gc
import ndlib.models.compartments.EdgeNumericalAttribute as ENA
from ndlib.viz.bokeh.DiffusionTrend import DiffusionTrend

N = 500  # 网络规模
K = 8  # 平均度
P = K / (N - 1)  # ER连边概率, k = p * (n - 1)
MU = 1  # 恢复概率μ
RHO_0 = 0.15  # 初始感染密度ρ0
ROUND = 250  # 模拟轮数，时间步
STEP = 0.001 # 感染率步长

# 可视化网络
# nx.draw(er)
# plt.show()

# ======= 全接触模式 ========
def simulation(graph, w, q, work_p=0.001, show_trends=False):
    """传染病模拟.
       graph: 底层的图结构--ER/WS/BA
       w：高权重边/朋友关系 的权重 >= 1
       q: 高权重/朋友 比例
       work_p: 工作关系感染率,即感染率λ
    """
    friend_p = 1 - (1 - work_p) ** w  # 朋友关系感染概率

    attr = {}  # 边属性
    edges = graph.edges()  # 所有的边
    len_edges = len(edges)  # 边的数量
    friend_e = random.sample(edges, int(q*len_edges))  # 朋友关系
    work_e = set(edges) - set(friend_e)  # 工作关系

    # 设置边属性
    for i in friend_e:
        attr[i] = {
            'weight': w,
        }
    for i in work_e:
        attr[i] = {
            'weight': 1,
        }
    nx.set_edge_attributes(graph, attr)

    # 组合模型
    model = gc.CompositeModel(graph)

    # 模型状态
    model.add_status('Susceptible')
    model.add_status('Infected')

    # compartment
    friend_com = ENA('weight', value=w, op='==', probability=friend_p,
            triggering_status='Infected')  # 朋友条件
    work_com = ENA('weight', value=1, op='==', probability=work_p,
            triggering_status='Infected')  # 工作条件
    recover_com = ENA('weight', value=0, op='!=',
            probability=MU)  # 恢复条件，用于任意边

    # 添加规则到模型
    model.add_rule('Susceptible', 'Infected', friend_com)
    model.add_rule('Susceptible', 'Infected', work_com)
    model.add_rule('Infected', 'Susceptible', recover_com)

    # 模型初始配置
    cfg = mc.Configuration()
    cfg.add_model_parameter('fraction_infected', RHO_0)
    model.set_initial_status(cfg)

    # 运行模拟，得到模拟的每一步结果
    iterations = model.iteration_bunch(ROUND)

    # 显示变化趋势
    if show_trends:
        trends = model.build_trends(iterations)
        viz = DiffusionTrend(model, trends)
        pic = viz.plot(width=800, height=800)
        bokeh.io.show(pic)

    return iterations

def threshold_simula(graph, w, q):
    """求模拟结果的爆发阈值.
    论文引用：在模拟过程中, 当w和q确定时, 随着感染率λ的增加最终稳定的平均感染密度ρ(t > 2500)
    将从0变为非0, 从而可以获得爆发阈值λ_c.
    graph: 底层的图结构--ER/WS/BA
    w：高权重边/朋友关系 的权重 >= 1
    q: 高权重/朋友 比例
    """
    work_p = 0  # 感染率λ
    while True:
        iterations = simulation(graph, w, q, work_p)  # 用当前参数进行模拟
        density = infected_density(graph, iterations)  # 平均感染密度
        if density > 0:
            break
        work_p += STEP

    thr = work_p / MU
    return thr

def threshold_formula(w, q):
    """爆发阈值公式.
       w: 高权重边的权重
       q: 高权重边比例
    """
    return K / ((1 - q + q * w) * (K ** 2))

def infected_density(graph, iterations):
    """根据模拟结果计算平均感染密度(忽略度为0的感染节点).
    graph: 底层的图结构--ER/WS/BA
    iterations: 模拟的每一步结果
    """
    infected_n = N + 1  # 感染节点数量
    index = 0  # 最小感染密度的位置
    zero_count = 0  # 记录度数为0的感染节点数
    zero_degree = set()  # 度为0的节点

    for i in graph.degree:  # 找到网络中所有度数为0的节点
        if i[1] == 0:
            zero_degree.add(i)

    # 取迭代结果中的最小平均感染密度
    print('infected_nodes: [', end='')
    for idx, itr in enumerate(iterations):
        temp = itr['node_count'][1]
        print(temp, end=', ')
        if temp < infected_n:
            infected_n = temp
            index = idx
    print(']')
    status = iterations[index]['status']

    # 统计度数为0的感染节点
    for s in status:
        if status[s] == 1 and s in zero_degree:
            zero_count += 1

    return (infected_n - zero_count) / N



# ======= 单点接触模式 ========


if __name__ == '__main__':
    er = nx.erdos_renyi_graph(N, P)  # ER随机图
    # ws = nx.watts_strogatz_graph(N, K, 0.3)  # WS小世界
    # ba = nx.barabasi_albert_graph(N, K)  # BA无标度网络

    # for w in range(1, 8):
    #     thr_simu = threshold_simula(er, w, 0.1)
    #     thr_form = threshold_formula(w, 0.1)
    #     print('w = %d: simulation = %f, formulation = %f' % (w, thr_simu, thr_form))

    for i in range(1, 9):
        q = i / 10
        thr_simu = threshold_simula(er, 2, q)
        thr_form = threshold_formula(2, q)
        print('q = %f: simulation = %f, formulation = %f' % (q, thr_simu, thr_form))
