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
K = 4  # 平均度
P = K / (N - 1)  # ER连边概率, k = p * (n - 1)
MU = 1  # 恢复概率μ
RHO_0 = 0.15  # 初始感染密度ρ0
ROUND = 200  # 模拟轮数
ZERO = 0.01 # 用一个非常小的数来表示0

# 可视化网络
# nx.draw(er)
# plt.show()

# ======= 全接触模式 ========
def reactive_process(graph, w, q, work_p=0.001, show_trends=False):
    """
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
    recover_com = ENA('weight', value=0, op='!=', probability=MU,
            triggering_status='Susceptible')  # 恢复条件，用于任意边

    # 添加规则到模型
    model.add_rule('Susceptible', 'Infected', friend_com)
    model.add_rule('Susceptible', 'Infected', work_com)
    model.add_rule('Infected', 'Susceptible', recover_com)

    # 模型初始配置
    cfg = mc.Configuration()
    cfg.add_model_parameter('fraction_infected', RHO_0)
    model.set_initial_status(cfg)

    # 运行模拟
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
    在实际模拟中，平均感染密度并不为0，因为网络中有的节点是孤立的
    """
    work_p = 0  # 感染率λ
    mean_r = 0  # 平均感染密度
    while mean_r <= ZERO:
        infect_rates = []
        iterations = reactive_process(graph, w, q, work_p)  # 用当前参数进行模拟

        for i in iterations:
            infect_rates.append(i['node_count'][1] / N)
        print(infect_rates)
        # 用最小感染密度作为平均感染密度
        mean_r = min(infect_rates)
        work_p += 0.002
    thr = work_p / MU
    return thr

def threshold_formula(w, q):
    """爆发阈值公式.
       w: 高权重边的权重
       q: 高权重边比例
       k: 平均度
    """
    return K / ((1 - q + q * w) * (K ** 2))


# ======= 单点接触模式 ========


if __name__ == '__main__':
    er = nx.erdos_renyi_graph(N, P)  # ER随机图
    # ws = nx.watts_strogatz_graph(N, K, 0.3)  # WS小世界
    # ba = nx.barabasi_albert_graph(N, K)  # BA无标度网络

    # iterations = reactive_process(er, 0.3, 1)
    thr_simu = threshold_simula(er, 2, 0.1)
    thr_form = threshold_formula(2, 0.1)
    print('simulation:', thr_simu)
    print('formulation:', thr_form)
