# 参考文档
# networkx: https://networkx.github.io/documentation/stable/tutorial.html
# NDlib: https://ndlib.readthedocs.io/en/latest/tutorial.html

import random
import bokeh
import networkx as nx
import matplotlib.pyplot as plt
import ndlib.models.ModelConfig as mc
import ndlib.models.CompositeModel as gc
import ndlib.models.compartments.EdgeNumericalAttribute as ENA
from ndlib.viz.bokeh.DiffusionTrend import DiffusionTrend

N = 10 ** 3  # 网络规模
K = 4  # 平均度
P = K / (N - 1)  # ER连边概率, k = p * (n - 1)
MU = 1  # 恢复概率μ
RHO_0 = 0.15  # 初始感染密度ρ0

# nx.draw(er)
# plt.show()

# ======= 全接触模式 ========
def reactive_process(graph, q, w, work_p=0.15):
    """
       graph: 底层的图结构--ER/WS/BA
       q: 高权重/朋友 比例
       w：高权重边/朋友关系 的权重 >= 1
       work_p: 工作关系感染率
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
    friend_com = ENA('weight', value=w, op='==',
            probability=friend_p)  # 朋友条件
    work_com = ENA('weight', value=1, op='==',
            probability=work_p)  # 工作条件
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

    # 运行模拟
    iterations = model.iteration_bunch(40)

    # 显示
    trends = model.build_trends(iterations)
    viz = DiffusionTrend(model, trends)
    pic = viz.plot(width=800, height=800)
    bokeh.io.show(pic)

    return iterations

def threshold_simula(iterations):
    """求模拟结果的爆发阈值."""

    return 0

def threshold_formula(q, w, k):
    """爆发阈值公式.
       q: 高权重边比例
       w: 高权重边的权重
       k: 平均度
    """
    return k / ((1 - q + q * w) * (k ** 2))


# ======= 单点接触模式 ========

if __name__ == '__main__':
    er = nx.erdos_renyi_graph(N, P)  # ER随机图
    # ws = nx.watts_strogatz_graph(N, K, 0.3)  # WS小世界
    # ba = nx.barabasi_albert_graph(N, K)  # BA无标度网络

    iterations = reactive_process(er, 0.3, 1)
    thr = threshold_simula(iterations)
