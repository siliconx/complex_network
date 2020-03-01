import os
import sys
import csv
import time
import random
import bokeh
import logging
import networkx as nx
import matplotlib.pyplot as plt
import ndlib.models.ModelConfig as mc
import ndlib.models.CompositeModel as gc
import ndlib.models.compartments.EdgeNumericalAttribute as ENA
from ndlib.viz.bokeh.DiffusionTrend import DiffusionTrend

N = 10 ** 3  # 网络规模(论文为10**4)
K = 8  # 平均度
P = K / (N - 1)  # ER连边概率, k = p * (n - 1)
MU = 1  # 恢复概率μ
RHO_0 = 0.15  # 初始感染密度ρ0
TIMES = 200  # 模拟轮数，时间步(论文为2500)
STEP = 0.1 # 感染率变化的初始步长
PRECISION = 0.0001  # 步长的变化精度


class ReactiveProcess(object):
    """全接触模式."""

    def __init__(self, graph_name, variable, w, q):
        """初始化.
           graph_name: 底层的图结构--ER/WS/BA
           variable: 变化的量(w/q). 如果为w，则q固定
           w：高权重边/朋友关系 的权重 >= 1
           q: 高权重/朋友 比例
        """
        super(ReactiveProcess, self).__init__()

        if graph_name == 'er':
            self.graph = nx.erdos_renyi_graph(N, P)  # ER随机图
        elif graph_name == 'ws':
            self.graph = nx.watts_strogatz_graph(N, K, 0.3)  # WS小世界
        elif graph_name == 'ba':
            self.graph = nx.barabasi_albert_graph(N, K)  # BA无标度网络
        else:
            raise ValueError('graph name = er/ws/ba')
        self.graph_name = graph_name

        if variable not in ['w', 'q']:
            raise ValueError('variable = w/s')
        self.variable = variable

        if not (1 <= w <= 10):
            raise ValueError('w = [1, 10]')
        self.w = w

        if not (0 <= q <= 1):
            raise ValueError('q = [0, 1]')
        self.q = q

        # 配置日志
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', 
            filename='%s.log' % self.graph_name)
        self.logger = logging.getLogger(__name__)

    def simulation(self, work_p, show_trends=False):
        """传染病模拟.
           work_p: 工作关系感染率,即感染率λ
        """
        friend_p = 1 - (1 - work_p) ** self.w  # 朋友关系感染概率

        attr = {}  # 边属性
        edges = self.graph.edges()  # 所有的边
        len_edges = len(edges)  # 边的数量
        friend_e = random.sample(edges, int(self.q*len_edges))  # 朋友关系
        work_e = set(edges) - set(friend_e)  # 工作关系

        # 设置边属性
        for i in friend_e:
            attr[i] = {
                'weight': 2,  # 边的标记，不是真正的权重，不用来计算
            }
        for i in work_e:
            attr[i] = {
                'weight': 1,  # 边的标记，不是真正的权重，不用来计算
            }
        nx.set_edge_attributes(self.graph, attr)

        # 组合模型
        model = gc.CompositeModel(self.graph)

        # 模型状态
        model.add_status('Susceptible')
        model.add_status('Infected')

        # compartment
        friend_com = ENA('weight', value=2, op='==', probability=friend_p,
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
        iterations = model.iteration_bunch(TIMES)

        # 显示变化趋势
        if show_trends:
            trends = model.build_trends(iterations)
            viz = DiffusionTrend(model, trends)
            pic = viz.plot(width=800, height=800)
            bokeh.io.show(pic)

        return iterations

    def threshold_simula(self):
        """求模拟结果的爆发阈值.
        论文引用：在模拟过程中, 当w和q确定时, 随着感染率λ的增加最终稳定的平均感染密度ρ(t > 2500)
        将从0变为非0, 从而可以获得爆发阈值λ_c.
        """
        work_p = STEP  # 初始感染率λ，给初值，减少迭代次数
        step_v = STEP / 2  # 变化的步长
        operation = None  # 控制步长变化的开关
        memory = dict()  # 记录已经计算过的值，减少计算量(震荡求解的时候某些值会重复)
        while True:
            den_m = memory.get(work_p)
            if den_m:  # 已经计算过此work_p
                # print('(old work_p)', end='')
                density = den_m
            else:  # 新的work_p
                iterations = self.simulation(work_p)  # 用当前参数进行模拟
                density = self.infected_density(iterations)  # 平均感染密度
                memory[work_p] = density  # 记录
            self.logger.info('w = %d, q = %.5f, density = %.5f, work_p = %.5f, step_v = %.5f' %\
                (self.w, self.q, density, work_p, step_v))

            # 震荡求解，用越来越小的步长逐步逼近实际的感染率，类似二分搜索
            if density > 0:
                if operation == None:
                    operation = '-'
                elif operation == '+':  # 方向变化，步长减半
                    operation = '-'
                    step_v /= 2
            elif density == 0:
                if operation == None:
                    operation = '+'
                elif operation == '-':  # 方向变化，步长减半
                    operation = '+'
                    step_v /= 2

            if step_v < PRECISION:  # 先判断精度再修改work_p
                break

            if operation == '+':
                work_p += step_v  # 增加感染率
            elif operation == '-':
                work_p -= step_v  # 减小感染率


        return work_p / MU  # 爆发阈值 = 感染率 / 恢复率

    def infected_density(self, iterations):
        """根据模拟结果计算平均感染密度(忽略度为0的感染节点).
        iterations: 模拟的每一步结果
        """
        infected_n = N + 1  # 感染节点数量
        index = 0  # 最小感染密度的位置
        zero_count = 0  # 记录度数为0的感染节点数
        zero_degree = set()  # 度为0的节点

        for i in self.graph.degree:  # 找到网络中所有度数为0的节点
            if i[1] == 0:
                zero_degree.add(i)

        # 取迭代结果中的最小平均感染密度
        for idx, itr in enumerate(iterations):
            temp = itr['node_count'][1]
            if temp < infected_n:
                infected_n = temp
                index = idx
        status = iterations[index]['status']

        # 统计度数为0的感染节点
        for s in status:
            if status[s] == 1 and s in zero_degree:
                zero_count += 1

        return (infected_n - zero_count) / N

    def threshold_formula(self):
        """爆发阈值公式.
        """
        return K / ((1 - self.q + self.q * self.w) * (K ** 2))

    def save2file(self, row):
        """保存结果."""
        file_name = '%s.csv' % (self.graph_name)
        exist = False
        if os.path.exists(file_name):  # 检测文件是否存在
            exist = True

        with open(file_name, 'a') as f:
            csv_wrt = csv.writer(f, delimiter='\t')
            if not exist:  # 新文件
                csv_wrt.writerow(['w', 'q', 'simula', 'formula', 'variable'])
            csv_wrt.writerow(row)

    def show(self):
        """显示网络拓扑结构."""
        nx.draw(self.graph)
        plt.show()

    def run(self):
        """运行程序."""
        start = time.time()
        thr_simu = self.threshold_simula()
        thr_form = self.threshold_formula()
        end = time.time()
        time_used = (end - start) / 60  # mins
        self.save2file((self.w, self.q, '%.5f' % thr_simu, '%.5f' % thr_form, self.variable))
        self.logger.info('w = %d, q = %f done, %.2f mins used' %\
            (self.w, self. q, time_used))


class ContactProcess(object):
    """单点接触模式."""
    pass


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('4 args are needed!:\n(1)graph name(er/ws/ba)\n' +\
            '(2)variable(w/q)\n(3)w[1-10]\n(4)q[0-10]')
        exit()
    graph_name = sys.argv[1]
    variable = sys.argv[2]
    w = int(sys.argv[3])
    q = float(sys.argv[4]) / 10  # bash不支持小数，故在python中转换
    rp = ReactiveProcess(graph_name, variable, w, q)
    rp.run()
