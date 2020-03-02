import os
import sys
import csv
import time
import random
import bokeh
import logging
import networkx as nx
import matplotlib.pyplot as plt
import ndlib.models.epidemics as ep
import ndlib.models.ModelConfig as mc
from ndlib.viz.bokeh.DiffusionTrend import DiffusionTrend

N = 10 ** 4  # 网络规模(论文为10**4)
K = 8  # 平均度
MU = 1  # 恢复概率μ
RHO_0 = 0.15  # 初始感染密度ρ0
TIMES = 2500  # 模拟轮数，时间步(论文为2500)
INIT_WORK_P = 0.1  # 感染率初值
INIT_STEP = INIT_WORK_P / 4  # 感染率变化的初始步长
PRECISION = 0.0001  # step_v的精度


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

        if graph_name == 'er':  # ER随机图
            edge_p = K / (N - 1)  # ER连边概率, k = p * (n - 1)
            self.graph = nx.erdos_renyi_graph(N, edge_p)
        elif graph_name == 'ws':  # WS小世界
            self.graph = nx.watts_strogatz_graph(N, K, 0.3)
        elif graph_name == 'ba':  # BA无标度网络
            # 经实验，m = <k> / 2
            self.graph = nx.barabasi_albert_graph(N, K / 2)
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
        # 对于网络中任意一条边, 利用平均场近似, 其在 τ = 1 时间步内成功传
        # 播病毒的概率。通过这个公式，把双关系网络转化为单一关系网络
        p = (1 - self.q) * work_p + \
                self.q * (1 - (1 - work_p) ** w)

        # 选择SIS模型
        model = ep.SISModel(self.graph)

        # 模型配置
        cfg = mc.Configuration()
        cfg.add_model_parameter('beta', p)
        cfg.add_model_parameter('lambda', MU)
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
        论文引用：在模拟过程中, 当w和q确定时, 随着感染率λ的增加最终稳定的
        平均感染密度ρ(t > 2500)将从0变为非0, 从而可以获得爆发阈值λ_c.
        """
        work_p = INIT_WORK_P  # 初始感染率λ，给初值，减少迭代次数
        step_v = INIT_STEP  # 变化的步长
        if self.graph_name == 'ba':  # ba无标度网络的初值要更小
            work_p /= 2
            step_v /= 2

        operation = None  # 控制步长变化的开关
        cache = dict()  # 记录已经计算过的值，减少计算量(震荡求解的时候某些值会重复)
        while True:
            start = time.time()
            den_m = cache.get(work_p)
            if den_m:  # 已经计算过此work_p
                density = den_m
            else:  # 新的work_p
                iterations = self.simulation(work_p)  # 用当前参数进行模拟
                density = self.infected_density(iterations)  # 平均感染密度
                cache[work_p] = density  # 记录
            end = time.time()
            time_used = (end - start) / 60  # mins
            self.logger.info(
                ('w = %d, q = %.5f, density = %.5f, work_p = %.5f, ' + \
                'step_v = %.5f, variable = %s, %.2f mins used') % \
                (self.w, self.q, density, work_p, step_v, self.variable,
                    time_used))

            # 震荡求解，用越来越小的步长逐步逼近实际的感染率，类似二分搜索
            if density > 0:  # 感染密度大于0，需降低感染率
                if operation == None:
                    operation = '-'
                elif operation == '+':  # 方向变化，步长减半
                    operation = '-'
                    step_v /= 2
            else:  # 感染密度等于0，需增大感染率
                if operation == None:
                    operation = '+'
                elif operation == '-':  # 方向变化，步长减半
                    operation = '+'
                    step_v /= 2

            # 先判断精度再修改work_p
            if step_v < PRECISION and density <= 0:
                break

            if operation == '+':
                work_p += step_v  # 增加感染率
            elif operation == '-':
                work_p -= step_v  # 减小感染率

        return work_p / MU  # 爆发阈值 = 感染率 / 恢复率

    def infected_density(self, iterations):
        """根据模拟结果计算平均感染密度.
        iterations: 模拟的每一步结果
        """
        infected_n = N + 1  # 感染节点数量

        # 取感染节点数最少的一个结果
        for idx, itr in enumerate(iterations):
            temp = itr['node_count'][1]
            if temp < infected_n:
                infected_n = temp

        return infected_n / N

    def threshold_formula(self):
        """爆发阈值公式.
        """
        d = 0  # 度之和
        d_2 = 0  # 度的平方之和
        for i in self.graph.degree:
            d += i[1]
            d_2 += i[1] ** 2

        k = d / N  # 平均度, <k>
        k_2 = d_2 / N  # 平均平方度, <k^2>

        return k / ((1 - self.q + self.q * self.w) * k_2)

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
        self.save2file((self.w, self.q, '%.5f' % thr_simu,
            '%.5f' % thr_form, self.variable))
        self.logger.info(('N = %d, TIMES = %d, w = %d, q = %.5f done' +\
            ', %.2f mins used') % (N, TIMES, self.w, self. q, time_used))


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
