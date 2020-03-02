import os
import pandas as pd
import matplotlib.pyplot as plt

def threshold_plot():
    """绘制爆发阈值."""
    graphs = ['er', 'ws', 'ba']  # 网络图
    line_style = [('ro:', 'r-'), ('b^--', 'b-'), ('gs-.', 'g-')]  # 对应的样式
    fig, axes = plt.subplots(2, 1)
    for ax in axes:
            ax.set_ylabel('λ_c')

    for g, ls in zip(graphs, line_style):
        file = g + '.csv'
        if not os.path.exists(file):  # 忽略不存在的文件
            continue

        df = pd.read_csv(file, sep='\t')
        for ax, v in zip(axes, ['w', 'q']):  # 分别对w和q画图
            var = df[df['variable'] == v]  # 取出w/q的行
            new_var = var[['simula', 'formula']]
            new_var.index = var[v]
            new_var = new_var.sort_index()  # 按w/q排序
            for idx, col in enumerate(new_var.columns):
                new_var[col].plot(ax=ax, style=ls[idx], label='%s_%s' % (g.upper(), col))
                if v == 'w':
                    ax.set_xlabel('w (q=0.1)')
                else:
                    ax.set_xlabel('q (w=2)')
                ax.legend()
    plt.show()

if __name__ == '__main__':
    threshold_plot()
