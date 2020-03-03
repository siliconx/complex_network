#!/bin/bash

# 因为ndlib不能利用多核cpu并行模拟，故这里为每一对graph-w-q参数建一个python进程
# 共有3*(7*5+9*5)=240个进程

for g in er ws ba
do
    # w为变量
    for w1 in {1..7}  # 个位数
    do
        for w2 in {0..9..2}  # 小数位
        do
            nohup python3 sis.py $g w $w1.$w2 0.1 &
        done
    done

    # q为变量
    for q1 in {0..8}  # 第一位小数
    do
        for q2 in {0..9..2}  # 第二位小数
        do
            nohup python3 sis.py $g q 2 0.$q1$q2 &
        done
    done
done
