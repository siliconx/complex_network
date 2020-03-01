#!/bin/bash

# 因为ndlib不能利用多核cpu并行模拟，故这里为每一对graph-w-q参数建一个python进程
# 共有3*(7+9)=48个进程
for g in er ws ba
do
    for w in {1..7}
    do
        nohup python3 sis.py $g w $w 1 &
    done

    for q in {0..8}  # bash不支持小数，在python中转换
    do
        nohup python3 sis.py $g q 2 $q &
    done
done
