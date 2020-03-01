#!/bin/bash

# 因为ndlib不能利用多核cpu并行模拟，故这里为每一对graph-w-q建一个python进程
for g in er ws ba
do
    for w in {1..7}
    do
        nohup python3 sis.py $g $w 0.1 &
    done

    for q in {0..8}
    do
        val=`echo "scale=2;$q/10"|bc`
        nohup python3 sis.py $g 2 $val &
    done
done