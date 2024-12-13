## tudatspace_example

行星际脉冲转移轨道算例

## tudat-space

1.tudat-space在window平台上install的时候老是包冲突，好在它官网上提供了一个在线binder环境，里边都配置好了，直接进。

2.mga DSM at  pe  fixtime 这个算例是固定leg之间飞行时间的，它的DSM是发生在飞跃之间的。

3.mga noDSM optimize use DE 这个算例不固定leg之间飞行时间，使用pygmo的差分进化算法进行迭代，求解最优，但是它只有出发时的一次点火，中途不进行深空机动。

4.我的pykep_example(https://github.com/Quantum-dogdog/pykep_example)里有基于pykep和pygmo的mga-1dsm模型优化算例，我们仍然需要基于tudat-space和pygmo的mga-1dsm模型优化算例，尚未发现有人给出。

## 2024-12-13更新

1.今天我通过 mga1dsm optimize use pygmo（sade）这个算例 自己回答了上边的第4个问题，在过去的三天里，我在智谱清言(https://chatglm.cn/)的帮助下，成功在全网率先给出开源版的tudat-space + pygmo版本的mga1dsm模型优化算例，与我的pykep_example(https://github.com/Quantum-dogdog/pykep_example)中提到的一样，也是算的EVEEJ算例。

2.这标志着我已经具备了自己解决问题的能力，是我学习轨道力学的重要里程碑！
