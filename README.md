# tudatspace_example
行星际脉冲转移轨道算例

#tudat-space

1.tudat-space在window平台上install的时候老是包冲突，好在它官网上提供了一个在线binder环境，里边都配置好了，直接进。

2.mga DSM at  pe  fixtime 这个算例是固定leg之间飞行时间的，它的DSM是发生在飞跃的近点的，所以又叫MPGA模型。

3.mga noDSM optimize use DE 这个算例不固定leg之间飞行时间，使用pygmo的差分进化算法进行迭代，求解最优，但是它只有出发时的一次点火，中途不进行深空机动。

4.我的pykep_example里有基于pykep和pygmo的mga-1dsm模型算例，我们仍然需要基于tudat-space和pygmo的mga-1dsm模型算例，尚未发现有人给出。
