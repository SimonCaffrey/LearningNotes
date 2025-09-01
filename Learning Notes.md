## 学习笔记

### 无线通信

#### · 《无线通信基础》学习笔记

***关于多天线技术***：我们知道，如果发射机已知信道，那么采用多副天线通过发射波束成形可以提供功率增益；同时多副发射天线还可用于产生信道波动，满足机会通信技术需求（解释为机会波束成形）。多天线技术采用多副发射天线和多副接收天线，提供了用于通信的额外的空间维数，并产生自由度增益，利用这些额外的自由度可以将若干数据流空间多路复用到MIMO信道中。多天线技术包含了MIMO、MISO、SIMO三种，能通过信号处理技术获得信噪比或传输速率的改善，具体改善可以通过**阵列增益**、**分集增益**、**复用增益**，三个方面来描述。

1. 阵列增益：通过天线阵列提高平均信噪比
2. 分集增益：增加传播路径数目，改善了接受信噪比的概率分布，改善误码率
3. 复用增益：也称“自由度”，改善可以同时传输的独立数据流数目

##### *第七章*

**章序**：研究实现空间多路复用的物理环境的属性，阐明如何获得这些属性。先通过容量分析确定MIMO信道多路复用容量的关键参数

###### 7.1 确定性MIMO信道的多路复用容量

​	假设有$n_t$副发射天线和$n_r$副接收天线，那么窄带时不变无线信道可以表示为一个$n_r \times n_t$的矩阵$\mathbf{H}$，从而可以得到时不变信道为：
$$
\mathbf{y}=\mathbf{Hx}+\mathbf{\omega}
$$
​	其中$\mathbf{x}\in\mathcal{C}^{n_t},\mathbf{y}\in\mathcal{C}^{n_r},\mathbf{\omega}$分别表示一个码元时刻的发射信号、接收信号和高斯白噪声，$\mathbf{H}$为确定性的信道模型，考虑到信道系统为时不变系统，该模型在所有时刻都保持不变。其中把$\mathbf{H}$做SVD分解可以有：
$$
\mathbf{H}=\mathbf{U} \Lambda \mathbf{V}^\mathrm{H}
$$
​	其中$\mathbf{U}\in\mathcal{C}^{n_r\times n_r},\mathbf{V}\in\mathcal{C}^{n_t\times n_t}$皆为酉矩阵（满足$\mathbf{A}^\mathrm{H}=\mathbf{A}^{-1} $or $ \mathbf{A}^\mathrm{H}\mathbf{A}=\mathbf{I}$的矩阵），$\mathbf{\Lambda}\in\mathcal{R}^{n_r\times n_t}$​是对角矩阵，且对角线元素非负实数，非对角线元素为0。对角线元素$\lambda_1\geq\lambda_2\geq\cdots\geq\lambda_{n_\min}$为$\mathbf{H}$的奇异值，其中$n_\min=\min\pqty{n_t,n_r}$。同时，由于$\mathbf{HH}^\mathrm{H}=\mathbf{U\Lambda V}^\mathrm{H}\mathbf{V\Lambda}^\mathrm{H}\mathbf{U}^\mathrm{H}=\mathbf{U\Lambda\Lambda}^\mathrm{H}\mathbf{U}^\mathrm{H}$（$\Lambda$是实数，所以$\Lambda^\mathrm{H}=\Lambda^\mathrm{T}$，采用$\Lambda^\mathrm{H}$更具有一致性），因此$\lambda^2$为$\mathbf{HH}^\mathrm{H}$和$\mathbf{H}^\mathrm{H}\mathbf{H}$的特征值。SVD的重新表示可以为
$$
\mathbf{H}=\sum_{i=1}^{n_\min}\lambda_i\mathbf{u}_i\mathbf{v}_i^\mathrm{H}
$$
即秩为1的矩阵$\lambda_i\mathbf{u}_i\mathbf{v}_i^\mathrm{H}$之和（每个$\lambda_i\mathbf{u}_i\mathbf{v}_i^\mathrm{H}$都是$\rank$-1，但是相加之后$\rank$就为$n_\min$，其中$\lambda_i$是一个标量从对角矩阵中取出，$\mathbf{u}_i$表示$\mathbf{U}$的一个列向量，$\mathbf{v}_i^*$表示$(\mathbf{V}^*)$的一个行向量）

​	将$\mathbf{H}=\mathbf{U} \Lambda \mathbf{V}^*$代入$\mathbf{y}=\mathbf{Hx}+\mathbf{\omega}$得到$\mathbf{y}=\mathbf{U} \Lambda \mathbf{V}^*\mathbf{x}+\mathbf{\omega}$，若定义
$$
\begin{cases}
\hat{\mathbf{y}}=U^*\times \mathbf{y}\\
\hat{\mathbf{x}}=V^*\times \mathbf{x}\\
\hat{\mathbf{\omega}}=U^*\times \mathbf{\omega}
\end{cases}
$$

​	则可以有$\hat{y}=\Lambda \hat{x} +\hat{w}$，其中，$\norm{\hat{\mathbf{x}}}^2=\norm{\mathbf{x}}^2$，且$\hat{\mathbf{\omega}}$与$\mathbf{\omega}$呈相同的分布。（考虑这里的$\Lambda$是一个对角矩阵，可以认为$\hat{y}$可以作为$\hat{x}$一一对应的信道，就是将原有的$y、x$经过处理后，即使在MIMO中也可以接收发送一一对应）从之类也可以看出$\text{DoF}_\max=\min(n_t,n_r)$。

​	等效的并行高斯信道（标量）形式可以表示为
$$
\hat{y}_i=\lambda_i\hat{x}_i+\hat{\omega}_i
$$
​	因此可以通过上述分析得到，$\hat{x}_i$经过平行的信道$\lambda_i$传输后得到的数据为$\hat{y}_i$。即为了满足传输的信号能近似平行信道的方式进行传输，则应该考虑到SVD的分解情况进行预处理
$$
\begin{align}
\mathbf{y}&=\mathbf{U\Lambda V}^\mathrm{H}\mathbf{x}+\mathbf{\omega}\\
\pqty{\mathbf{U}\hat{\mathbf{y}}}&=\mathbf{U\Lambda V}^\mathrm{H}\pqty{\mathbf{V\hat{x}}}+\pqty{\mathbf{U\hat{\omega}}}\\
\hat{\mathbf{y}}&=\mathbf{\Lambda}\hat{\mathbf{x}}+\hat{\mathbf{\omega}}
\end{align}
$$
可以容易发现$\hat{\mathbf{x}}$就可以沿等效平行信道传输，最终得到$\hat{\mathbf{y}}$。（为便于理解，可以尝试从上到下和从下到上两种方式进行理解和推导）

​	这里可以相似地看待这里空间维起的作用和在正交频率OFDM多频域通信下，频率维起的作用，这两种情况都是利用变换将矩阵信道转换为一组并行的独立子信道，那么实际的信道容量是多个并行信道通信容量的总和
$$
C=\sum_{i=1}^{n_\min}\log\pqty{1+\frac{P_i^*\lambda_i^2}{N_0}}
$$
其中，$P_1^*,\cdots,P_{n_\min}^*$是注水功率分配，表示为$P_i^*=\pqty{\mu-\frac{N_0}{\lambda_i^2}}^+$，$\mu$使其满足总功率约束$\sum_{i=1}^{n_\min}P_i^*=P$。而$\lambda_i$就对应于信道的一个特征模式/特征信道，每个非零特征信道都能支撑一路独立的数据流，因此再次说明验证了$\text{DoF}=\min\pqty{n_t,n_r}$。

###### · 注水功率算法

**目的**：为了更合理分配MIMO系统中各独立信道的发射功率，使其达到信道容量上限，并且满足总发射功率约束。

已知上述信道容量表达形式，则用拉格朗日乘数法有
$$
J\pqty{P_1^*,\cdots,P_n^*}=\sum_{i=1}^{n_\min}\log\pqty{1+\frac{P_i^*\lambda_i^2}{N_0}}+\beta\sum_{i=1}^{n_\min}P_i^*
$$
为求上式的最大值，则对分别其求$P_i^*$的偏导，得到
$$
\frac{\partial J}{\partial P_i^*}=\frac{\lambda_i^2}{N_0+P_i^*\lambda_i^2}+\beta=0
$$
求解上式得到
$$
P_i^{*opt}=-\frac{1}{\beta}-\frac{N_0}{\lambda_i^2}=\mu-\frac{N_0}{\lambda_i^2}
$$
由于发射功率之和受到总功率的约束$\sum_{i=1}^{n_\min}P_i=P$，则可以通过这一约束求得$\mu$
$$
\sum_{i=1}^{n_\min}\pqty{\mu-\frac{N_0}{\lambda_i^2}}=n_{\min}\mu-N_0\sum_{i=1}^{n_\min}\frac{1}{\lambda_i^2}=P
$$
得到$\mu=\frac{P+N_0\sum_{i=1}^{n_\min}\frac{1}{\lambda^2_i}}{n_\min}$，从而可以方便地有
$$
\begin{align}
P_i^{*opt}
&=\frac{P+N_0\sum_{k=1}^{n_\min}\frac{1}{\lambda_k^2}}{n_\min}-\frac{N_0}{\lambda_i^2}\\
&=\frac{P+N_0\sum_{k=1}^{n_\min}\pqty{\frac{1}{\lambda_k^2}-\frac{1}{\lambda_i^2}}}{n_\min}
\end{align}
$$
​	上式就是最终的表达形式，如果一个特征信道的$\lambda_i$较小，则$\frac{1}{\lambda_i^2}$较小，得到的$P_i^{*opt}$就小，这样就不将功率浪费到信道环境较差的特征信道；相反，如果$\lambda_i$较大，则说明信道质量较好，最终得到的$P_i^{*opt}$就更大。

​	需要注意的是，由于上式可能得到小于0的结果，所以需要将其直接置0，然后重新迭代计算充分进行功率分配，直至每一路特征信道的分配功率大于0。

​	接下来，将分**高信噪比**和**低信噪比**两种情况讨论决定MIMO系统性能的关键参数。

​	· **高信噪比时**，给非零特征模式分配等量功率的策略是渐进最优，即，规定$k$为非零$\lambda_i^2$的数量（$\mathbf{H}$的秩/该MIMO系统每秒每赫兹的空间自由度/传输信号经过MIMO信道之后信号（$\mathbf{Hx}$）的维度)，则$P_i^*=P/k$，代入有
$$
C\approx\sum_{i=1}^k\log\pqty{1+\frac{P\lambda_i^2}{kN_0}}\approx k\log\text{SNR}+\sum_{i=1}^k\pqty{\frac{\lambda_i^2}{k}}
$$
其中$\text{SNR}\triangleq P/N_0$。秩$k$或$n_\min$提供了关于MIMO系统空间自由度的简单度量。

​	进一步，由詹森不等式可以得到处理后的如下不等式
$$
\frac{1}{k}\sum_{i=1}^k\log\pqty{1+\frac{P}{kN_0}\lambda_i^2}\leq\pqty{1+\frac{P}{kN_0}\pqty{\frac{1}{k}\sum_{i=1}^k\lambda_i^2}}
$$
上式的等号仅在$\lambda_i=\lambda_j,\forall i,j$满足。考虑到$\sum_{i=1}^k\lambda_i^2=\tr\pqty{\mathbf{HH}^\mathrm{H}}\sum_{i,j}\abs{h_{i,j}}^2$可以视为信道矩阵的总功率增益，则通过上述结果可以说明，当SNR较大时，在相同的信道总功率增益下，该信道矩阵的所有特征值相同时能达到信道容量的上限，同时，若各特征信道不完全相同，则这些特征值越相近，能达到的信道容量就越大。$\pqty{\max_i(\lambda_i)/\min_i(\lambda_i)}$称为矩阵$\mathbf{H}$的状态数，可以衡量信道质量，当该数接近1时为最优。由此得到结论：***高信噪比时，信道质量较好能促进通信性能***。

​	· **低信噪比时**，功率分配的最优策略是只分配功率给较大特征值的特征信道，最终的信道容量为$C=\frac{P}{N_0}\pqty{\max_i\lambda_i^2}\log_2e$，此时只考虑MIMO信道最大特征值$\max_i\lambda_i^2$的一路特征信道，因此状态数和秩也不再重要，唯一重要的是有多少发射机能量成功到达了接收机。

###### 7.2 MIMO信道的物理建模

​	本节展示了MIMO的多路复用性能对物理环境的依赖，分析了一些理想状况下的信道矩阵秩和条件数。本节研究的都是均匀性天线阵列（天线均匀间隔分布于一条直线上）。

**视距SIMO信道**：单天线的发射机和多天线的接收机间只有一条无散射、无反射的自由空间直射信道，接收机各天线间的间距为$\Delta_r\lambda_c$，其中$\Delta_r$为归一化的接收天线间距系数（归一化为单位载波波长），$\lambda_c$为载波波长。天线规模远小于收发机之间的距离。

​	第$i$个接收机和发射机之间的信道相应可以表示为$h_i(\tau)=a\delta(\tau-d_i/c),\forall i$，其中$d_i$是第$i$个接收机到发射机之间的距离，$c$是光速，$a$为路径衰减（这里假设各条路径衰减相同）。假设$d_i/c\ll1/W$，$W$是传输带宽，则基带信道增益通过傅里叶变换可从时域变换到频域为
$$
h_i=a\exp\pqty{-j\omega\frac{d_i}{c}}=a\exp\pqty{-\frac{j2\pi f_cd_i}{c}}=a\exp\pqty{-\frac{j2\pi d_i}{\lambda_c}}
$$
其中$f_c$为载波频率。那么SIMO信道可以被重新表示为$\mathbf{y}=\mathbf{h}x+\mathbf{w}$，$x$是传输信号，$\mathbf{w}$是噪声，$\mathbf{y}$是接收信号向量，$\mathbf{h}=[h_1,\cdots,h_{n_r}]^\mathrm{H}$是信道增益向量，同时也可称为信号方向或发射信号再接收天线阵列上感应出的空间特征图。

​	又由于收发机之间的距离远大于阵列大小，则可以有$d_i\approx d+(i-1)\Delta_r\lambda_c\cos\phi,\forall i$，其中$d$是发射天线到第一个接收天线的距离，$\phi$是其入射角（也是相对于第一个接收天线），而$(i-1)\Delta_r\cos\phi$被称为第$i$个天线相较于第一个天线的距离增量，并且$\Omega\triangleq\cos\phi$被称为相对于接收天线的方向余弦。则代入$d_i$可以有
$$
\mathbf{h}=a\exp\pqty{-\frac{j2\pi d}{\lambda_c}}
\begin{bmatrix}
1\\ \exp\pqty{-j2\pi\Delta_r\cos\phi}\\ \vdots \\ \exp\pqty{-j2\pi\pqty{n_r-1}\Delta_r\cos\phi}
\end{bmatrix}
$$
为方便表示，这里定义单位空间特征图
$$
e_r\pqty{\Omega}\triangleq\frac{1}{\sqrt{n_r}}
\begin{bmatrix}
1\\ \exp\pqty{-j2\pi\Delta_r\cos\phi}\\ \vdots \\ \exp\pqty{-j2\pi\pqty{n_r-1}\Delta_r\cos\phi}
\end{bmatrix}
$$
这里$\sqrt{1/n_r}$的意义是

​	最佳接收机只是将含有噪声的接收信号投影到信号方向上，即最大比合并或接收波束成形。这样能够适应不同的接收时延从而使接受信号相长合并，提供$n_r$倍功率增益，最终得到信道容量为
$$
C=\log\pqty{1+\frac{P\norm{\mathbf{h}}^2}{N_0}}=\log\pqty{1+\frac{Pa^2n_r}{N_0}}\ \ \text{bits/s/Hz}
$$
​	需要注意的是，SIMO只提供了功率增益没有提供自由度增益。

**视距MISO信道**：与视距SIMO信号对称相似，有相同的$\mathbf{h}$和$e_r\pqty{\Omega}$，但是$y=\mathbf{h}^\mathrm{H}\mathbf{x}+w$。同样只有$n_r$的功率增益不存在自由度增益。

**单视距信道的MIMO阵列**：这里考虑只有单直射信道的MIMO天线阵列（收发天线阵列皆为线性阵列，且收发天线阵列中各天线间的间距分别被归一化为$\Delta_r$和$\Delta_t$），根据之前在研究SIMO信道时得到的结论，第$k$个发射天线到第$i$个接收天线间的信道增益可以被表示为
$$
h_{ik}=a\exp\pqty{\frac{-j2\pi d_{ik}}{\lambda_c}}
$$
其中$d_{ik}$是各天线间的距离，$a$是衰减因子（这里也假设各天线对之间的衰减因子相同）。不妨同样认为天线阵列尺寸明显小于收发天线间的距离，则$d_{ik}$可以被表示为
$$
d_{ik}=d+(i-1)\Delta_r\lambda_c\cos\phi_r-(k-1)\Delta_t\lambda_c\cos\phi_t
$$
其中$d$是线性收发天线阵列中的第一个发射天线到第一个接收天线间的距离，而$\phi_t$和$\phi_r$则分别表示为发射天线阵列和接收天线阵列上的视距信道入射角。为简化式子，分别定义$\Omega_t\triangleq\cos\phi_t$和$\Omega_r\triangleq\cos\phi_r$并代入$h_{ik}$得到
$$
h_{ik}=a\exp\pqty{\frac{-j2\pi d}{\lambda_c}}\cdot\exp\pqty{-j2\pi\pqty{i-1}\Delta_r\Omega_r}\cdot\exp\pqty{j2\pi\pqty{k-1}\Delta_t\Omega_t}
$$
进而可以导出信道矩阵$\mathbf{H}$为
$$
\mathbf{H}=a\sqrt{n_tn_r}\exp\pqty{\frac{-j2\pi d}{\lambda_c}}\mathbf{e}_r\pqty{\Omega_r}\mathbf{e}_t\pqty{\Omega_t}^*
$$
其中$\mathbf{e}_r\pqty{\cdot}$和$\mathbf{e}_t\pqty{\cdot}$则参考SIMO和MIS中的单位空间特征图，$\sqrt{n_tn_r}$是为了解决$\mathbf{e}_r$和$\mathbf{e}_t$自带的$\sqrt{n_r}$和$\sqrt{n_t}$。通过$\mathbf{e}_r\pqty{\Omega_r}\mathbf{e}_t\pqty{\Omega_t}^*$则可构造MIMO的单位空间特征图，便于理解可以认为是多个MISO或SIMO系统的叠加。显然$\mathbf{H}$秩为1，且有一个非零特征值为$\lambda_1=a\sqrt{n_tn_r}$，则该信道的信道容量为
$$
C=\log\pqty{1+\frac{Pa^2n_tn_r}{N_0}}\ \ \text{bits/s/Hz}
$$
​	需要注意的是，即使有多副收发天线组成天线阵列，该通信系统的DoF仍旧是一维的（由于信号最终投影在一维空间）。$n_tn_r$是功率增益因子，单发射天线时，接收机采用接收波束成形即可获得增益；单接收天线时，采用发射波束成形即可获得增益；多发射天线和多接收天线下，系统可同时受益于发射波束成形和接收波束成形。最终，视距信道下，MIMO系统仅提供了功率增益，未提供DoF增益。

**地理上离散的天线阵列**：之前较为密集的MIMO阵列都无法取得DoFs增益，这里研究地理上离散的天线阵列以获得DoFs的增益。假设两发射天线，其间距与其与接收天线阵列的距离可比拟，且他们到接收天线阵列都只有一条直射（LoS）信道，两条直射信道分别有衰减因子和入射角为$a_1$，$\phi_{r1}$和$a_2$，$\phi_{r2}$。进一步假设发射天线的信号时延扩展远小于$1/W$，则可以继续使用上文的建模，得到发射天线$k$发射信号到接收天线阵列的空间特征图为
$$
\mathbf{h}_k=a_k\sqrt{n_r}\exp\pqty{\frac{-j2\pi d_{1k}}{\lambda_c}}\mathbf{e}_r\pqty{\Omega_{rk}},k=1,2,
$$
其中$d_{1k}$是第$k$个发射天线到接收阵列中的第1个接收天线的距离，$\Omega_{rk}$和$\mathbf{e}_r\pqty{\cdot}$参考之前的定义即可。

观察$\mathbf{e}_r\pqty{\Omega_{rk}}$的原函数向量，除第一个元素外，每一个元素的周期应该为$\mathrm{T}=2\pi/2\pi\pqty{n_{rk}-i}\Delta_r=1/\pqty{n_{rk}-i}\Delta_r$，因此，其中周期最大的为$\mathrm{T}_\max=1/\Delta_r$，最小的为$\mathrm{T}_\min=1/\pqty{n_{rk}-1}\Delta_r$，由于要考虑$\mathbf{e}_r\pqty{\cdot}$整个的周期，函数向量内各元素的周期最小公约数为$1/\Delta_r$，则$\mathbf{e}_r\pqty{\cdot}$以$\Omega_{rk}$为变量的函数的周期应该为$1/\Delta_r$（这里强调以$\Omega_{rk}$为变量是为便于后续求解两个阵列的入射角之差是否为$1/\Delta_r$周期的倍数，排除空间特征向量重叠的可能性）。当满足$\Omega_r\triangleq\Omega_{r2}-\Omega_{r1}\neq0,\text{mod}\frac{1}{\Delta_r}$时， $\mathbf{H}=[\mathbf{h}_1,\mathbf{h}_2]$则各列不同且线性无关，从而可以得到$\mathbf{H}$有两个非零特征向量$\lambda_1^2$，$\lambda_2^2$，产生了两个自由度。

尽管满足$\Omega_r\neq0$后得到的$\mathbf{H}$已经是满秩矩阵，但是该矩阵仍可能为病态矩阵（矩阵系数的微小变化对最后的求解影响很大），因此还需从方向角的距离进一步研究实现$\mathbf{H}$良矩阵的必要条件，以实现空间自由度的高效利用。

$\mathbf{H}$矩阵的好坏取决于两个发射天线的空间特征图的一致性，一致性越低则$\mathbf{H}$状态越好。两个空间特征图间的差角$\theta$满足
$$
\abs{\cos\theta}\triangleq\abs{\mathbf{e}_r\pqty{\Omega_{r1}}^*\mathbf{e}_r\pqty{\Omega_{r2}}}
$$
注意到$\mathbf{e}_r\pqty{\Omega_{r1}}^*\mathbf{e}_r\pqty{\Omega_{r2}}$只取决于$\Omega_r$，从而可以定义
$$
\begin{align}
f_r\pqty{\Omega_r}
&\triangleq\mathbf{e}_r\pqty{\Omega_{r1}}^*\mathbf{e}_r\pqty{\Omega_{r2}}\\
&=\frac{1}{n_r}\sum_{k=0}^{n_r-1}\exp\pqty{j2\pi\Delta_r\Omega_rk}\\
&=\frac{1}{n_r}\exp\pqty{j\pi\Delta_r\Omega_r\pqty{n_r-1}}\sum_{k=1,3,5,\cdots}^{2n_r-1}\exp\pqty{j\pi\Delta_r\Omega_r\pqty{k-n_r}}\\
&=\frac{1}{n_r}\exp\pqty{j\pi\Delta_r\Omega_r\pqty{n_r-1}}\frac{\exp\pqty{j\pi\Delta_r\Omega_r\pqty{1-n_r}}\pqty{1-\exp\pqty{2j\pi\Delta_r\Omega_rn_r}}}{1-\exp\pqty{2j\pi\Delta_r\Omega_r}}\\
&=\frac{1}{n_r}\exp\pqty{j\pi\Delta_r\Omega_r\pqty{n_r-1}}\frac{\exp\pqty{-j\pi\Delta_r\Omega_rn_r}-\exp\pqty{j\pi\Delta_r\Omega_rn_r}}{\exp\pqty{-j\pi\Delta_r\Omega_r}-\exp\pqty{j\pi\Delta_r\Omega_r}}\\
&=\frac{1}{n_r}\exp\pqty{j\pi\Delta_r\Omega_r\pqty{n_r-1}}\frac{\sin\pqty{\pi\Delta_r\Omega_rn_r}}{\sin\pqty{\pi\Delta_r\Omega_r}}\\
&=\frac{1}{n_r}\exp\pqty{j\pi\Delta_r\Omega_r\pqty{n-1}}\frac{\sin\pqty{\pi L_r\Omega_r}}{\sin\pqty{\pi L_r\Omega_r/n_r}}
\end{align}
$$
其中$L_r\triangleq n_r\Delta_r$是归一化的接收天线阵列的长度。因此
$$
\abs{\cos\theta}=\abs{f_r\pqty{\Omega_r}}=\frac{1}{n_r}\abs{\frac{\sin\pqty{\pi L_r\Omega_r}}{\sin\pqty{\pi L_r\Omega_r/n_r}}}
$$
信道矩阵$\mathbf{H}$的状态取决于$\cos\theta$，这里不妨假设衰减因子相同（或接近）$a_1=a_2=a$，则矩阵$\mathbf{H}$的奇异平方为
$$
\lambda_1^2=a^2n_r\pqty{1+\abs{\cos\theta}},\lambda_2^2=a^2n_r\pqty{1+\abs{\cos\theta}}
$$
条件数则可相应地表示为$\frac{\lambda_1}{\lambda_2}=\sqrt{\frac{1+\abs{\cos\theta}}{1-\abs{\cos\theta}}}$，从中可以看出，当$\abs{\cos\theta}\approx1$时，信道矩阵$\mathbf{H}$为病态矩阵，否则为良矩阵。进一步分析，当$\cos\theta=1$时，$\abs{\frac{\sin\pqty{\pi L_r\Omega_r}}{\sin\pqty{\pi L_r\Omega_r/n_r}}}$需要为$n_r$，根据取极限的方式，该式只有在$\pi L_r\Omega_r$趋近于0，或者直接的$\Omega_r=\Omega_{r2}-\Omega_{r1}$趋近于0时满足该式为$n_r$，考虑到实际中$\Omega_{r1}$和$\Omega_{r2}$的取值范围有限，所以可以基本说明，一定范围内，$\Omega_{r1}$和$\Omega_{r2}$差异越大时，信道矩阵$\mathbf{H}$状态越好。





#### · 《6G无线传输技术》学习笔记

##### 第七章 大规模MIMO

###### 1. XL-MIMO容量分析

对于一个MIMO传输系统（基站天线数量$M$，单用户单天线数量$K$），其上行传输的信号模型分别为
$$
\mathbf{y=Gx+z}
$$
其中$\mathbf{G}\in\mathbb{C}^{M\times K}$表示的是MIMO信道矩阵，$\mathbf{x}=[x_1,x_2,\cdots,x_K]^\mathrm{T}\in\mathbb{C}^{K\times 1}$，$\mathbf{z}\in\mathbb{C}^{M\times 1}$即是噪声。

根据香农公式
$$
C=\underbrace{\log_2\det\pqty{\mathbf{I}_K+\frac{\rho_d}{M}\mathbf{G}^\mathrm{H}\mathbf{G}}}_{下行通信速率}=\underbrace{\log_2\det\pqty{\mathbf{I}_M+\frac{\rho_d}{M}\mathbf{GG}^\mathrm{H}}}_{上行通信速率}
$$
$\rho_d$表示发射端SNR。因此，容量上下界可以表示为
$$
\log_2\pqty{1+\frac{\rho_d\tr\pqty{\mathbf{GG}^\mathrm{H}}}{M}}\leq C\leq\min\pqty{M,K}\log_2\pqty{1+\frac{\rho_d\tr\pqty{\mathbf{GG}^\mathrm{H}}}{M\min\pqty{M,K}}}
$$



---

### 计算机基础

#### · Markdown

​	Markdown语法是一种非常实用的标记性语言，对于日常科研和论文撰写来说，能够非常良好地充当论文、教材阅读和论文写作的桥梁，相较正常的文字输入，Markdown向我们提供了更加多样且便捷的文本输入方式，同时还能兼容LaTeX的公式语法，这里以Typora为例展示Markdown的部分实用语法（需要注意其中一些快捷键无法在Typora之外的软件正常使用）。

1. 特殊文字：**加粗**、*斜体*、<u>下划线</u>，***<u>综合</u>***
2. 上下标和高亮：正文^上标^、正文~下标~、==高亮==
3. 无序列表：-/+/*（空格）对应分别不同的无序列表，用Tab键区分无序列表分级
4. 有序列表：数字+“.”+（空格），同样使用Tab分级
5. 任务列表："- [x]" 任务一（此为勾选状态）
6. 行内代码和公式：代码`int`，公式$E=mc^2$
7. 段落代码和公式：可以通过行内公式进行引用$(1.1)$

```c++
## 代码
int main()
{
    cout<<"Hello Wordl!"<<endl;
    return 0;
}
```

$$
\text{公式}\\
C=\log\pqty{1+\frac{P}{\sigma^2}}\tag{1.1}
$$

8. 表格：使用Ctrl+T快捷键打开表格设置行列即可，每一列同样可以设置对齐方式
9. 引用和注释：段落引用">"+（空格）

> 这一段作为引用展示[^1]

[^1]: 引用对应注释

10. 链接：[Markdown语法学习视频]([12分钟学会Markdown｜附Typora使用方法_哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1Fg411j7CW/))
11. 目录：在任意地方单开一段，然后输入“[TOC]”

#### · Missing Class

**Shell**

​	在Shell中使用echo可以打印出参数中的内容，例如

```shell
echo Hello
>>> Hello
```

​	Shell中的路径是一组被分割的目录，在Linux和macOS上使用`/`分割，且根目录为`/`；而Windows上使用`\`分割，且每个盘都有一个根目录`C:\`。若为Linux系统，则以`/`开头是为绝对路径，其他情况下都是相对路径（相对于当前工作目录下的路径，当前工作目录可以由`pwd`命令获取）。切换目录时使用`cd`命令，`.`代表本级目录（当前工作目录），`..`代表上级目录。查看当前工作目录下包含哪些文件使用`ls`命令。Shell中都会实时显示当前工作目录。

**注意**：*当未明确指定目录时，目录通常在当前工作目录下运行，因此，在执行程序时应当确认是否希望在当前目录下执行。*

```shell
pwd
>>> C:\Users\Simon\documents\CppProjects\task0
------------------------------------
ls
>>>	helloworld.cpp
	helloworld.exe

ls -l [path]
drwxr-xr-x 1 missing  users  4096 Jun 15  2019 missing
------------------------------------
cd ..
C:\Users\Simon\documents\CppProjects
```

​	`ls -l [path]`命令可以给出指定路径下的文件列表，且文件信息 也更加详细。`d`代表是一个目录，后边紧接着为9个字符，每三个字符`rwx`一组，分别代表读取权限、写入权限和搜索权限，若为`-`则对应用户无该权限。读取权限，能否获得文件夹内文件列表，读取文件夹内信息；写入权限，能否操作文件夹内文件（增删）；搜索权限，能否看到并进入该文件夹进行访问。

​	在对文件整体进行操作的过程中，还有一些命令需要关注

```shell
rm -i -r -f [-path]	//-i逐一确认，-r逐一删除，-f强制删除
mv [-path] [-destination]	//可以移动文件到指定位置，或完成改名
cp [-path] [-destination]	//复制文件到指定位置（-r用于目录复制）
mkdir -p [-path]	//创建指定路径的目录（-p判断是否存在，若不存在即创建）
```

​	其中具体的文件可以用通配符`*`代替用以批量执行命令。其他不常用的命令行可以通过`man`命令获取具体使用方法和意义。

​	在Shell中主要有输入流和输出流，而在终端界面中，键盘是输入，显示器是输出，但是可以通过`<`和`>`分别重定向



#### · Git

​	一个开源的分布式版本控制软件，可以接入GitHub、Gitee等远程项目托管网站。可以从中建立仓库(repository)，并管理仓库中各个文件的版本，可以方便地存储和回退。需要下载专门的Git关键并完成部署（可以通过终端、IDE端和Git自带的GUI进行控制和管理）

**初始化配置**

```
git config --global user.name/email Name/@gmail.com	//配置用户名/邮箱
git config --global --list	//查看配置信息
```

​	`git config`只需在最开始的时候设置一次作为作者信息记录，后续不需要重复设置，但是一定要在开始的时候记得设置作者签名。

**Repo创建**

```
//本地repo创建
git init	//将当前文件夹作为repo
git init [foldername]	//在当前目录创建一个文件夹foldername，并作为repo

//Clone远程repo
git clone [website]/[SSHcode]	//下载对应项目成为子文件夹，并作为一个repo
```

​	运行上述代码的同时还会在repo所在的文件夹中生成一个隐藏的`.git`文件夹，文件夹中为repo的必要软件，不能随意删除（删除`.git`文件夹后，该repo随即退化为普通文件夹）。

**添加和提交**

```
git status	//可以获取该仓库的状态（哪些仓库未被Git管理，哪些在暂存中）
git log	//查看提交记录（包括作者和邮箱）
git add [filename]	//将指定文件filename添加进入暂存中
git commit -m "Operation"	//提交到暂存区，并命名操作为"Operation"以区分
git commit -am/(-a -m) "Operation"	//同时添加暂存并完成提交到仓库
git rm --cached [filename]	\\将暂存区中的指定文件filename从暂存区中移除
```

​	对于`git add`命令，可以用文件后缀作为通配符取代`filename`将文件批量添加进入暂存区，同时还可以添加一整个文件夹（例如`git add .`可以将当前文件夹下的所有文件都添加进入暂存区）。

**版本回退**

```
git reset [harsh number]/HEAD^	//指定回退的版本号，HEAD^表示上一个版本
git reset --hard/soft/mixed	//分别代表三种模式的版本回退
```

​	`git reset`是进行版本控制回退的重要代码，可以设置三个不同的参数，分别有不同的作用

· `--mixed`：可以缺省，默认模式的reset，不删除本地的工作区，但清空本地的暂存区

· `--soft`：不可缺省，同时保留本地的工作区和暂存区，只将版本回退到指定版本

· `--hard`：不可缺省，同时删除指定版本之后在工作区和暂存区产生的一切文件

​	这一部分的关键是要理解工作区和暂存区的意义，工作区就是文件夹的区域，文件夹中能看到的所有文件都位于工作区中。`--mixed`和`--soft`在版本回退后都会保留工作区中的内容，也就是只将Git的repo中的版本完成回退（保存的内容进行删除或退回删除），实际的本地文件夹（工作区）不会受到影响，由于这一过程执行的版本回退可以再重新提交。而`--hard`则会删除暂存区和工作区中多余的文件，使其与repo中保持严格一致。因此，最常用的就是默认的`--mixed`（或直接缺省）

**查看差异**

```
git diff	//默认状态下查看工作区和暂存区中对应文件的差异
git diff HEAD	//加入关键字HEAD，比较工作区和repo中最后一个版本中文件差异
git diff --cached	//比较暂存区和repo版本库中的差异
git diff [version1] [version2]	//比较version1和version2两版本的文件差异
git diff HEAD~(n) HEAD	//比较当前版本和上n个版本间的文件差异
git diff HEAD^ HEAD [filename]	//查看指定文件filename在两版本之间的差异
```

**删除文件**

```
git rm [filename]	//同时删除工作区和暂存区的文件filename（但是还需提交）
git rm --cached [filename]	//仅从暂存区中删除
```

**忽略文件**

`.gitignore`文件用于将那些不用上传到repo中的文件在上传中忽略（例如账户信息、密码、日志或可以由其中文件生成的文件）。在Linux中使用以下代码完成所需忽略的文件的添加

```
echo [filename] > .gitignore	//添加单个文件
echo *.[suffix] > .gitignore	//利用通配符*添加一类文件
```

`.gitignore`本质上是一个文本文档，可以向其写入所需忽略的文档名或利用通配符规则表述所需忽略的文档特征。*匹配不定长字符串，?匹配单个字符，[]匹配括号中出现的任意字符（[]中可以用“1-9”或“1/2/3”表示），!表示取反。

**关联本地和远程仓库**

可以将本地已经有的和远程已经有的仓库关联起来，不必再每次都使用`git clone`。

```
git remote add [ProjName] [Github_SSH_address]	//将该目录与远程进行关联
git remote -v	//查看远程仓库的别名和地址
git branch -M main/[BranchName]	//指定分支名称
git push -u [ProjName] [BranchName1]:[BranchName2]	//关联本地和远程分支
git pull [ProjName] [BranchName1]:[BranchName2]	//拉取
```

​	关联本地和远程的分支意味着上传或拉取，这里的`[ProjName]`指的是在添加远程仓库时的命名的别名，用以找到对应的远程仓库，`[BranchName]`用于本地分支和远程分支对齐。

**在VsCode中使用Git**

​	VsCode是常用的代码编辑器，可以通过VsCode的Git管理功能管理文件，同样可以分为工作区、暂存区和仓库区，还可以直接上传到线上。VsCode可以直接自动跟踪文件的修改和添加。

**分支管理和操作**

​	分支管理便于大团队项目协作的场景，可以解决多用户冲突的问题，每个分支都可以形成独立的项目。

```
git branch	//查看当前项目的分支（*指示当前所在的分支）
git branch [BranchName]	//创建名为[BranchName]的分支
git switch/checkout [BranchName]	//切换分支（优先switch)
git merge [BranchName]	//合并[BranchName]到当前分支（合并前确认当前分支）
git branch -d [BranchName]	//删除已经完成合并的[BranchName]分支
git branch -D [BranchName]	//强制删除未合并的[BranchName]分支
git merge --abort	//遇到合并冲突后放弃合并
```

​	当合并分支时，如果合并的文件中有冲突（两个分支中对部分文件进行了不同的修改）需要首先解决冲突（打开冲突文件，手动修改文件内容进行合并）再完成合并，可以使用`git status`和`git diff`查看冲突的具体情况（也可以选择放弃合并）。

​	需要注意的是，完成分支合并之后，原来的分支依旧存在，只是在分支中创建和修改的文件或内容被批量提交到了主线中。

```
git rebase [BranchName]	//将当前分支变基到[BranchName]分支上
```

​	变基操作是找到当前分支和目标分支的共父子节点，然后将当前分支嫁接到目标分支的HEAD指针后。这里需要注意是哪个分支嫁接到哪个分支上。

​	使用`merge`可以很好保留项目的历史情况，说明提交记录和各分支原本的情况，但是分支图会比较复杂；使用`rebase`能够让项目更加线性和直观，但是会改变历史关系，应避免在共享分支时使用。

#### · GitHub



#### · SSH

​	SSH是一种远程安全数据传输协议，需要现在Windows系统上安装OpenSSH客户端和服务器端。检验是否安装成功并且正确运行ssh

```
Get-Service -Name *ssh*	//抓取是否正在运行SSH（确保自动启动）
netstat -an | findstr : 22	//确认是否正确监听22端口
```

**内网连接**

​	内网连接适用于处于同一局域网或者目标对象为公网IP的连接设备

```
ssh "UserName"@IP_address	//连接该IP地址中的用户"UserName"
```

这一代码可以在终端中直接使用这样就可以跳转到对应主机的终端，然后利用命令行完成操作；也可以在VsCode上安装Remote -ssh插件后，添加远程主机输入上述代码，可访问的远程主机可在VsCode中统一管理（还可修改主机在本地的命名），通过VsCode打开对应主机，就可以更加清晰直观借助GUI对文件进行操作。

​	为了实现便捷的免密访问，可以利用ssh的密钥机制。

```
ssh-keygen -t rsa -b 4096	//生成-t类型-b长度的一对密钥
ssh-copy-id -i [file_address] "UserName"@IP_address	//-i指定公钥文件地址
```

运行第一行代码会生成一个公钥一个私钥密钥对，同时会提示你进行命名，完成命名后会设置密码用于私钥的使用。最终`.pub`文件为公钥，另一个无后缀的为私钥。可以通过`cat`（windows）或`vi`（Linux）完成公钥的读取。注意，私钥是不能随意泄露的。部署好了ssh公钥到目标主机后就可以免密登录访问目标主机文件并进行操作。

**非内网连接**

​	适用于进行连接的两台主机不在同一局域网下，且双方均无公网IP，但是依旧有远程进行文件管理和操作的需求。这里采用frp进行内网穿透。

```
sudo vim etc/hosts	//MacOS超级管理员模式下的Linux（打开后加入桥服务器IP）
echo [ip_address] >> C:\Windows\System32\drivers\etc\hosts	//Windows
```

根据请求联网主机的系统版本和类型不同，运行上述代码打开不同系统下的`hosts`文件，在`hosts`文件中添加桥服务器的公网IP地址（还可在添加桥服务器公网IP地址时通过`[ip_address] alias`设置该IP地址的别名，便于之后使用）

​	*1.* 从github中下载开源软件frp，根据桥服务器和目标主机（服务器）的系统类型下载对应版本的frp，其中frps部分为桥服务器的程序及配置文件；frpc为目标主机的程序及配置文件（在不同设备上保留必要的文件，其余可删去）。*2.* 在配置文件中，大部分保留原有设置即可，需要时可以再读配置代码修改，这一过程需要注意和修改的只有frpc.ini文件中server_addr改为桥服务器的公网IP。*3.* 将完成修改的配置文件分别拷贝到对应设备，并根据配置文件开放云服务器上的对应端口（设置为TCP协议）

​	接下来需要分别启动桥服务器和目标主机上的frp程序。在分别执行以下代码之前需要将命令行转至frps或frpc文件夹

```
./frps -c ./frps.ini	//启动桥服务器端的frps程序
./frpc -c ./frpc.ini	//启动目标主机端的frpc程序
```

上述命令可能由于不同系统有细微区别，但是大致相同。配置完成后，在使用中，如果希望（经由桥服务器）连接目标主机进行操作，则需要执行以下代码

```
ssh -p [in_port] "UserName"@IP_address	//连接目标主机
```

运行上述代码即可完成连接。需要注意的是，其中`in_port`是指定端口，是在桥服务器端设置的入口（也就是桥服务器的配置文件说明的数据接收端口），`UserName`是目标主机的用户名，而`IP_address`是桥服务器的公网IP。

​	上述过程只满足了手动启动，通过systemd中的文件可以完成自动启动。找到frps文件夹中的frps.service复制到`/lib/systemd/system/frps.service`，并将`User`一项改为当前用户的用户名（桥服务器的用户名）；`ExectStart`一项修改为`frps -c frps.ini`命令，且其中的文件应为文件的绝对路径（`pwd`可查看绝对路径）。之后运行以下代码

```
systemctl daemon-reload
systemctl enable frps
systemctl start frps
systemctl status frps
```

运行上述代码期间可能需要输入多次密码，待最后一次代码执行完成后，若找到`Active:active`即配置桥服务器自启动成功。

如果目标主机为Linux系统，则将上述frps改为frpc即可。

如果目标主机为Windows系统，则建立一个内容如下的`.bat`文件

```
@echo off

if not defined TAG{
	set TAG=1
	start wt -p "cmd" %0
	exit
}
:home
frpc -c frpc.ini
goto home
```

并将上述`.bat`文件的快捷方式放置在`C:\Users\[Your_User_Name]\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup`中即完成Windows系统的目标主机的自启动。



---

### 深度学习

#### · 《李沐动手学深度学习》学习笔记

一个完整的故事：数据科学家（数据$\rightarrow$模型训练）+AI专家（更好的模型）+领域专家（分析用户行为）

找paper的经验

 
