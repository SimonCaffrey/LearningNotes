### Paper Reading Notes

#### · CAPA

##### · *Beamforming Optimization for Continuous Aperture Array (CAPA)-based Communications*

###### ·Abstract

​	Under LoS condition, [25] investigated the SINR of uplink multi-user CAPA systems and  proposed an adaptive interference mitigation method. While, [27] explored a FT-based beamforming optimization in downlink multi-user CAPA system to maximize the sum rate. [26] also studied the capacity of the uplink multi-uesr CAPAs equipped at both transceivers and the downlink ones duality.

###### ·System Model

​	Define the BS as a surface with area of $S_t$, and the $\pmb{J}(s,\omega)$ represent the FT of source current density at point s with an angle of $\omega=2\pi f/c=2\pi /\lambda$, where $s\sub S_t$ and $f, \lambda$ denote frequency and wavelength respectively. 

​	If only single carrier, $\omega$ can be omitted

​	if only vertically polarized transmitter, vectors of others will be also omitted, except $\pmb{y}$.

​	After a series of simplification, we finally know
$$
\pmb{J}(s)=J_y(s)\hat{u_y}=J(s)\hat{u_y}
$$
​	considered the linear addition of k users in signals, the scaler
$$
J(s)=\sum_{k=1}^\boldsymbol{K} J_k(s)x_k
$$
​	where $J(s)$ is current and $x$ is communication symbol (0/1) for each users.

​	Having known $J(s)$ is the source current density of BS surface, the EF at $r_k$ is given due to Maxwell function:
$$
\pmb{E_k}=\int_{S_T}\pmb{G}(r_k, s)\pmb{J}(s)ds\\
\mathbf{G}(\mathbf{r}, \mathbf{s})=-\frac{j \eta e^{-j \frac{2 \pi}{\lambda}\|\mathbf{r}-\mathbf{s}\|}}{2 \lambda\|\mathbf{r}-\mathbf{s}\|}\left(\mathbf{I}_{3}-\frac{(\mathbf{r}-\mathbf{s})(\mathbf{r}-\mathbf{s})^{T}}{\|\mathbf{r}-\mathbf{s}\|^{2}}\right)
$$
​	where $\eta$ is intrinsic impedance, $\pmb{G}(r,s)$ can be modeled as a stochastic process in scatter env.

​	Then, there is  a uni-polarized antenna for each user, which means only single direction signal can be captured, if the direction vector is $\hat{u_k^T}$, the EF (having considered the noise) received will be:
$$
E_k=\hat{\pmb{u}_k^T}\pmb{E_k}+n_k=\int_{S_T}\hat{\pmb{u}_k^T}\pmb{G}(r_k, s)\pmb{J}(s)ds+n_k\\
$$
​	Naturally, consider the transmission, extracting a system function for each user $k$ from the formula above
$$
{H}_k(s)=\hat{\pmb{u}_k^T}\pmb{G}(r_k,s)\hat{\pmb{u}_y}
$$
Then, we can know that EF received by each user is:
$$
\begin{flalign}
E_k=
&\int_{S_T}H_k(s)\sum_{i=1}^KJ_i(s)x_ids+n_k\\
=&\underbrace{\int_{S_T}H_k(s)J_k(s)x_kds}_\text{interest signal}+\underbrace{\sum_{i=1,i\neq k}^K\int_{S_T}H_k(s)J_i(s)x_ids}_{\text{interference}}+\underbrace{n_k}_\text{noise}
\end{flalign}
$$
​	Apparently, $E_k$ is divided into three parts: the signal user $k$​ want, the signal from BS to others, the noise.

​	Now, we directly investigate the SINR.

​	Firstly, we calculate the Power received
$$
P_k=\mathbb{E}\{\varepsilon_k\times \frac{|E_k|^2}{2}\times\frac{1}{\eta}\}
$$
​	where $\eta$ means intrinsic impedance and $\varepsilon_k$ means absorption efficiency at user $k$.

​	We already know $E_k$, the $P_k$ can be calculated from $E_k$.
$$
\begin{flalign}
P_k=
&\mathbb{E}\{\frac{\varepsilon_k|E_k|^2}{2\eta}\}\\
=&\mathbb{E}\{\frac{\varepsilon_k}{2\eta}(\sum_{i=1}^K\int_{S_T}H_k(s)J_i(s)x_ids+n_k)\times (\sum_{j=1}^K\int_{S_T}H_k(s)J_j(s)x_jds)+n_k\}\\
=&\frac{\varepsilon_k}{2\eta}\sum_{i=1}^K\sum_{j=1}^K\mathbb{E}\{x_ix_j^*\}(\int_{S_T}H_k(s)J_i(s)ds)\times (\int_{S_T}H_k^*(s)J_j^*(s)ds)+\frac{\varepsilon_k}{2\eta}\sigma_k^2\\
=&\frac{\varepsilon_k}{2\eta}\sum_{i=1}^K|\int_{S_T}H_k(s)J_i(s)ds|^2+\frac{\varepsilon_k}{2\eta}\sigma_k^2
\end{flalign}
$$
​	Thus, the SINR is:
$$
\gamma_{k}=\frac{\left|\int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{k}(\mathbf{s}) d \mathbf{s}\right|^{2}}{\sum_{i=1, i \neq k}^{K}\left|\int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{i}(\mathbf{s}) d \mathbf{s}\right|^{2}+\sigma_{k}^{2}}
$$
​	And the rate will be $R_k=log_2(1+\gamma_k)$​

​	Finally, the problem can be expressed as:
$$
\max_{\{J_k(s)\}_{k=1}^K}\sum_{k=1}^K\alpha_kR_k=\alpha_klog_2(1+\gamma_k)\\
s.t \ \ \sum_{k=1}^K\int_{\mathcal{S}_T}|J_k(s)|^2ds=P_T
$$

​	where $\alpha_k$ means the weight of each user, and we intend to maximize the sum-rate (which means maximize this formula) under the limited power cost in BS.

​	After all, such a problem is a continuous optimization problem different from the one of MIMO (discrete) and massive formula contains integrals. FT do cope with this problem but not very well, which transfer this continuous problem to a discrete one, with the cost of massive basic function and low degree of approximation. So, we propose another approach combined with CoV and Corr-ZF, which directly optimize continuous problem.

###### ·CoV-BCD Approach

​	To better deal with this optimization problem, we rewrite this problem as:
$$
\max_{\{J_k(s)\}_{k=1}^K}\ \sum_{k=1}^K\alpha_klog_2(1+\gamma_k')
$$
​	where $\gamma_k'$ is a more convenient expression make our optimize more easily. Meanwhile, we denote $\{\overline{J_k'(s)}\}_{k=1}^K$ as the optimal solution to $\gamma_k'$ and above problem as well.
$$
\gamma_k'=\frac{|\int_{\mathcal{S}_T}H_k(s)J_k'(s)ds|^2}{(\sum_{i=1,i\neq k}^K|\int_{\mathcal{S}_T}H_k(s)J_i'(s)ds|^2+\sigma_k^2\ \frac{\sum_{i=1}^K\int_{\mathcal{S}_T}|J_i'(s)|^2ds}{P_T})}
$$

​	And the optimal solution to this rewritten problem can be expressed as:
$$
\overline{\gamma_k'}=\frac{|\int_{\mathcal{S}_T}H_k(s)\overline{J_k'(s)}ds|^2}{(\sum_{i=1,i\neq k}^K|\int_{\mathcal{S}_T}H_k(s)\overline{J_i'(s)}ds|^2+\sigma_k^2\ \frac{\sum_{i=1}^K\int_{\mathcal{S}_T}|\overline{J_i'(s)}|^2ds}{P_T})}
$$
​	ps: why we adopt such a formula as our optimization formula?   To find the optimal ratio rather than higher power demand. In this formula,the value of $\gamma_k'$ will not change with the ratio of $J_i(s)$, which means easier to find the optimal ratio of each $J_i(s)$ rather than demanding power without considering limits.

​	Having know the optimal ratio of $J_i(s)$ for each user, we can multiple it with the appropriate coefficient under power constraint:
$$
\overline{J_k(s)}=\sqrt{\frac{P_T}{\sum_{k=1}^K\int_{\mathcal{S}_T}|\overline{J_k'(s)}|^2ds}}\overline{J_k'(s)}
$$
​	where $J_k(s)$ and $J_k'(s)$ represent the current density of original and rewritten problem, respectively. And in such condition satisfy $\overline{\gamma_k'}=\overline{\gamma_k}$, which prove formula $\overline{\gamma_k'}$​ is a better approach to find optimal ratio.

​	With the guide of ***Fractional programming for communication systems—part I: Power control and beamforming*** and ***part II: Uplink scheduling via matching*** or [无线通信中的分式规划（Fractional Programming）与Matlab实现](https://www.zhihu.com/tardis/zm/art/599204238?source_id=1005), the objective function is eventually derived as:
$$
\max _{\left\{\mu_{k}, \lambda_{k}, J_{k}(\mathbf{s})\right\}_{k=1}^{K}}  \sum_{k=1}^{K} \alpha_{k}\left(2 \mu_{k} \Re\left\{\lambda_{k}^{*} \int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{k}(\mathbf{s}) d \mathbf{s}\right\}\right.-\left|\lambda_{k}^{2}\right|\left(\sum_{i=1}^{K}\left|\int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{i}(\mathbf{s}) d \mathbf{s}\right|^{2}\right.\left.\left.+\frac{\sigma_{k}^{2}}{P_{\mathrm{T}}} \sum_{i=1}^{K} \int_{\mathcal{S}_{\mathrm{T}}}\left|J_{i}(\mathbf{s})\right|^{2} d \mathbf{s}\right)\right)
$$
​	Having eliminated the constraints and fraction and even decoupled the optimization variables, BCD is used to alternately optimize two blocks: $\{\mu_k,\ \lambda_k\}_{k=1}^K$ and $\{J_k(s)\}_{k=1}^K$​​.

1. $\{\mu_k,\ \lambda_k\}_{k=1}^K$: derived by above two papers into:

$$
\begin{flalign}
&\mu_k=\sqrt{1+\gamma'_k}\\
&\lambda_k=\frac{\mu_k\int_{\mathcal{S}_T}H_K(s)J_k(s)ds}{\sum_{i=1}^K(|\int_{\mathcal{S}_T}H_k(s)J_i(s)ds|^2+\frac{\sigma_k^2}{P_T}\int_{\mathcal{S}_T}|J_i(s)|^2ds)}
\end{flalign}
$$

2. $\{J_k(s)\}_{k=1}^K$: to manage this problem with CoV, we rewrite it as:

$$
g\left(J_{k}\right)=2 \Re\left\{A_{k} \int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{k}(\mathbf{s}) d \mathbf{s}\right\}-\sum_{i=1}^{K}\left(B_{i}\left|\int_{\mathcal{S}_{\mathrm{T}}} H_{i}(\mathbf{s}) J_{k}(\mathbf{s}) d \mathbf{s}\right|^{2}+C_{i} \int_{\mathcal{S}_{\mathrm{T}}}\left|J_{k}(\mathbf{s})\right|^{2} d \mathbf{s}\right),
$$

​	This function is a separable optimization for each $J_k(s)$, which is maximized by finding optimal $J_k(s)$ through CoV and where $A_k$, $B_i$, $C_i$​ satisfy following formula:
$$
\begin{cases}
A_k=\alpha_k\mu_k\lambda_k^*\\
B_i=\alpha_i|\lambda_i|^2\\
C_i=\frac{\alpha_i|\lambda_i|^2\sigma_i^2}{P_T}
\end{cases}
$$

​	***Lemma 3***: For all functions $U$ define on an open set $\mathcal{S}$ in complex space and satisfy $(U(s)=0, \forall s \in\partial\mathcal{S})$, if $\Re\{\int_{\mathcal{S}}U^*(s)V(s)ds\}=0$, then $V(s)=0, \forall s\in\mathcal{S}$.

​	***Propostion 1***: the optimal $J_k(s)$ has the following structure:
$$
J_{k}(\mathbf{s})=A'_{k} H_{k}^{*}(\mathbf{s})-\sum_{i=1}^{K} B'_{i} H_{i}^{*}(\mathbf{s}) \int_{\mathcal{S}_{\mathrm{T}}} H_{i}(\mathbf{z}) J_{k}(\mathbf{z}) d \mathbf{z}
$$
​	Where $A'_k=A_k/(\sum_{i=1}^KC_k)$, $B'_k=B_k/(\sum_{i=1}^KC_k)$ and the optimal $J_k(s)$ is also the solution to the above formula. What needs to be noticed is $\mathbf{z}$ in$\int_{\mathcal{S}_{\mathrm{T}}} H_{i}(\mathbf{z}) J_{k}(\mathbf{z}) d \mathbf{z}$ rather than $\mathbf{s}$.

​	Rewrite it to better solve it:
$$
J_{k}(\mathbf{s})=A'_{k} H_{k}^{*}(\mathbf{s})-\sum_{i=1}^{K}w_{k,i}\ B'_{i} H_{i}^{*}(\mathbf{s})
$$
​	Where we define the $w_{k,i}$ as $ \int_{\mathcal{S}_{\mathrm{T}}} H_{i}(\mathbf{z}) J_{k}(\mathbf{z}) d \mathbf{z}$:
$$
w_{k,i} \triangleq \int_{\mathcal{S}_{\mathrm{T}}} H_{i}(\mathbf{z}) J_{k}(\mathbf{z}) d \mathbf{z}
$$
​	Apparently, $A'_k(s), B'_k(s)$ for each user $k$ can be easily derived from the corresponding user $k$ block $\{\mu_k,\ \lambda_k\}_{k=1}^K$, and $w_{i,k}$ can also be calculated through following process:

1. multiply both sides by another new part $H_j(s)$ and integrate it over $ds$: 

$$
\begin{flalign}
\int_{\mathcal{S}_T}H_j(s)\times J_k(s)ds
&=\int_{\mathcal{S}_T}A'_kH_j(s)H_k^*(s)ds\ -\ \int_{\mathcal{S}_T}\sum_{i=1}^Kw_{k,i}\ B'_iH_j(s)H_i^*(s)ds\\
w_{k,j}&=A'_kh_{k,j}\ - \sum_{k=1}^K \ w_{k,i} B'_ih_{i,j}
\end{flalign}
$$

​	Having observed the form of above equation, we define a matrix $\mathbf{W}\in\mathbb{C}^{K\times K}$, where $\mathbf{W}(k,i)=w_{k,i}$, the above equation can be transferred into:
$$
\begin{flalign}
&\mathbf{W=HA-HBW}\\ 
\Rightarrow&\mathbf{(I}_K+\mathbf{HB)W=HA}\\
\Rightarrow&\mathbf{W=}\mathbf{(I}_K+\mathbf{HB)}^{-1}\mathbf{HA}
\end{flalign}
$$

​	Where the matrices respectively formed from following rules:
$$
\begin{cases}
\mathbf{A}=diag\{A'_1, ...A'_K\}\\
\mathbf{B}=diag\{B'_1,...B'_K\}\\
\mathbf{H}=[\mathbf{h}_1,..\mathbf{h}_K]\\
\mathbf{h}_k=[h_{k,1}\ ,...,\ h_{k,K}]^T
\end{cases}
$$

​	(ps: where the column and row is different from the normal way (actually is in a reverse way) that we think, so theoretically we can easily transfer it into a normal row and column way later.)

​	Having know the matrix of $\mathbf{W}$ (means you know every $w_{i,j}$) we can easily know the optimal $J_k(s)$ from $J_{k}(\mathbf{s})=A'_{k} H_{k}^{*}(\mathbf{s})-\sum_{i=1}^{K}w_{k,i}\ B'_{i} H_{i}^{*}(\mathbf{s})$ with known $A'_k$ and $B'_k$.

​	Considering too many steps contains integral part to heavy the cost of calculating, we also propose another similar approach without integral.

###### 		·Low-complexity Approach

​	Considering the relationship between $w_{k,i}$ and the block of $\{\mu_k,\ \lambda_k\}^K_{k=1}$, we easily find the elements of block is both the function of $w_{k,i}$ and reformulate the total power as:
$$
\begin{flalign}
\sum_{k=1}^K \int_{\mathcal{S}_T}|J_k(s)|^2ds
&=\sum_{k=1}^K|A'_k|^2 h_{k,k} - \sum_{k=1}^K \sum_{i=1}^K 2\Re{\{w_{k,i}^* \ A'_kB^{'*}_kh_{k,i}\}}+\sum_{k=1}^K\sum_{i=1}^K\sum_{j=1}^Kw_{k,i}\ w_{k,j}^*B'_i B^{'*}_jh_{i,j}\\
&=\operatorname{tr}\left(\mathbf{A}^{H} \mathbf{H A}-2 \Re\left\{\mathbf{A} \mathbf{W}^{H} \mathbf{B}^{H} \mathbf{H}\right\}+\mathbf{W}^{H} \mathbf{B}^{H} \mathbf{H B W}\right)\\
\end{flalign}
$$
​	Namely:
$$
\operatorname{tr}\left(\mathbf{A}^{H} \mathbf{H A}-2 \Re\left\{\mathbf{A} \mathbf{W}^{H} \mathbf{B}^{H} \mathbf{H}\right\}+\mathbf{W}^{H} \mathbf{B}^{H} \mathbf{H B W}\right)\triangleq \rho(\mathbf{W})
$$


​	Thus, apparently, the each elements of block $\{\mu_k, \ \lambda_k\}^K_{k=1}$ can be retell as:
$$
\begin{flalign}
\mu_{k}(\mathbf{W})&=\sqrt{1+\frac{\left|w_{k, k}\right|^{2}}{\sum_{i=1, i \neq k}^{K}\left|w_{i, k}\right|^{2}+\frac{\sigma_{k}^{2}}{P_{\mathbf{T}}} \rho(\mathbf{W})}}, \\
\lambda_{k}(\mathbf{W})&=\frac{\mu_{k}(\mathbf{W}) w_{k, k}}{\sum_{i=1}^{K}\left|w_{i, k}\right|^{2}+\frac{\sigma_{k}^{2}}{P_{\mathbf{T}}} \rho(\mathbf{W})}
\end{flalign}
$$
​	Thanks to this approach, we can hugely avoid too many calculation of integrals, except the initialization stage to calculate $h_{i,j}$ in matrix $\mathbf{H}$, only block $\{\mu_k, \ \lambda_k\}^K_{k=1}$ and $\mathbf{A,B,W}$ need to be iteratively updated.

###### ·CoV-BCD Approach/Algorithm Conclusion

​	Having roughly known the approach and low-complexity approach, we try to clarify each stages and the concrete sequence of it.

1. Calculate the channel correlation matrix $\mathbf{H}$ by $\mathbf{H}=[\mathbf{h}_1,..\mathbf{h}_K]$, where $\mathbf{h}_k=[h_{k,1}\ ,...,\ h_{k,K}]^T$ and $h_{k,j}=\int_{\mathcal{S}_T}H_j(s)H_k^*(s)ds$

2. Choose appropriate initial numbers of $J_k^0(s)$, $\mu_k^0$, $\lambda_k^0$

3. Repeat following steps:

   a. Update $\mathbf{A}^n$ and $\mathbf{B}^n$ with $\mu_k^n$ and $\lambda_k^n$ by:
   $$
   \mathbf{A}^n=diag\{A_1^n,A_2^n,...A_K^n\},\ \ \ A_k^n=\frac{(\alpha_k\mu_k^n\lambda_k^n)P_T}{\sum_{i=1}^K\alpha_i|\lambda_i^n|^2\sigma_i^2}\\
   \mathbf{B}^n=diag\{B_1^n,B_2^n,...B_K^n\},\ \ \ B_k^n=\frac{\alpha_k|\lambda_k^n|^2P_T}{\sum_{i=1}^K\alpha_i|\lambda_i^n|^2\sigma_i^2}
   $$

​	b. Update $\mathbf{W}^{n+1}$ with $\mathbf{A}^n$and $\mathbf{B}^n$ by $\mathbf{W=}\mathbf{(I}_K+\mathbf{HB)}^{-1}\mathbf{HA}$

​	c. Update $\mu_k^{n+1}$ and $\lambda_k^{n+1}$ with $\mathbf{W}^{n+1}$, $\mathbf{A}^n$, $\mathbf{B}^n$ by:
$$
\begin{flalign}
\mu_{k}(\mathbf{W})&=\sqrt{1+\frac{\left|w_{k, k}\right|^{2}}{\sum_{i=1, i \neq k}^{K}\left|w_{i, k}\right|^{2}+\frac{\sigma_{k}^{2}}{P_{\mathbf{T}}} \rho(\mathbf{W})}}, \\
\lambda_{k}(\mathbf{W})&=\frac{\mu_{k}(\mathbf{W}) w_{k, k}}{\sum_{i=1}^{K}\left|w_{i, k}\right|^{2}+\frac{\sigma_{k}^{2}}{P_{\mathbf{T}}} \rho(\mathbf{W})}\\
\end{flalign}
$$
​	Where $\rho(\mathbf{W})\triangleq \operatorname{tr}\left(\mathbf{A}^{H} \mathbf{H A}-2 \Re\left\{\mathbf{A} \mathbf{W}^{H} \mathbf{B}^{H} \mathbf{H}\right\}+\mathbf{W}^{H} \mathbf{B}^{H} \mathbf{H B W}\right)$.

4. Terminate this circle until the increase of objective value falls below a threshold.
5. Calculate $\bar{J}'_k(s)$ by $\bar{J}'_{k}(\mathbf{s})=A^n_{k} H_{k}^{*}(\mathbf{s})-\sum_{i=1}^{K}w_{k,i}\ B^n_{i} H_{i}^{*}(\mathbf{s})$, select $w_{k,i}$ from final $\mathbf{W}$.
6. Calculate $\bar{J}_k(s)$ by $\overline{J_k(s)}=\sqrt{\frac{P_T}{\sum_{k=1}^K\int_{\mathcal{S}_T}|\overline{J_k'(s)}|^2ds}}\overline{J_k'(s)}$ under power limit.

###### ·Corr-ZF Approach

​	In this part, inspired by conventional discrete MIMO systems, we introduce Corr-ZF approach to reduce the interference between users. It is Corr-ZF offers us an path to asymptotically approach the optimal capacity in multi-user systems and closed-form expression.

​	Based on the constraint of WSR and power, under ZF, we introduce another constraint (ZF constraint) to alleviate interference.
$$
\int_{\mathcal{S}_T}H_i(\mathbf{s})J_k(\mathbf{s})d\mathbf{s}=0, \forall i\neq k
$$

​	***Proposition 2***: the current distribution $\qty{J_k\pqty{\mathbf{s}}}_{k=1}^K$ satisfying the ZF constraint will be expressed as
$$
\begin{align}
J_k\pqty{\mathbf{s}}&=\sqrt{\rho_k}J_k^{\text{ZF}}(\mathbf{s})\\
J_k^\text{ZF}(\mathbf{s})&=\sum_{j=1}^Ku_{k,j}H_j^*(\mathbf{s})
\end{align}
$$
​	Where $\rho_k$ is power scaling factor, while $u_{k,j}$ is the entry of $\mathbf{H}^{-1}$ in $k$-th column and $j$-th row.

​	*Proof*: substitude above $J_k(\mathbf{s})$ into ZF constraint we can reorganize as follows
$$
\begin{align}
\int_{\mathcal{S}_T}H_i(\mathbf{s})J_k(\mathbf{s})d\mathbf{s}
&=\int_{\mathcal{S}_T}H_i(\mathbf{s})\sqrt{\rho_k}\sum_{j=1}^Ku_{k,j}H_j^*(\mathbf{s})d\mathbf{s}\\
&=\sqrt{\rho_k}\sum_{j=1}^Ku_{k,j}\int_{\mathcal{S}_T}H_i(\mathbf{s})H_j^*(\mathbf{s})d\mathbf{s}\\
&=\sqrt{\rho_k}\sum_{j=1}^Ku_{k,j}h_{i,j}\\
&=\sqrt{\rho_k}\mathbf{h}_i^T\mathbf{u}_k
\end{align}
$$
​	If we need to satisfy the ZF constraint, i.e., $\mathbf{h}_i^T\mathbf{u}_k=0,\forall i\neq k$, then it is not hard to derive $u_{k,j}$ is the entry of $\mathbf{H}^{-1}$, Proposition proved.

​	Based on this proposition, we further simplify the power scaling factor $\rho_k$ with $P_k$ as the power allocated to $k$-th user, satisfying $\sum_{k=1}^KP_k=P_T$.
$$
\rho_k=\frac{P_k}{\int_{\mathcal{S}_T}\vqty{J_k^\text{ZF}(\mathbf{s})}^2d\mathbf{s}}=\frac{P_k}{\mathbf{u}_k^H\mathbf{Hu}_k}=\frac{P_k}{u_{k,k}}
$$

​	Meanwhile, substitute the ZF form $J_k^\text{ZF}(\mathbf{s})=\sum_{j=1}^Ku_{k,j}H^*_j(\mathbf{s})$ to $\gamma_k$ can derive as
$$
\begin{align}
\gamma_k^\text{ZF}
&=\frac{\rho_k\vqty{\int_{\mathcal{S}_T}H_k(\mathbf{s})J_k^\text{ZF}(\mathbf{s})d\mathbf{s}}^2}{\rho_k\sum_{i=1,i\neq k}^K\vqty{\int_{\mathcal{S}_T}H_i(\mathbf{s})J_k^\text{ZF}(\mathbf{s})d\mathbf{s}}^2+\sigma_k^2}\\
&=\frac{P_k\vqty{\mathbf{h}_k^T\mathbf{u}_k}^2}{P_k\sum_{i=1,i\neq k}^K\vqty{\mathbf{h}_i^T\mathbf{u}_k}^2+u_{k,k}\sigma_k^2}\\
&=\frac{P_k}{u_{k,k}\sigma_k^2}
\end{align}
$$
​	Thereby, the problem is converted into a non-interference problem as follows
$$
\begin{align}
\max_{\qty{P_k}_{k=1}^K}&\sum_{k=1}^K\alpha_k\log_2\pqty{1+\frac{P_k}{u_{k,k}\sigma_k^2}}\\
\text{s.t.}&\sum_{k=1}^KP_k=P_T
\end{align}
$$
​	***Proposition 3***: The optimal solution to above problem is $P_k=\pqty{\nu\alpha_k-u_{k,k}\sigma_k^2}^+$, where
$$
\nu=\frac{P_T+\sum_{k=1}^Mu_{k,k}\sigma_k^2}{\sum_{k=1}^M\alpha_k}
$$
​	where $M$ is the number of non-zero $P_k$, and the above optimal solution is a classical water-filling solution, which is no need to prove.

###### ·Comparison with Fourier-Based Approach

​	In this section, we compare our scheme with the Fourier-based approach, which aim to approximate the continuous function of $H_k(\mathbf{s})$ and $J_k(\mathbf{s})$ for each user $k$. In the paper, the current pattern can be expressed as a Fourier series like $J_k(\mathbf{s})=\sum_{\mathbf{n}=-\infty}^{\infty}v_{k,n}\Phi_{\mathbf{n}}(\mathbf{s})$, where $\mathbf{n}$ means three dimension of $x,y,z$ and the $\sum$ is defined in the area of $x,y,z$. Meanwhile, the Fourier coefficients and basis functions are given by
$$
\begin{array}{ll}
&\text{coefficient:}& v_{k,n}=\frac{\int_{\mathcal{S}_T}J_k(\mathbf{s})\Phi_{\mathbf{n}}^*(\mathbf{s})d\mathbf{s}}{\sqrt{S_T}}\\
&\text{basis functions:}& \Phi_\mathbf{n}=\frac{e^{j2\pi\pqty{\frac{n_x}{L_x}\pqty{s_x-\frac{L_x}{2}}+\frac{n_y}{L_y}\pqty{s_y-\frac{L_y}{2}}+\frac{n_z}{L_z}\pqty{s_z-\frac{L_z}{2}}}}}{\sqrt{S_T}}
\end{array}
$$
​	where $\mathbf{s}=[s_x,s_y,s_z]^T$ and $L_x$, $L_y$ and $L_z$ are the projection lengths of $S_T$ on $x,y,z$, respectively. Then the finite-dimension approximation are typically as follows
$$
J_k(\mathbf{s})\approx\sum_{\mathbf{n=-N}}^\mathbf{N}v_{k,n}\Phi_\mathbf{n}(\mathbf{s})
$$
​	Where the $\mathbf{n}$, $\mathbf{N}$ and $\sum$ are similarly defined and $N_x,N_y,N_z$ is with
$$
\begin{array}{ccc}
N_x=\lceil\frac{L_x}{\lambda}\rceil&N_y=\lceil\frac{L_y}{\lambda}\rceil&N_z=\lceil\frac{L_z}{\lambda}\rceil
\end{array}
$$
​	Then the approximation can be substitute into transmit power and received power as
$$
\begin{align}
\int_{\mathcal{S}_T}\vqty{J_k(\mathbf{s})}^2d\mathbf{s}&\approx\sum_{\mathbf{n=-N}}^\mathbf{N}\vqty{v_{k,\mathbf{n}}}^2=\norm{\mathbf{v}_k}^2\\
\int_{\mathcal{S}_T}H_i(\mathbf{s})J_k(\mathbf{s})d\mathbf{s}&\approx\sum_{\mathbf{n=-N}}^\mathbf{N}g_{i,\mathbf{n}}v_{k,\mathbf{n}}=\mathbf{g}_i^T\mathbf{v}_k
\end{align}
$$
​	Where $g_{i,\mathbf{n}}=\int_{\mathcal{S}_T}H_i(\mathbf{s})\Phi_\mathbf{n}(\mathbf{s})d\mathbf{s}$ is the Fourier Transform of $H_i(\mathbf{s})$, and $\mathbf{v}_k,\mathbf{g}_i\in\mathbb{C}^{N_\mathbf{F}\times 1}$ with $N_\mathbf{F}=\pqty{2N_x+1}\pqty{2N_y+1}\pqty{2N_z+1}$. Naturally, the WSR problem is approximated as
$$
\begin{align}
\max_{\qty{\mathbf{v}_k}_{k=1}^K}&\sum_{k=1}^K\alpha_k\log_2\pqty{1+\frac{\vqty{\mathbf{g}_k^T\mathbf{v}_k}^2}{\sum_{i=1,i\neq k}^K\vqty{\mathbf{g}_k^T\mathbf{v}_i}^2+\sigma_k^2}}\\
\text{s.t.}&\sum_{k=1}^K\norm{\mathbf{v}_k}^2=P_T
\end{align}
$$
​	Apparently, the above problem is a classical WSR problem in spatially-discrete arrays, normally solved by FP or WMMSE. However, it is the Fourier-Transform that makes the computational complexity hugely improved with a loss of approximation. The detail analysis will be as follows





###### ·Numerical Results

​	In this section, numerical results are provided to verify the effectiveness of our proposed CoV-based and Corr-ZF scheme for CAPA beamforming by Monte-Carlo simulations. 

​	First, the setup is exploited as following unless notice specified. We assume that CAPA is deployed within a $x-y$ plane as
$$
\mathcal{S}_T=\qty{[s_x,s_y,0]^T\vert\vqty{s_x}\leq\frac{L_x}{2},\vqty{s_y}\leq\frac{L_y}{2}}
$$
​	and $L_x=L_y=\sqrt{S_T}$ and $S_T=0.25\text{m}^2$ and $K=8$ communication users randomly located within the region as following
$$
\mathcal{U}=\qty{[r_x,r_y,r_z]^T\vert\vqty{r_x}\leq U_x,\vqty{r_y}\leq U_y,U_{z,\min}\leq r_z \leq U_{z,\max}}
$$
​	Where $U_x=U_y=5\ \text{m}$ and $U_{z,\min}=15\ \text{m}$, $U_{z,\max}=30\ \text{m}$ and the polarization direction of BS are set to $\mathbf{\hat{u}}_y=[0,1,0]^T$, while users is also $\mathbf{\hat{u}}_k=[0,1,0]^T$. The signal frequency is set to $f=2.4\text{GHz}$ and intrinsic impedance is set to $\eta=120\pi\Omega$. The transmit power is to $P_T=100\text{mA}^2$ and $\sigma_k^2=5.6\times 10^{-3}\text{V}^2/\text{m}^2$, and the weight of each user is set to $\alpha_k=1/K$, temporarily. The Gauss-Legendre point is $M=20$, and all results are obtained by averaging over $200$ random channel realizations unless specified.

​	Secondly, to compare, we mainly consider the following benchmark schemes

1. ***Fourier-based Approach***. This approach is already described before
2. ***Conventional MIMO***. The conventional spatially discrete antenna array is exploited to compare, where the continuous surface $\mathcal{S}_T$ is replaced with discrete antennas with $A_d=\frac{\lambda^2}{4\pi}$ effective aperture area of antenna and $d=\frac{\lambda}{2}$. The location of the $\pqty{n_x,n_y}$-th discrete antenna is given by

$$
\mathbf{\bar{s}}_{n_x,n_y}=\bqty{\pqty{n_x-1}d-\frac{L_x}{2},\pqty{n_y-1}d-\frac{L_y}{2}}^T
$$

so there will be totally $N_d=\lceil\frac{L_x}{d}\rceil\times\lceil\frac{L_y}{d}\rceil$ discrete antennas, and the channel between the $\pqty{n_x,n_y}$-th discrete antenna and user $k$-th can be derived as
$$
h_{n_x,n_y}=\sqrt{A_d}\hat{\mathbf{u}}_k^T\mathbf{G}\pqty{\mathbf{r}_k,\bar{\mathbf{s}_{n_x,n_y}}}\hat{\mathbf{u}}_y
$$
Furthermore, the conventional MIMO beamforming optimization problem for maximize WSR will be solved by using traditional methods, like FP and WMMSE methods, and also include ZF to reduce the beamforming complexity.

***Convergence and Complexity of CoV Approach***

​	Fig 3 illustrates the advantages on the performance between our BCD-CoV approach and the Fourier-based scheme, the Fourier-based schemes will hugely increase the calculation resource when the improvement of aperture size and frequency, which is owing to the highly require on the reserved Fourier series items $N_\text{F}$ when improving the size of aperture and the frequency. However, our BCD-CoV approach always performance well in any condition, no matter which size and frequency are equipped with the BS.

​	By the way, our scheme is also satisfied with the command on convergence.










##### · *Multi-User Continuous-Aperture Array Communications:  How to Learn Current Distribution?*

###### ·Abstract

​	Aiming to cope with a non-convex functional optimization problem without closed-form objective function and constraint, this paper propose a L-CAPA deep learning structure to optimize the distribution strategy. In such a DL structure, we find a finite-dimensional variable of channel functions and current distributions as the input of neural network, and two more DNN specifically for integrals in loss function without closed-form expression to better address training.

###### ·System Model

​	We define $\mathcal{A}$ and $\mathcal{A}_k\subseteq\mathbb{R}^{3\times 1}$ as CAPAs at BS and $k$-th user in 3D space respectively, while $k$-th user centered at $\mathbf{s}_k\in\mathbb{R}^{3\times 1}$. Subsequently, we can express the transmitted signal at point $\mathbf{r}\in\mathcal{A}$ at BS as:
$$
x(\mathbf{r})=\sum_{k=1}^K\mathrm{V}_k(\mathbf{r})s_k
$$
​	Where $s_k$ is the symbol to be transmitted to the $k$-th user with $\mathbb{E}\{|s_k|^2\}=1$, namely orthogonal and $\mathrm{V}_k(\mathbf{r})$ is the continuous current distribution to convey $s_k$, thus, the SINR can be easily expressed as:
$$
\gamma_{k} \approx \frac{\left|\mathcal{A}_{k}\right| \cdot\left|\int_{\mathcal{A}} \mathrm{H}_{k}(\mathbf{r}) \mathrm{V}_{k}(\mathbf{r}) d \mathbf{r}\right|^{2}}{\sum_{j=1, j \neq k}^{K}\left|\mathcal{A}_{j}\right| \cdot\left|\int_{\mathcal{A}} \mathrm{H}_{k}(\mathbf{r}) \mathrm{V}_{j}(\mathbf{r}) d \mathbf{r}\right|^{2}+\sigma_{k}^{2}}
$$
​	Where $\mathrm{H}_k(r)$ represents the channel response from $\mathbf{r}$ to the center of $k$-th user, $|\mathcal{A}_k|$ is the aperture size of $k$-th user (smaller enough compared to BS and distance) and $\sigma_k^2$ means noise. Considering the LoS channels, $\mathrm{H}_k(\mathbf{r})$ can be expressed as:

$$
\mathrm{H}_k(\mathbf{r})=\sqrt{\frac{\mathbf{e}_r^T(\mathbf{s}_k-\mathbf{r})}{\norm{\mathbf{r-s}_k}}}\cdot\frac{jk_0\eta e^{-jk_0\norm{\mathbf{r-s}_k}}}{4\pi\norm{\mathbf{r-s}_k}}(1+\frac{j/k_0}{\norm{\mathbf{r-s}_k}}-\frac{1/k_0^2}{\norm{\mathbf{r-s}_k}^2})
$$

​	Where $\mathbf{e}_r\in\mathbb{R}^{3\times 1}$ is normal vector of BS surface, $\eta=120\pi$ is impedance of free space and $k_0=\frac{2\pi}{\lambda}$ is wavenumber with wavelength $\lambda$.

​	Thus, the current distribution optimization problem can be reformulated and policy can be clearly represented. The optimization can be regard as the vector mapping from channel space vector $\{\mathbf{H}_k(\cdot)\}_{k=1}^K$ to the distribution vector $\{\mathbf{V}_k(\cdot)\}_{k=1}^K$.

​	ps: We use $K$ set and $\{\mathbf{V}_k(\cdot),\mathbf{H}_k(\cdot)\}$ due to reducing misunderstanding of mapping, actually, the distribution $\mathbf{V}_k(\cdot)$ of any $k$-th user and at any point $\mathbf{r}$ on BS is influenced by all $\mathbf{H}$.

###### ·Proposed Deep-Learning Framework

​	Now, we already know the mapping from $\mathbf{H}_k(\cdot)$ to $\mathbf{V}_k(\cdot)$, but we are still facing two problem

​		·Issue 1: Both $\mathbf{H}$ and $\mathbf{V}$ are continuous (means infinite dimension) can not be input of DNN

​		·Issue 2: The integrals of objective function hinder the training.

​	Thus, we deal with these problems respectively.

​	Firstly, we observe the expression of $\mathbf{H}_k(\mathbf{r})$, the parameter $\mathbf{e}_r^T,\eta,k_0$ are fixed, and $\mathbf{s}_k$ can also be fixed after knowing the concrete location of each $k$-th user, only point $\mathbf{r}$ on BS as the variable of $\mathbf{H}_k(\mathbf{r})$, which means every $\mathbf{H}$ can be uniquely represented by $\mathbf{s}_k$ and as the finite vector input for neural network.

​	Similarly, having proved that $\mathbf{V}_k(\mathbf{r})=\sum_{j=1}^Ka_{jk}\mathbf{H}_j(\mathbf{r})$, naturally $\mathbf{V}$ can also be represented by $a_{jk}$ due to the fixed $\mathbf{H}_j(\mathbf{r})$ and as the output for neural network.

​	Having found the substitude of $\mathbf{H}_k(\mathbf{r}),\mathbf{V}_k(\mathbf{r})$ to finish mapping, we can learn the distribution policy by learning $\mathbf{A}\triangleq [\mathbf{a}_1,...,\mathbf{a}_K]\in\mathbb{C}^{K\times K}$, where $\mathbf{a}_k\triangleq[a_{1k},...,a_{Nk}]^T$ (which can be derived into $\mathbf{V}$, $\mathbf{V=HA}$, where $\mathbf{V,H}\in \mathbb{C}^{1\times K}$) from $\mathbf{S}\triangleq[\mathbf{s}_1,...,\mathbf{s}_K]^T\in \mathbb{C}^{K\times 3}$, which is oringinated from $\mathbf{H}$. Eventually, we learn the optimal $\mathbf{A}^*$ corresponding to optimal distribution.

##### 	· *Holographic Integrated Data and Energy  Transfer (Conference)*

###### ·Abstract

​	Meta-materials bring a higher spatial diversity gain by generate any current distribution on surface, it is the aid of EM manipulation of H-MIMO, that H-IDET is yielded, who can realize energy focusing and inter-use interference elimination by fully exploitation of the EM channels.

​	This paper proposed a BCD based scheme to solve the non-convex functional programing in optimizing H-IDET system（H-IDET, where maximize the sum-rate of data users by guaranteeing the energy of energy users）, together with considering the FT and equivalence between SNIR and MSE, and enhance the robustness by SCA and an initialization scheme.

​	The numerical result is better than benchmark, especially the one adopting discrete antennas, while the near-field focusing using EM channel also better than traditional ones, especially for WPT.

#### · PASS

##### · Exploiting Pinching-Antenna Systems  in Multicast Communications

###### · Abstract

​	The PASS is a novel antenna structure for reconfiguring wireless by pinching beamforming, where the activated locations of PA are optimized. This paper investigate application of PASS in multicast to maximize the multicast rate. ***In the single-waveguide scenario*** with single PA and linearly user location, a closed-form optimal activated location is derived. Then, a closed-form expression for achievable multicast rate is obtained and proven to outperform the fixed antenna. For the general multiple-PA system with arbitrary user distributions, an element-wise AO-based algorithm is proposed to beamform. ***In the multiple-waveguide scenario***, an AO-based scheme is design to jointly optimize the transmit and pinching beamformer. The transmit beamformer is updated by MM framework and SOCP, while the pinching beamformer is optimized by element-wise sequential refinement. Numerical results demonstrate that PASS can achieve higher multicast rate, especially when larger users and spatial coverage, while more PAs can further improves the performance of system.

###### · System Model

​	In this section, we propose a PASS-enabled multicast system, the BS is equipped with $M$ waveguides, each equipped with $N$ PAs to serve $K$ users. The location of users $K\in\mathcal{K}$ is given by $\psi_k=\bqty{\hat{x}_k,\hat{y}_k,0}^\mathrm{T}$, and distributed in a square region of $D=D_x\times D_y\text{m}^2$ randomly. 

​	Assume that all waveguides are with fixed height of $h$, and the PASS spans across region $D$ and $\mathcal{M}$, $\mathcal{N}_m$ denote the sets of waveguides and PAs on $m$-th waveguide, respectively. The waveguides are spaced along the $y$-axis with equal interval of $d_y=D_y/\pqty{M-1}\text{m}$

#### · RSMA

##### · An Efficient Max-Min Fair Resource Optimization Algorithm for Rate-Splitting Multiple Access

###### · Abstract

​	The max-min fairness (MMF) problem in rate-splitting multiple access (RSMA) is challenging for its non-convex and non-smooth nature, with coupled beamforming and common rate variables. In this work, we propose a optimization algorithm named extragradient-fractional programming (EG-FP) to address the MMF problem of downlink RSMA without high computational complexity and MMF rate performance degradation. The proposed framework first use FP to convert the original problem into a block-wise convex problem, whose Lagrangian dual is equal to a variational inequality problem and solved by EG, as a subproblem. Based on the above, we discover the optimal beamforming structure and a low-dimensional EG-FP with a transmit antenna computational complexity is introduced then, which is benefit to scenarios with a large scale of antennas, and we also extend this to the imperfect CSI condition. 

#### ·Multicast

##### · *Intelligent Reflecting Surface Aided Multigroup  Multicast MISO Communication Systems*

​	 This paper considered a multicast transmission system (with multiple antenna BS and single antenna users) by introducing RIS to maximize the minimum SINR of each multigroup.

###### ·System Model

​	we can naturally know the signal we sent is:
$$
x=\sum_{g=1}^Gx_g=\sum_{g=1}^G \pmb{f}_gs_g
$$
​	where $G$ means the set of multigroup, $s_g$ means the data that group $g$ want and follows $\mathbb{E}\{|s_g|^2=1\}$, $\pmb{f}_g$ means corresponding precode matrix, which still be limited by power of BS.

​	Define the diagonal matrix $\pmb{E}=diag([e_1,...,e_M]^T)\in \mathbb{C}^{M\times M}$ as the reflection coefficient matrix, and each $e_i$ satisfy $|e_i|^2=1$.

​	Due to the BS with N antenna, while RIS with M antenna, we can simply define the channel from BS to user $k$ as $\mathbf{h}_{d,k}\in\mathbb{C}^{N\times 1}$, from BS to RIS as $\mathbf{H}_{d,r}\in\mathbb{C}^{M\times N}$, (while the reflection coefficient matrix is $\mathbf{E}=\text{diag}([e_1,...e_M]^T)\in\mathbb{C}^{M\times M}$) from RIS to user $k$ as $\mathbf{h}_{r,k}\in\mathbb{C}^{M\times 1}$, respectively.

​	Assumed BS knows CSI and send to RIS, the signal received by each user of  group $g$ is:
$$
y_{k}=\left(\mathbf{h}_{\mathrm{r}, k}^{\mathrm{H}} \mathbf{E} \mathbf{H}_{\mathrm{dr}}+\mathbf{h}_{\mathrm{d}, k}^{\mathrm{H}}\right) \sum_{g=1}^{G} \mathbf{f}_{g} s_{g}+n_{k}
$$

​	Where $\mathbf{h}_{d,h}^H$ represents the direct channel, while $\mathbf{h}_{r,k}^{\mathrm{H}}\mathbf{EH}_{d,r}$ represents the indirect channel through RIS.

​	*Note*: We already know the $\mathbf{h}_{d,k}^\mathrm{H}\in\mathbb{C}^{1\times N}$, so $diag(\mathbf{h}_{d,k}^\mathrm{H})\in\mathbb{C}^{M\times M}$, while $\mathbf{E}\in\mathbb{C}^{M\times M}$ and only have diagonal elements will easily become to $\mathbf{e}'=[e_1,...e_M]^T\in\mathbb{C}^{1\times M}$.

​	We can transfer and gain the $\mathbf{H}_k=[diag(\mathbf{h}_{r,k}^\mathrm{H})\mathbf{H}_{dr}, \mathbf{h}_{d,k}^\mathrm{H}]^T$, and $\mathbf{e}=[e_1,...,e_M,1]^T$, so:
$$
y_k=\mathbf{e}^\mathrm{H}\mathbf{H}_k\pmb{x}+n_k=\mathbf{e}^\mathrm{H}\mathbf{H}_k\sum_{g=1}^G\mathbf{f}_gs_g+n_k
$$

$$
\gamma_k=\frac{|(\mathbf{h}_{r,k}^\mathrm{H}\mathbf{EH}_{dr}+\mathbf{h}_{d,k}^\mathrm{H})\mathbf{f}_g|^2}{\sum_{i=1,i\neq g}^G|(\mathbf{h}_{r,k}^\mathrm{H}\mathbf{EH}_{dr}+\mathbf{h}_{d,k}^\mathrm{H})\mathbf{f}_i|^2+\sigma^2}=\frac{|\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g|^2}{\sum_{i=1,i\neq g}^G|\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_i|^2+\sigma^2}
$$

​	Thus, having known the SINR, we easily calculate the ratio as following:
$$
\begin{flalign}
\begin{array}{ll}
&R_k=\log_2(1+\gamma_k), \ \ \ &\forall k\in \mathcal{K}\\
&R_g=\min_{k\in \mathcal{K}_g}\{R_k(\mathbf{F,e})\}, &\forall g\in\mathcal{G}
\end{array}
\end{flalign}
$$
​	Because the ratio of any group is always limited by the worst channel (the minimum SINR).Thus, we can also derive the practical ratio of each group. Where $R_k$ is the theoretical ratio of each user $k$ in their corresponding group, while $R_g$​ is the practical ratio of each user in group $g$. Formulate the problem as:
$$
\max_{\{\mathbf{F,e}\}}\{F(\mathbf{F,e}) \}=\max_{\{\mathbf{F,e}\}}\{\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\{R_k(\mathbf{F,e})\}\}
$$
​	Then, we try to optimize it by finding an easier surrogate problem, under MM framework, referring to  ***Majorization-Minimization Algorithms in Signal Processing, Communications, and Machine Learning*** and ***A Tutorial on MM Algorithms*** or [非凸优化问题的大杀器：Majorization-Minimization 算法 - 知乎](https://zhuanlan.zhihu.com/p/642700446)

###### ·SOCP MM Method

​	This section, we propose a SCOP-based MM method. Under this framework, we transfer the originally complex non-convex objective function into a concave surrogate one and alternately optimize $\mathbf{F}$ and $\mathbf{e}$ with AO method.

​	Considering that the composition function $F(\mathbf{F,e})$ is a  linear combination of minimum subfunction $R_k(\mathbf{F,e})$ in same group like $k\in\mathcal{K}_g$, we tackle this non-concave function first.

​	If $\{\mathbf{F}^n,\mathbf{e}^n\}$ is the $(n-1)$-th solution obtained by iteration, then the $R_k(\mathbf{F,e})$ can be surrogated by $\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)$ as MM method and defined by
$$
\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)=\text{const}_k+2\Re\{a_k\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g\}-b_k\sum_{i=1}^G|\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_i|^2\leq R_k(\mathbf{F,e})
$$
​	Where we define some parameter under fixed point $\{\mathbf{F}^n,\mathbf{e}^n\}$ as following
$$
\begin{align}
a_{k}&=\frac{\left(\mathbf{f}_{g}^{n}\right)^{\mathrm{H}} \mathbf{H}_{k}^{\mathrm{H}} \mathbf{e}^{n}}{\sum_{i \neq g}^{G}\left|\left(\mathbf{e}^{n}\right)^{\mathrm{H}} \mathbf{H}_{k} \mathbf{f}_{i}^{n}\right|^{2}+\sigma_{k}^{2}}\\
b_{k}&=\frac{\left|\left(\mathbf{e}^{n}\right)^{\mathrm{H}} \mathbf{H}_{k} \mathbf{f}_{g}^{n}\right|^{2}}{\left(\sum_{i \neq g}^{G}\left|\left(\mathbf{e}^{n}\right)^{\mathrm{H}} \mathbf{H}_{k} \mathbf{f}_{i}^{n}\right|^{2}+\sigma_{k}^{2}\right)\left(\sum_{i=1}^{G}\left|\left(\mathbf{e}^{n}\right)^{\mathrm{H}} \mathbf{H}_{k} \mathbf{f}_{i}^{n}\right|^{2}+\sigma_{k}^{2}\right)} \\
\text {const}_{k}&=R_{k}\left(\mathbf{F}^{n}, \mathbf{e}^{n}\right)-b_{k} \sigma_{k}^{2}-b_{k}\left(\sum_{i=1}^{G}\left|\left(\mathbf{e}^{n}\right)^{\mathrm{H}} \mathbf{H}_{k} \mathbf{f}_{i}^{n}\right|^{2}+\sigma_{k}^{2}\right)
\end{align}
$$
​	Tackled above equation, we have already surrogated the original objective optimization function $R_k(\mathbf{F,e})$ by $\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)$ (which is concave and easier to solve the maximum point), the optimization problem is transformed into
$$
\begin{align}
\max_{\{\mathbf{F,e}\}}\{\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\{R_k(\mathbf{F,e})\}\}&\Rightarrow\max_{\{\mathbf{F,e}\}}\{\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\{\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)\}\}\\
\max_{\{\mathbf{F,e}\}}\{\tilde{F}(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n) \}&=\max_{\{\mathbf{F,e}\}}\{\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\{\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)\}\}
\end{align}
$$
​	Then, we commence on optimizing this bio-concave function, in such a condition, to clarify our problem we try to solve the two concave functions respectively. Concretely, we consider optimize one of them under another variable fixed.

**·Matrix $\mathbf{F}$**

​	With given fixed $\mathbf{e}$, we try to optimize the precoding matrix $\mathbf{F}$, the $\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)$ can be shown to be a quadratic function of $\mathbf{F}$, i.e., $\tilde{R}_k(\mathbf{F|F}^n)$.
$$
\begin{align}
\tilde{R}_k(\mathbf{F|F}^n)
&=\text{const}_k+2\Re\{a_k\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g\}-b_k\sum_{i=1}^G|\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_i|^2\\
&=\text{const}_k+2\Re\{\Tr{[\mathbf{C}_k^\mathrm{H}\mathbf{F}]}\}-\Tr{[\mathbf{F}^\mathrm{H}\mathbf{B}_k\mathbf{F}]}
\end{align}
$$
​	Where $\mathbf{B}_k=b_k\mathbf{H}_k^\mathrm{H}\mathbf{ee}^\mathrm{H}\mathbf{H}_k$, $\mathbf{C}_k^\mathrm{H}=a_k\mathbf{t}_g\mathbf{e}^\mathrm{H}\mathbf{H}_k$, and where $\mathbf{t}_g\in\R^{G\times 1}$ is a specified selection vector to choose $g$-th element/group. Then the problem thus reorganized to
$$
\max_\mathbf{F}\{\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\{\text{const}_k+2\Re\{\Tr{[\mathbf{C}_k^\mathrm{H}\mathbf{F}]}\}-\Tr{[\mathbf{F}^\mathrm{H}\mathbf{B}_k\mathbf{F}]}\}\}
$$
​	To better tackle with this objective function, we introduce auxiliary variables $\gamma=[\gamma_1,\cdots,\gamma_G]^T$ and the problem as follows
$$
\begin{align}
\max_{\mathbf{F},\gamma}\ \ &\sum_{g=1}^G\gamma_g\\
\text{s.t.}\ \ &\mathbf{F}\in\mathcal{S}_F\\
&\text{const}_k+2\Re\{\Tr{[\mathbf{C}_k^\mathrm{H}\mathbf{F}]}\}-\Tr{[\mathbf{F}^\mathrm{H}\mathbf{B}_k\mathbf{F}]}\geq\gamma_g\\
&\forall k\in\mathcal{K}_g,\ \ \forall g\in\mathcal{G}
\end{align}
$$
​	Which is a SOCP problem can be easily solved by the CVX solver, like MOSEK.

**·Vector $\mathbf{e}$**

​	With given fixed $\mathbf{F}$, we try to optimize the reflection coefficient vector $\mathbf{e}$, and the $\tilde{R}_k(\mathbf{F,e}|\mathbf{F}^n,\mathbf{e}^n)$ can be rewritten as $\tilde{R}_k(\mathbf{e|e}^n)$, i.e.
$$
\tilde{R}_k(\mathbf{e|e}^n)=\text{const}_k+2\Re\{a_k\mathbf{H}_k\mathbf{f}_g\mathbf{e}\}-\mathbf{e}^\mathrm{H}\mathbf{A}_k\mathbf{e}
$$
​	Where $\mathbf{A}_k=b_k\mathbf{H}_k\sum_{i=1}^G\mathbf{f}_i\mathbf{f}_i^\mathrm{H}\mathbf{H}_k^\mathrm{H}$. Thus, replace the optimization problem by above equation, we gain
$$
\max_\mathbf{e}\{\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\{\text{const}_k+2\Re\{a_k\mathbf{H}_k\mathbf{f}_g\mathbf{e}\}-\mathbf{e}^\mathrm{H}\mathbf{A}_k\mathbf{e}\}\}
$$
​	Introduce auxiliary variables $\mathbf{\kappa}=[\kappa_1,\cdots,\kappa_G]^T$ and the problem as follows
$$
\begin{aligned}
\max _{\mathbf{e}, \boldsymbol{\kappa}} & \sum_{g=1}^{G} \kappa_{g} \\
\text {s.t.} & \mathbf{e} \in \mathcal{S}_{e} \\
& \text {const}_{k}+2 \operatorname{Re}\left\{\mathbf{a}_{k}^{\mathrm{H}} \mathbf{e}\right\}-\mathbf{e}^{\mathrm{H}} \mathbf{A}_{k} \mathbf{e} \geq \kappa_{g}, \\
& \forall k \in \mathcal{K}_{g}, \forall g \in \mathcal{G}.
\end{aligned}
$$
​	Which is still non-concave due to non-concave variable $\mathbf{e}\in\mathcal{S}_e$, thus, we replace it with a relaxed convex one $\mathbf{e'}=\{e^2_1,\cdots,e^2_M,1\},\forall i\in[1,M], e_i<1$, and $\mathcal{S}_{e-relax}$ is  the set of $\mathbf{e'}$. 
$$
\begin{aligned}
\hat{\mathbf{e}}=\arg\max _{\mathbf{e}} & \sum_{g=1}^{G} \gamma_{g} \\
\text {s.t.} & \mathbf{e} \in \mathcal{S}_{e-relax} \\
& \text {const}_{k}+2 \operatorname{Re}\left\{\mathbf{a}_{k}^{\mathrm{H}} \mathbf{e}\right\}-\mathbf{e}^{\mathrm{H}} \mathbf{A}_{k} \mathbf{e} \geq \kappa_{g}, \\
& \forall k \in \mathcal{K}_{g}, \forall g \in \mathcal{G}.
\end{aligned}
$$
​	Where $\hat{\mathbf{e}}$ means the optimal solution of the above relaxed version of problem, and the originally optimal solution $\mathbf{e}$ in the $n$-th will be derived by
$$
\mathbf{e}^{n+1}=
\begin{cases}
\exp\{j\angle\left(\frac{\hat{\mathbf{e}}}{\hat{\mathbf{e}}_{M+1}}\right)\},\text{if}\ F\left(\mathbf{F}^{n+1},\exp\{j\angle\left(\frac{\hat{\mathbf{e}}}{\hat{\mathbf{e}}_{M+1}}\right)\}|\mathbf{F}^n,\mathbf{e}^n\right)\geq F\left(\mathbf{F}^{n+1},\mathbf{e}^n|\mathbf{F}^n,\mathbf{e}^n\right)\\
\mathbf{e}^n, \text{otherwise}
\end{cases}
$$
​	Where $\hat{\mathbf{e}}_{M+1}$ means the $m$-th element of $\hat{\mathbf{e}}$ and all operations is made based on elements.

###### ·Low-Complexity MM Method

​	As seen in above, it is too complex and waste to tackle two SOCP problems in each iteration, thus, we commence on deriving a closed-form and low-complexity algorithm subsequently.

​	Since $\min _{k\in\mathcal{K}_g}\{\tilde{R}_k(\mathbf{F,e|F}^n,\mathbf{e}^n)\}$ is non-differentiable, we approximate it as a smooth/differentiable function by using the following smooth log-sum-exp lower bound
$$
\min_{k\in\mathcal{K}_g}\{\tilde{R}_k(\mathbf{F,e|F}^n,\mathbf{e}^n)\}\approx f_g(\mathbf{F,e})\triangleq-\frac{1}{\mu_g}\log_2\left(\sum_{k\in\mathcal{K}_g}\exp\{-\mu_g\tilde{R}_k(\mathbf{F,e|F}^n,\mathbf{e}^n)\}\right)
$$
​	Where $f_g$ surrogate the $\min\{\}$ function, directly representing the rate of each multi-group $g$, and $\mu_g>0$ is a smooth parameter satisfying
$$
f_g(\mathbf{F,e})\leq\min_{k\in\mathcal{K}_g}\{\tilde{R}_k(\mathbf{F,e|F}^n,\mathbf{e}^n)\}\leq f_g(\mathbf{F,e})+\frac{1}{\mu_g}\log_2\left(|\mathcal{K}_g|\right)
$$
​	Proved that $f_g(\mathbf{F,e})$ is bio-concave of $\mathbf{F}$ ad $\mathbf{e}$ by [29] for its Hessian matrix is semi-negative.

​	Larger $\mu_g$, higher approximation accuracy, but causing problem nearly ill-condition. Only choose $\mu_g$ appropriately, the optimization problem can be approximated as
$$
\max_{\mathbf{F,e}}\ \ \sum_{g=1}^Gf_g(\mathbf{F,e})\\
\text{s.t.}\ \ \mathbf{F}\in\mathcal{S}_F,\mathbf{e}\in\mathcal{S}_e
$$
​	Which is still a bio-convex problem of $\mathbf{F,e}$, we can update $\mathbf{F,e}$ alternately by AO method.

**·Optimize matrix $\mathbf{F}$**

​	Similarly, given fixed $\mathbf{e}$ the subproblem of $\mathbf{F}$ optimization is
$$
\max_{\mathbf{F}} \ \ \sum_{g=1}^Gf_g(\mathbf{F})\\
\text{s.t.}\ \ \ \mathbf{F}\in\mathcal{S}_F
$$
​	Which is concave and continuous, but also hard, we give another surrogate function of $f_g(\mathbf{F})$ in MM as
$$
\tilde{f}_g(\mathbf{F|F}^n)=2\Re\{\Tr[\mathbf{U}_g^\mathrm{H}\mathbf{F}]\}+\alpha_g\Tr[\mathbf{F}^\mathrm{H}\mathbf{F}]+\text{consF}_g
$$
​	Where
$$
\begin{array}{l}
\mathbf{U}_g=\sum_{k\in\mathcal{K}_g} g_k\left(\mathbf{F}^{n}\right)\left(\mathbf{C}_{k}-\mathbf{B}_{k}^{\mathrm{H}} \mathbf{F}^{n}\right)-\alpha_{g} \mathbf{F}^{n}, \\
g_{k}\left(\mathbf{F}^{n}\right)=\frac{\exp \left\{-\mu_{g} \widetilde{R}_{k}\left(\mathbf{F}^{n}\right)\right\}}{\sum_{k \in \mathcal{K}_{g}} \exp \left\{-\mu_{g} \widetilde{R}_{k}\left(\mathbf{F}^{n}\right)\right\}}, k \in \mathcal{K}_{g}, \\
\alpha_{g}=-\max_{k \in \mathcal{K}_g}\left\{b_{k} \mathbf{e}^{\mathrm{H}} \mathbf{H}_{k} \mathbf{H}_{k}^{\mathrm{H}} \mathbf{e}\right\}-2 \mu_g\max_{k \in \mathcal{K}_{g}}\left\{t p_{k}\right\}, \\
t p_{k}=P_{\mathrm{T}} b_{k}^{2}\left|\mathbf{e}^{\mathrm{H}} \mathbf{H}_{k} \mathbf{H}_{k}^{\mathrm{H}}\right|^{2}+\left\|\mathbf{C}_{k}\right\|_{F}^{2}+2 \sqrt{P_{\mathrm{T}}}\left\|\mathbf{B}_{k} \mathbf{C}_{k}\right\|_{F}, \\
\operatorname{consF}_{g}=f_{g}\left(\mathbf{F}^{n}\right)+\alpha_{g} \operatorname{Tr}\left[\left(\mathbf{F}^{n}\right)^{\mathrm{H}} \mathbf{F}^{n}\right]-2 \operatorname{Re}\left\{\operatorname{Tr}\left[\mathbf{D}_{g}^{\mathrm{H}} \mathbf{F}^{n}\right]\right\} .
\end{array}
$$
​	With above definition, the problem transferred into
$$
\max_\mathbf{F}\sum_{g=1}^G(2\Re\{\Tr[\mathbf{U}_g^\mathrm{H}\mathbf{F}]\}+\alpha_g\Tr[\mathbf{F}^\mathrm{H}\mathbf{F}]+\text{consF}_g)\\
\text{s.t.}\ \ \mathbf{F}\in\mathcal{S}_F
$$
​	The optimal $\mathbf{F}^{n+1}$ gained by introducing a Lagrange multiplier $\tau\geq0$ associated with the power constraint, then the function as
$$
\mathcal{L}(\mathbf{F},\tau)=2\Re\left\{\Tr\left[\sum_{g=1}^G\mathbf{U}_g^\mathrm{H}\mathbf{F}\right]\right\}+\sum_{g=1}^G\alpha_g\Tr[\mathbf{F}^\mathrm{H}\mathbf{F}]+\sum_{g=1}^G\text{consF}_g-\tau(\Tr[\mathbf{F}^\mathrm{H}\mathbf{F}]-P_T)
$$
​	Let the first-order derivative of above equation as $\mathbf{F'}$ to 0, i.e. $\partial\mathcal{L}(\mathbf{F})/\partial\mathbf{F'}=0$. Then the globally optimal solution of $\mathbf{F}$ in iteration $n$ is derived as

$$
\mathbf{F}^{n+1}=\frac{1}{\tau-\sum_{g=1}^G\alpha_g}\sum_{g=1}^G\mathbf{U}_g
$$
​	And substituting it into the power constraint as
$$
\frac{\Tr\left[\left(\sum_{g=1}^G\mathbf{U}_g\right)^\mathrm{H}\left(\sum_{g=1}^G\mathbf{U}_g\right)\right]}{\left(\tau-\sum_{g=1}^G\alpha_g\right)^2}\leq P_T
$$
​	Obviously, the above function is a decreasing function of $\tau$, then we give two possibility

1. If this function (inequality equation) holds when $\tau=0$, then $\mathbf{F}^{n+1}=\frac{-1}{\sum_{g=1}^G\alpha_g}\sum_{g=1}^G\mathbf{U}_g$.
2. Otherwise, there must be a $\tau>0$ makes the equation held, then $\mathbf{F}^{n+1}=\sqrt{\frac{P_T}{\Tr\left[\left(\sum_{g=1}^G\mathbf{U}_g\right)^\mathrm{H}\left(\sum_{g=1}^G\mathbf{U}_g\right)\right]}}\sum_{g=1}^G\mathbf{U}_g$.

###### ·Appendix

​	In this part, we aim to prove some above  theorem to better understand the logistic behind the process of surrogating and parameters forming.

**·How the surrogate function $\tilde{R}_k(\mathbf{F,e|F}^n\mathbf{,e}^n)$ form?** 

​	We perform transformation of $R_k$, showing its hidden convexity as following
$$
\begin{align}
R_k(\mathbf{F,e})
&=\log_2\pqty{1+\frac{\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_{g'}}+\sigma_k^2}}\\
&=\log_2\pqty{1+r_{k,-g}^{-1}\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g}^2}\\
&=-\log_2\pqty{1-\pqty{r_{k,-g}+\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g}^2}^{-1}\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g}^2}\\
&=-\log_2\pqty{1-r_{k}^{-1}\vqty{t_k}^2}\\
&=-\log_2\pqty{1-r_{k}^{-1}t_kt^*_k}
\end{align}
$$
​	Where $r_{k,-g}=\sum_{g'\neq g}^G\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_{g'}}^2+\sigma_k^2$ , thus express $r_k=r_{k,-g}+\vqty{t_k}^2$, and $t_k=\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g$. Which split the fraction respectively.

​	If we consider the $r_k$, $r_{k,-g}$ and $t_k$ as a unitary variable, and replace the original variable $\mathbf{F,e}$ by these variables, transform the $R_k(\mathbf{F,e})$ into $R_k(t_k,r_k)$, this $R_k(t_k,r_k)$ is finally jointly convex in $\{t_k,r_k\}$.

​	According to Taylor Formula, the lower bound surrogate function will be derived by the first-order approximation as (I still try to figure out why differential for $t_k^*$)
$$
\begin{align}
R_k(t_k,r_k)\geq &R_k(t_k,r_k)+\underbrace{\frac{\partial R_k}{\partial t_k}|_{t_k=t_k^n}(t_k-t_k^n)}_{t_k}+\underbrace{\frac{\partial R_k}{\partial t^*_k}|_{t^*_k=t_k^{*,n}}(t_k^*-t_k^{*,n})}_{t_k^*}+\underbrace{\frac{\partial R_k}{\partial r_k}|_{r_k=r_k^n}(r_k-r_k^n)}_{r_k}\\
=&R_k(t_k,r_k)+2\Re\qty{\frac{t_k^{*,n}(t_k-t_k^n)}{r_k^n-\vqty{t_k^n}^2}}-\frac{\vqty{t_k^n}^2(r_k-r_k^n)}{r_k^n(r_k^n-\vqty{t_k^n}^2)}\\
=&R_k(t_k,r_k)+2\Re\qty{\frac{t_k^{*,n}}{r_k^n-\vqty{t_k^n}^2}t_k}-\frac{\vqty{t_k^n}^2(r_k+r_k^n)}{r_k^n(r_k^n-\vqty{t_k^n}^2)}
\end{align}
$$

​	Replace the auxiliary variable $\{t_k,r_k\}$ with original variable $\mathbf{F,e}$, we can easily derive the function as:
$$
R_k(\mathbf{F,e})\geq R_k(\mathbf{F}^n,\mathbf{e}^n)+2\Re\qty{a_k\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g}-\frac{\vqty{t_k^n}^2}{r_k^n-\vqty{t_k^n}^2}-b_k\sum_{g=1}^G\vqty{\mathbf{e}^\mathrm{H}\mathbf{H}_k\mathbf{f}_g}^2-b_k\sigma_k^2
$$






##### · *Multigroup Multicast Precoding in Massive MIMO (SDR-based scheme)*

###### ·System Model

​	We assume $G$ data streams to $G$ mutlicast groups and denote the set of groups as $\mathcal{G}=\{1,\cdots,G\}$, and also denote $j\in\mathcal{G}$ as the $j$-th stream, is of interest for $K_j$ and say these $K_j$ belong to the $j$-th group, i.e. the set $\mathcal{K}_j=\{1,\cdots,K_j\}$ is a subset of $\mathcal{K}$ and represent those users belonging to $j$-th groups.

​	Then, we model the channel by uplink pilot transmision (the modeling process can be omitted in our intention) and eventually gaining the channel response function as $\hat{\mathbf{g}}_{jk}$ for user $k$ in group $j$. And the total channel response is $\hat{\mathbf{G}}=[\hat{\mathbf{G}}_1,\cdots,\hat{\mathbf{G}}_G]$, among $\hat{\mathbf{G}}_j=[\hat{\mathbf{g}}_{j1},\cdots,\hat{\mathbf{g}}_{jK_j}]$.

​	Denote the $\mathbf{s}=[s_1,\cdots,s_G]^T$ as the signal for each groups, among $s_j$ means the interest of $j$-th group,  and $\mathbf{W}=[\mathbf{w_1,\cdots,w_G}]$ means precoding matrix, among $\mathbf{w}_j$ means joint precoding matrix vector of $j$-th group. The received signal of user $k$ in $g$ group is thus given by $y_{jk}=\mathbf{g}_{jk}^\mathrm{H}\mathbf{Ws}+n$.

###### ·SDR-based Two-Layer Scheme

​	In **[15]** a computationally efficient precoder have been proposed, but only for ideal and known CSI, so we propose another SDR-based method to solve the problem. 

​	Firstly, our precoding matrix is a $N\times G$ matrix $\mathbf{W}=[\mathbf{w_1,\cdots,w_G}]$ as above mentioned, where each $\mathbf{w}_j$ for $j$-th multicast group is 
$$
\mathbf{w}_j=\underbrace{\pqty{\mathbf{I}_N-\hat{\mathbf{G}}_{-j}\pqty{\hat{\mathbf{G}}^\mathrm{H}_{-j}\hat{\mathbf{G}}_{-j}}^{-1}\hat{\mathbf{G}}^\mathrm{H}_{-j}}}_{\text{outer layer}}\underbrace{\sum_{k=1}^{K_j}\sqrt{\mu_{jk}}\hat{\mathbf{g}}_{jk}}_{\text{inner layer}}
$$
​	Where $\hat{\mathbf{G}}_{-j}=[\hat{\mathbf{G}}_1,\cdots,\hat{\mathbf{G}}_{j-1},\hat{\mathbf{G}}_{j+1},\cdots,\hat{\mathbf{G}}_G]$ is a $N\times (K_{total}-K_j)$ matrix, and $\mu_{jk}=\sqrt{\frac{p_{jk}^{dl}}{(N-K_{total}+K_j)\gamma_{jk}}}$ with $p_{jk}^{dl}$ is the downlink power of user $k$ in group $j$ and $P_T=\sum_{g=1}^G\sum_{k=1}^{K_j}p_{jk}^{dl}$. The outer cancels the inter-group interference, while the inner aims to optimize the MMF problem.

​	In this paper, the authors do not offer us the concrete methods how they gain the above solution, thus I give up reading the following contents and read the detail paper [15] for figuring out how the $\mathbf{w}_j$ are calculated.





##### · *Reducing the Computational Complexity of Multicasting in Large-Scale Antenna Systems*

###### ·System Model

​	We denote the set of all groups as $\mathcal{G}=\{1,\cdots,G\}$ and the set of the users as $\mathcal{K}=\{1,\cdots,K\}$, and set of the users associated with group $g$ as $\mathcal{K}_g=\{1,\cdots,K_g\}$ who satisfying $\mathcal{K}_i\cap\mathcal{K}_j=\emptyset,\forall i\neq j$. Besides, in this paper, we further adapt the double index notation like $K_{g,k}$ to refer to each user $k$ in group $g$.

​	We subsequently express the channel between user $k$ in group $g$ and the BS by $\mathbf{g}_{g,k}=\sqrt{\beta_{g,k}}\mathbf{h}_{g,k}$, where $\mathbf{h}_{g,k}\sim\mathcal{CN}(\mathbf{0}_N,\mathbf{I}_N)$ is small scale fading channel and $\beta_{g,k}$ is large scale channel attenuation/path loss. Consequently, the signal received by $k$ in $g$ will be
$$
y_{g,k}=\mathbf{g}_{g,k}^\mathrm{H}\mathbf{w}_gs_g+\sum_{g'\neq g}^G\mathbf{g}_{g,k}^\mathrm{H}\mathbf{w}_{g'}s_{g'}+n_{g,k}
$$
​	Where $s_g\sim\mathcal{CN}(0,1)$ is interest of group $g$ and $n_{g,k}\sim\mathcal{CN}(0,\sigma_{g,k}^2)$ is additive Gaussian Noise, then the SINR $\gamma_{g,k}$ of $k$ in $g$ is given by
$$
\gamma_{g,k}=\frac{\vqty{\mathbf{g}_{g,k}^\mathrm{H}\mathbf{w}_j}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{g}_{g,k}^\mathrm{H}\mathbf{w}_{g'}}^2+\sigma_{g,k}^2}
$$
​	The total transmit power constraint is $\sum_{g=1}^G\norm{\mathbf{w}_g}^2$, thus the QoS and MMF problem can be respectively formulated as
$$
\begin{align}\label{QoF}
\mathcal{Q}(\eta):\ &\min_{\qty{\mathbf{w}_g}}\ \ \ \sum_{g=1}^G\norm{\mathbf{w}_g}^2\\
&\text{s.t.} \ \ \ \ \ \gamma_{g,k}\geq\eta_{g,k},\forall g,k
\end{align}
$$
​	Where $\eta_{g,k}$ are prescribed SINR for $k$ in $g$.
$$
\begin{align}
\mathcal{F}(\eta,P):&\max_{\qty{\mathbf{w}_g}}\ \ \ \min_g\min_k\alpha_{g,k}\gamma_{g,k}\\
&\text{s.t.}\ \ \ \ \ \ \sum_{g=1}^G\norm{\mathbf{w}_g}^2\leq P
\end{align}
$$
​	Where $\alpha_{g,k}$ represent the weight of $\gamma_{g,k}$ or user $k$ in group $g$, problem $\mathcal{F}$ aims to maximize the minimum SINR among all users of all groups.

###### ·Proposed Two-Layer Precoding Matrix

​	We try to eliminate the interference and solve MMF problem by using two-layer precoding matrix, respectively, among inner one $\mathbf{c}_g\in\mathbb{C}^{N-\tau_g}$ for MMF and outer one $\mathbf{F}_g\in\mathbb{C}^{N\times (N-\tau_g)}$ for interference, then the structure of our  matrix is
$$
\mathbf{w}_g=\mathbf{F}_g\mathbf{c}_g, \forall g\in\mathcal{G}
$$
**·Outer Layer (Interference Elimination)**

​	We denote $\mathbf{G}_g\in\mathbb{C}^{N\times K_g}$ as the channel response matrix of user $k$ in group $g$. We subsequently use the block-diagonalization zero-forcing to eliminate the interference $\sum_{g'\neq g}^G\mathbf{g}_{k}^\mathrm{H}\mathbf{w}_{g'}s_{g'}$.

​	**The $\mathbf{F}_g$ can be gained by singular value decomposition by $\mathbf{G}_{-g}$**, which means the columns of $\mathbf{F}_g$ can form the basis for null space of $\mathbf{G}_{-g}$, and proved in [31,32], i.e., the space spaned by $\mathbf{F}_g$ should be totally orthogonal to space spaned by $\mathbf{G}_{-g}$ in case of SINR. And the $\mathbf{F}_g$ can be obtained by following QR-based decomposition 

$$
\mathbf{G}_{-g}=\mathbf{Q}_g\mathbf{R}_g=
\begin{bmatrix}
\mathbf{Q'}_g&\mathbf{Q''}_g
\end{bmatrix}
\begin{bmatrix}
\mathbf{R'}_g\\
\mathbf{0}
\end{bmatrix}
=\mathbf{Q'}_g\mathbf{R'}_g
$$
​	Where $\mathbf{Q''}_g\in\mathbb{C}^{N\times (N-\tau_g)}$ gives the null space of $\mathbf{G}_{-g}$, then, the $\mathbf{F}_g$ can be directly derived as $\mathbf{Q''}_g$,i.e.$\mathbf{F}_g=\mathbf{Q''}_g$, satisfying $\mathbf{F}_g^\mathrm{H}\mathbf{F}_g=\mathbf{I}_{N-\tau_g}$, meanwhile. 

​	Therefore, the equivalent channel precoded by $\mathbf{F}_g$ can be $\mathbf{g'}_{g,k}=\mathbf{Q''}_g^\mathrm{H}\mathbf{g}_{g,k}$ and interference-eliminated SINR/SNR can be $\gamma_{g,k}=\frac{\vqty{\mathbf{g'}_{g,k}^\mathrm{H}\mathbf{c}_g}^2}{\sigma_{g,k}^2}$. And the power constraint can be $\sum_{g=1}^G\norm{\mathbf{w}_g}^2=\sum_{g=1}^G\norm{\mathbf{c}_g}^2$, due to above $\mathbf{F}_g^\mathrm{H}\mathbf{F}_g=\mathbf{I}_{N-\tau_g}$.

​	Having eliminated the inter-group interference of the original QoS and MMF problem can be simplified into a single group multicasting problem, correspondingly. Then, the problem is rewritten as
$$
\begin{align}
\mathcal{Q'}_g(\eta_g):\ \ \ &\min_{\{\mathbf{c}_g\}}\norm{\mathbf{c}_g}^2\\
&\text{s.t.}\frac{\vqty{\mathbf{g'}_{g,k}^\mathrm{H}\mathbf{c}_g}^2}{\sigma_{g,k}^2}\geq\eta_{g,k}
\end{align}
$$

​	Where $\eta_g\in\mathbb{C}^{K_g}$ is the vector composed by $\{\eta_{g,k}\}_{k=1}^{K_g}$, and we thus express our problem as a set of problems ,i.e.,$\mathcal{Q'}(\eta)=\qty{\mathcal{Q'}_g(\eta_g)}_{g=1}^G$. Accordingly, we also reduce the MMF problem to
$$
\begin{align}
\mathcal{F'}(\eta,P):&\max_{\qty{\mathbf{c}_g}}\ \ \ \min_g\min_k\frac{\alpha_{g,k}}{\sigma_{g,k}^2}\vqty{\mathbf{g'}_{g,k}^\mathrm{H}\mathbf{c}_g}^2\\
&\text{s.t.}\ \ \ \ \ \ \sum_{g=1}^G\norm{\mathbf{c}_g}^2\leq P
\end{align}
$$





##### · *Secure Beamforming Design for RIS-Assisted Integrated  Sensing and Communication Systems*

###### ·System Model

​	Our BS equipped with $L$ antenna serves $N$ users assisted with RIS equipped with $K$ reflection surface, meanwhile, they are placed near the users and far away from detective target to enhance the communication and appropriately reduce the effect of detective wave on communication. We assume the detective target *Eve* potentially intend to eavesdrop our communication. And CSI is known by BS.

​	Consequently, the emitted signal is given by
$$
\mathbf{x}=\sum_{i=1}^N\mathbf{w}_is_i+\mathbf{w}_\vartheta s_\vartheta
$$
​	Where $\mathbf{w}_i$ and $s_i$ is precoding matrix and intended signal of user $i$, respectively, and $\mathbf{w}_\vartheta$ and $s_\vartheta$ is precoding matrix and signal of additive AN. Then, the received signal at user $U_i$ can be expressed as
$$
y_i=\pqty{\mathbf{h}_{s,u_i}^\mathrm{H}+\mathbf{h}_{r,u_i}^\mathrm{H}\mathbf{\Theta H}_{rs}}\mathbf{x}+n_i
$$
​	Where $\mathbf{h}_{s,u_i},\mathbf{h}_{r,u_i}$ and $\mathbf{H}_{rs}$ represent channel response function from BS to user $i$, from RIS to user $i$ and from BS to RIS, respectively, while $\mathbf{\Theta}$ is reflection matrix of RIS.

​	Actually, from the radar detective respect, the performance is proportional to the echo SINR, but considering the detective target is a potential *Eve*, the echo signal should be minimum SINR. Based on above analysis, we give the following formulations to represent SINR of user $i$, target, and echo signal BS received, among target's SINR is the how much target eavesdrop communication between BS and user $i$.
$$
\begin{align}
\gamma_{i}&=\frac{\left|\left(\mathbf{h}_{s, u_{i}}^{H}+\mathbf{h}_{r, u_{i}}^{H} \boldsymbol{\Theta} \mathbf{H}_{r s}\right) \mathbf{w}_{i}\right|^{2}}{\sum_{j \neq i}^{N}\left|\left(\mathbf{h}_{s, u_{i}}^{H}+\mathbf{h}_{r, u_{i}}^{H} \boldsymbol{\Theta} \mathbf{H}_{r s}\right) \mathbf{w}_{j}\right|^{2}+\left|\left(\mathbf{h}_{s, u_{i}}^{H}+\mathbf{h}_{r, u_{i}}^{H} \boldsymbol{\Theta} \mathbf{H}_{r s}\right) \mathbf{w}_{\vartheta}\right|^{2}+\delta_{i}^{2}} \\
\gamma_{e, i}&=\frac{\left|\mathbf{h}_{e}^{H} \mathbf{w}_{i}\right|^{2}}{\sum_{j \neq i}^{N}\left|\mathbf{h}_{e}^{H} \mathbf{w}_{j}\right|^{2}+\left|\mathbf{h}_{e}^{H} \mathbf{w}_{\vartheta}\right|^{2}+\delta_{e}^{2}}\\
\Upsilon_{\mathrm{echo}}&=\frac{\sum_{j=1}^{N}\left\|\beta_{0} d_{t}^{-\sigma} \boldsymbol{\alpha}\left(\varphi_{t}\right) \boldsymbol{\alpha}^{H}\left(\varphi_{t}\right) \mathbf{w}_{j}\right\|^{2}+\left\|\beta_{0} d_{t}^{-\sigma} \boldsymbol{\alpha}\left(\varphi_{t}\right) \boldsymbol{\alpha}^{H}\left(\varphi_{t}\right) \mathbf{w}_{\vartheta}\right\|^{2}}{\delta_{t}^{2}} \\
\end{align}
$$
​	Naturally, the secrecy rate (defined to measure secure transmission at user $i$) is given by
$$
C_s^i=\bqty{\log_2 \pqty{1+\gamma_i}-\log_2\pqty{1+\gamma_{e,i}}}^+
$$
​	We aim to optimize the matrices of precoding and reflection to maximize the secrecy rate. Eventually, having set the threshold $\Upsilon_\mathrm{th}$ as the SINR of echo signal, the problem is reorganized as
$$
\begin{align}
P1:\ \ \ \max_{\qty{\mathbf{w}_j}_{j=1}^N,\mathbf{w}_\vartheta,\mathbf{\Theta}}&\sum_{i=1}^NC_s^i\\
\text{s.t.}\ \ \
&\Upsilon_{\mathrm{echo}}\geq \Upsilon_\mathrm{th}\\
&\sum_{j=1}^N\norm{\mathbf{w}_j}^2\leq P_s\\
&\theta_n\in[0,2\pi],\forall n
\end{align}
$$

###### ·Joint Beamforming Optimization

​	We first convert the original problem into **two subproblems** and simplify them by reformulate objective function and constraints aided with Taylor expansion and SDR. Then, we further solve them by SCA and AO strategy and finally gain the suboptimal solution.

**·Optimizing $\qty{\mathbf{w}_j}_{j=1}^N$ and $\mathbf{w}_\vartheta$ With Given $\mathbf{\Theta}$**

​	To reduce the complexity of formulations, we replace $\pqty{\mathbf{h}_{s, u_{i}}^{H}+\mathbf{h}_{r, u_{i}}^{H} \boldsymbol{\Theta} \mathbf{H}_{r s}}$ with $\mathbf{h}_i$ as the unit channel response function and $\mathbf{w}_{N+1}=\mathbf{w}_\vartheta$, the objective function is converted into following form
$$
C_s^i=\underbrace{\log_2\pqty{\frac{\sum_{j=1}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2}{\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2}}}_{\text{original rate}}-\underbrace{\log_2\pqty{\frac{\sum_{j=1}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2}{\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2}}}_{\text{eavesdrop rate}}
$$
​	We introduce exponential auxiliary variables to simplify the $\log$ form objective function as following
$$
e^{\tau_i}=\sum_{j=1}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2\\
e^{\varepsilon_i}=\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2\\
e^{u_i}=\sum_{j=1}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2\\
e^{v_i}=\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2
$$
​	Then, the objective function will be $\max_{\qty{\tau_i,\varepsilon_i,u_i,v_i}_{i=1}^N}\sum_{i=1}^N\pqty{\tau_i-\varepsilon_i-u_i+v_i}$, which is obviously a convex function, and the constraints of $\tau_i,\varepsilon_i,u_i,v_i$ are expressed as
$$
e^{\tau_i}\leq\sum_{j=1}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2\\
e^{\varepsilon_i}\geq\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2\\
e^{u_i}\geq\sum_{j=1}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2\\
e^{v_i}\leq\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2
$$
​	The inequalities only hold equality at the optimum point, proved by the monotonicity of functions. However, the formulation of $e^{\varepsilon_i},e^{u_i}$ are still non-convex, then, we makes it convex by Taylor expansion at iterate point $k$ (means iterate $k$ times)
$$
\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2\leq e^{\varepsilon_i^k}\pqty{\varepsilon_i-\varepsilon_i^k+1}\\
\sum_{j=1}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2\leq e^{\varepsilon_i^k}\pqty{u_i-u_i^k+1}
$$
​	Furthermore, to convert the quadratic form of echo signal SINR, we define $\mathbf{A}(\varphi_t)=\beta_0d_t^{-\sigma}\alpha(\varphi_t)\alpha^\mathrm{H}(\varphi_t)$, then the SINR will be rewritten as
$$
\sum_{j=1}^{N+1}\Tr\pqty{\mathbf{A}(\varphi_t)\mathbf{W}_j\mathbf{A}^\mathrm{H}(\varphi_t)}\geq \delta_t^2\Upsilon_{\mathrm{th}}
$$
​	Eventually, the problem is reformulated as
$$
\begin{align}
P2:\ \ \max_{\qty{\tau_i,\varepsilon_i,u_i,v_i}_{i=1}^N,\qty{\mathbf{W}_j}_{j=1}^{N+1}}\ \ \ &\sum_{i=1}^N\pqty{\tau_i-\varepsilon_i-u_i+v_i}\\
\text{s.t.}\ \ \ \ &\sum_{j=1}^{N+1}\Tr(\mathbf{W}_j)\leq P_s,\\
&\mathbf{W}_j\succeq0,\ \ \ \text{Rank}(\mathbf{W}_j)=1\\
&e^{\tau_i}\leq\sum_{j=1}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2,\\
&e^{v_i}\leq\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2,\\
&\sum_{j\neq i}^{N+1}\Tr(\mathbf{H}_i\mathbf{W}_j)+\delta_i^2\leq e^{\varepsilon_i^k}\pqty{\varepsilon_i-\varepsilon_i^k+1},\\
&\sum_{j=1}^{N+1}\Tr(\mathbf{H}_e\mathbf{W}_j)+\delta_e^2\leq e^{\varepsilon_i^k}\pqty{u_i-u_i^k+1},\\
&\sum_{j=1}^{N+1}\Tr\pqty{\mathbf{A}(\varphi_t)\mathbf{W}_j\mathbf{A}^\mathrm{H}(\varphi_t)}\geq \delta_t^2\Upsilon_{\mathrm{th}}
\end{align}
$$
​	Among these constraints, the non-convex constraint $\text{Rank}(\mathbf{W}_j)=1$ can be strategically relaxed, and this problem $P2$ can be solved by CVX solver.

##### · *Joint AP Selection and Power Allocation for Unicast-Multicast Cell-Free Massive MIMO*

###### · Abstract

​	This paper focus on a cell-free massive MIMO system enabled to simultaneously unicast and multi-group multicast. The precise achievable downlink spectral efficiency (SE) of unicast and multicast users are derived from ZF and maximum precoding design. This paper formulates a WSE maximization problem with constraints like transmit power, fronthual capacity and QoS to jointly optimize the APs selection and power allocation, then, we reformulate this non-convex problem into an easier structure and an accelerate projected gradient (APG)-based algorithm is adopted to obtain the near-optimal solutions. As comparison, SCA-based algorithm is implemented, and the results indicate the proposed joint optimization approach enhances the WSE in various environments. Moreover, APG-based algorithm reduce computational complexity hugely with considerable performance.

###### · System Model

​	We consider a Cell-Free massive MIMO (CF-mMIMO) system with joint unicast and multi-group multicast transmission, this system consists of $N$ APs, each equipped with $L$ antennas, simultaneously serving $U$ unicast users and $M$ multicast groups, where $m$-th group includes $K_m$ users. The sets of $N$ APs, $U$ unicast users, $M$ multicast groups and $K_m$ users of $m$-th are denoted by $\mathcal{N,U,M}$ and $\mathcal{K}_m$, respectively. The channel vector between the $u$-th unicast user and the $n$-th AP, $u\in\mathcal{U},n\in\mathcal{N}$ is
$$
\mathbf{c}_{n,u}=\beta_{n,u}^{1/2}\mathbf{h}_{n,u}\in\mathbb{C}^{L\times1},
$$
and the channel between the $k_m$-th multicast user of the $m$-th multicast group and $n$-th AP $k_m\in\mathcal{K}_m,m\in\mathcal{M},n\in\mathcal{N}$ is
$$
\mathbf{t}_{n,m,k}=\bar{\beta}_{n,m,k}^{1/2}\mathbf{h}_{n,m,k}\in\mathbb{C}^{L\times1}
$$
where $\beta_{n,u}$ and $\bar{\beta}_{n,m,k}$ are the according large-scale fading coefficients, and $\mathbf{h}_{n,u}\sim\mathcal{CN}\pqty{\mathbf{0,I}_L}$ and $\mathbf{h}_{n,m,k}\sim\mathcal{CN}\pqty{\mathbf{0,I}_L}$ are the according small-scale fading vectors.

​	This system is assumed to work under the reciprocity-based TDD protocol, where the channels remains unchanged during a coherence interval $T$. CSI is needed by APs, and the pilot to unicast is orthogonal and to all users in each multicast group is shared, due to limited coherence interval and sharing mechanism of multicast, thus, $U+M$ orthogonal pilots are needed. Let $\phi_u\in\mathbb{C}^{\tau\times1},\norm{\phi_u}^2=1$ be the pilot of $u$-th unicast user and $\varphi_m\in\mathbb{C}^{\tau\times1},\norm{\varphi_m}^2=1$ be the pilot of multicast users in $m$-th multicast group, while $U+M\leq\tau\leq T$. Therefore, $\phi_u^\mathrm{H}\phi_{u'}=0,\phi_u^\mathrm{H}\varphi_m=0,\varphi_m^\mathrm{H}\varphi_{m'}=0,\forall u\neq u',m\neq m'$. The received signal at $n$-th AP during uplink training can be
$$
\mathbf{Y}_{n,p}=\sqrt{\tau p_\text{ul}}\sum_{u=1}^U\mathbf{c}_{n,u}\phi_u^\mathrm{H}+\sqrt{\tau p_\text{ul}}\sum_{m=1}^M\sum_{k=1}^{K_m}\mathbf{t}_{n,m,k}\varphi_m^\mathrm{H}+\Psi_{n,p}
$$
where $\Psi_{n,p}\in\mathbb{C}^{L\times\tau}$ is the addictive white Gaussian noise and $p_\text{ul}$ is the uplink transmit power. To estimate $\mathbf{c}_{n,u}$, we can correspondingly project the received signal $\mathbf{Y}_{n,p}$ onto the pilot of $u$-th unicast user
$$
\check{\mathbf{y}}_{n,p,u}=\mathbf{Y}_{n,p}\phi_u=\sqrt{\tau p_\text{ul}}\mathbf{c}_{n,u}+\psi'_{n,p}
$$
where $\psi'_{n,p}=\Psi_{n,p}\phi_u\sim\mathcal{CN}\pqty{\mathbf{0,I}_L}$ can be easily gained. Each AP can estimate the unicast user channel and further minimize the backhaul signaling. The MMSE estimate of $\mathbf{c}_{n,u}\in\mathbb{C}^{L\times1}$ is
$$
\begin{align}
\hat{\mathbf{c}}_{n,u}&=\frac{\mathbb{E}\qty{\mathbf{c}_{n,u}\check{\mathbf{y}}_{n,p,u}^\mathrm{H}}}{\mathbb{E}\qty{\check{\mathbf{y}}_{n,p,u}\check{\mathbf{y}}_{n,p,u}^\mathrm{H}}}\cdot\check{\mathbf{y}}_{n,p,u}=\frac{\mathbb{E}\qty{\sqrt{\tau p_\text{ul}}\mathbf{c}_{n,u}\mathbf{c}_{n,u}^\mathrm{H}}}{\mathbb{E}\qty{\tau p_\text{ul}\mathbf{c}_{n,u}\mathbf{c}_{n,u}^\mathrm{H}+\psi_{n,p}'\psi_{n,p}'^\mathrm{H}}}\cdot\check{\mathbf{y}}_{n,p,u}\\
&=\frac{\sqrt{\tau p_\text{ul}}\beta_{n,u}\mathbb{E}\qty{\mathbf{I}_L}}{\pqty{\tau p_\text{ul}\beta_{n,u}+1}\mathbb{E}\qty{\mathbf{I}_L}}\cdot\check{\mathbf{y}}_{n,p,u}=\frac{\sqrt{\tau p_\text{ul}}\beta_{n,u}}{\tau p_\text{ul}\beta_{n,u}+1}\cdot\check{\mathbf{y}}_{n,p,u}
\end{align}
$$

due to $\mathbf{c}_{n,u}\psi_{n,p}'^\mathrm{H}=\psi_{n,p}'\mathbf{c}_{n,u}^\mathrm{H}=0$. Then, the variance of the estimated channel $\hat{\mathbf{c}}_{n,u}$ will be
$$
\gamma_{n,u}=\mathbb{E}\qty{\vqty{\bqty{\hat{\mathbf{c}}_{n,u}}_l}^2}=\frac{\tau p_\text{ul}\beta_{n,u}^2}{\tau p_\text{ul}\beta_{n,u}+1},
$$
and the estimation error of $\mathbf{c}_{n,u}$ is $\tilde{\mathbf{c}}_{n,u}=\mathbf{c}_{n,u}-\mathbf{\hat{c}}_{n,u}\sim\mathcal{CN}\pqty{\mathbf{0},\pqty{\beta_{n,u}-\gamma_{n,u}}\mathbf{I}_L}$.

Similarly, the MMSE estimate of the $k$-th multicast user $\mathbf{t}_{n,m,k}$ is given by
$$
\hat{\mathbf{t}}_{n,m,k}=\frac{\sqrt{\tau p_\text{ul}}\bar{\beta}_{n,m,k}}{\tau p_\text{ul}\sum_{t=1}^{K_m}\bar{\beta}_{n,m,t}+1}\check{\mathbf{y}}_{n,p,m}
$$
where $\check{\mathbf{y}}_{n,p,m}$ is the received signal projection onto the pilot of $m$-th multicast group like $\check{\mathbf{y}}_{n,p,u}$, and $\check{\mathbf{y}}_{n,p,m}=\mathbf{Y}_{n,p}\varphi_m=\sqrt{\tau p_\text{ul}}\sum_{k=1}^{K_m}\mathbf{t}_{n,m,k}+\psi''_{n,p}$, while $\psi''_{n,p}=\Psi_{n,p}\varphi_m\sim\mathcal{CN}\pqty{\mathbf{0,I}_L}$. Thus, the estimation error of $\mathbf{t}_{n,m,k}$ will similarly be $\tilde{\mathbf{t}}_{n,m,k}=\mathbf{t}_{n,m,k}-\hat{\mathbf{t}}_{n,m,k}\sim\mathcal{CN}\pqty{\mathbf{0},\pqty{\bar{\beta}_{n,m,k}-\bar{\gamma}_{n,m,k}}\mathbf{I}_L}$, where
$$
\bar{\gamma}_{n,m,k}=\mathbb{E}\qty{\vqty{\bqty{\hat{\mathbf{t}}_{n,m,k}}_l}^2}=\frac{\tau p_\text{ul}\bar{\beta}^2_{n,m,k}}{\tau p_\text{ul}\sum_{t=1}^{K_m}\bar{\beta}_{n,m,t}+1}
$$

Considering the co-pilot strategy in [32], we gain
$$
\hat{\mathbf{t}}_{n,m}=\sum_{k=1}^{K_m}\hat{\mathbf{t}}_{n,m,k}=\hat{\mathbf{t}}_{n,m,k}=\frac{\sqrt{\tau p_\text{ul}}\sum_{k=1}^{K_m}\bar{\beta}_{n,m,k}}{\tau p_\text{ul}\sum_{t=1}^{K_m}\bar{\beta}_{n,m,t}+1}\check{\mathbf{y}}_{n,p,m}
$$
and regard this as the channel estimate of the $m$-th multicast group, with mean square as
$$
\zeta_{n,m}=\mathbb{E}\qty{\vqty{\bqty{\hat{\mathbf{t}}_{n,m}}_l}^2}=\hat{\mathbf{t}}_{n,m,k}=\frac{\pqty{\sqrt{\tau p_\text{ul}}\sum_{k=1}^{K_m}\bar{\beta}_{n,m,k}}^2}{\tau p_\text{ul}\sum_{t=1}^{K_m}\bar{\beta}_{n,m,t}+1}
$$



#### · New Mathematical Tools

##### · Precoder Design for User-Centric Network Massive MIMO: A Symplectic Optimization Approach

###### · Abstract

​	In this paper, we design a precoder for a massive MIMO-based UCN by symlectic optimization, where the system serves users by partial BSs rather than all of them. The dimension of precoders is reduced compared to conventional massive MIMO, simplifying the precoders in practice. To avoid computational assumption on matrix inversion, symplectic optimization framework is adopted, which is based on dissipative Hamiltonian dynamical systems. However, to better fit this framework, we transform the received model into the real field and reformulate the WSR problem. The objective function is regarded as the potential energy of dynamical system. Due to the energy dissipation, the continuous dynamical system converges to a minimal potential energy state. By discretizing the continuous system while preserving the symplectic structure, we gain a iterative precoder design method. Finally, our proposed method is highly computationally efficient, while simulations prove the proposed symplectic optimization-based method outperform the WMMSE precoder in massive MIMO-based UCN.

###### · Introduction

​	Massive MIMO 

###### · System Model

​	Consider a massive MIMO-based UCN consisting of $K$ UTs and $B$ BSs, each BS equipped with a uniform planar array with $M_t=M_v\times M_h$ antennas, while each UT equipped with $M_r$ antennas. This communication system work in the TDD mode. Assume that BSs are synchronized and interconnected by backhual, thus, the BSs can jointly coherently transmit. Let $\mathcal{S}_B=\qty{1,2,\cdots,B}$ denotes the set of BSs, while $\mathcal{S}_U=\qty{1,2,\cdots,K}$ denotes that of UTs, and the set $\mathcal{B}_k=\qty{k_1,k_2,\cdots,k_{B_k}}$ is the subset of $\mathcal{S}_B$ to serve specific user $k$, the set $\mathcal{U}_l=\qty{1,2,\cdots,l_{U_l}}$ is also the subset of $\mathcal{S}_U$ and served by specific BS $l$. This user-centric rule enables each UT to receive signals from the BS with best channel quality instead of the traditional cell concept. Fig 1 illustrates the massive MIMO-based UCN system.

​	Let $x_k$ denotes the symbol transmitted to the $k$-th UT with $\mathbb{E}\qty{\vqty{\mathbf{x}}^2}=\mathbf{I}_K$, while $\mathbf{x}=\bqty{x_1,x_2,\cdots,x_K}^\mathrm{T}$. The channel vector from the $l$-th BS to the $k$-th user is denoted as $\mathbf{h}_{l,k}\in\mathbb{C}^{M_t\times M_r}$, and $\mathbf{p}_{l,k}$ denoted as the precoder vector for the transmission from $l$-th BS to $k$-th UT, among $l\in\mathcal{B}_k$. Then the received signal at $k$-th UT will be
$$
y_k=\sum_{l\in\mathcal{B}_k}\mathbf{h}_{l,k}^\mathrm{H}\mathbf{p}_{l,k}x_k+\sum_{l\in\mathcal{B}_k}\sum_{t\in\mathcal{U}_l,t\neq k}\mathbf{h}_{l,k}^\mathrm{H}\mathbf{p}_{l,t}x_t+\sum_{m\notin\mathcal{B}_k}\sum_{t\in\mathcal{U}_m}\mathbf{h}_{m,k}^\mathrm{H}\mathbf{p}_{m,t}x_t+z_k
$$
where $z_k\sim\mathcal{CN}\pqty{0,\sigma_z^2}$ is the AWGN. Let $z_k'$ denote the interference-plus-noise of $k$-th UT defines as
$$
\begin{align}
z_k'
&=\sum_{l\in\mathcal{B}_k}\sum_{t\in\mathcal{U}_l,t\neq k}\mathbf{h}_{l,k}^\mathrm{H}\mathbf{p}_{l,t}x_t+\sum_{m\notin\mathcal{B}_k}\sum_{t\in\mathcal{U}_m}\mathbf{h}_{m,k}^\mathrm{H}\mathbf{p}_{m,t}x_t+z_k\\
&=\sum_{t\neq k}\sum_{m\in\mathcal{B}_t}\mathbf{h}_{m,k}^\mathrm{H}\mathbf{p}_{m,t}x_t+z_k
\end{align}
$$
(Here, the author transfers the view from BSs to UTs, considering all the UT and the served BSs set except $k$-th UT)

whose covariance is given by
$$
r_k=\sum_{t\neq k}\pqty{\sum_{m\in\mathcal{B}_t}\mathbf{h}_{m,k}^\mathrm{H}\mathbf{p}_{m,t}}\pqty{\sum_{m\in\mathcal{B}_t}\mathbf{h}_{m,k}^\mathrm{H}\mathbf{p}_{m,t}}^\mathrm{H}+\sigma_z^2
$$
then, the rate of $k$-th UT is obtained as 
$$
\mathcal{R}_k=\log_2\pqty{1+r_k^{-1}\pqty{\sum_{l\in\mathcal{B}_k}\mathbf{p}_{l,k}^\mathrm{H}\mathbf{h}_{l,k}}\pqty{\sum_{l\in\mathcal{B}_k}\mathbf{h}_{l,k}^\mathrm{H}\mathbf{p}_{l,k}}}.
$$
​	For simplicity, let $\mathbf{p}_l=\bqty{\mathbf{p}_{l,1}^\mathrm{T},\mathbf{p}_{l,2}^\mathrm{T},\cdots,\mathbf{p}_{l,K}^\mathrm{T}}^\mathrm{T}$ and $\mathbf{p}=\bqty{\mathbf{p}_1^\mathrm{T},\mathbf{p}_2^\mathrm{T},\cdots,\mathbf{p}_L^\mathrm{T}}^\mathrm{T}$. Since each BS has its individual power constraint in massive MIMO-based UCN, we organize a WSR-maximization precoder design method
$$
\begin{align}
\max_{\mathbf{p}}\ \ \ \ &\sum_{k\in\mathcal{S}_U}\omega_k\mathcal{R}_k,\\
\text{s.t.}\ \ \ \ &\sum_{k\in\mathcal{U}_l}\mathbf{p}_{l,k}^\mathrm{H}\mathbf{p}_{l,k}\leq\rho_l,\forall l\in\mathcal{L}
\end{align}
$$
where $\omega_k$ is the weight of $k$-th user and $\rho_l$ is the power constraints of $l$-th BS.

​	This WSR-maximization problem can be directly solved by the iterative WMMSE.

​	For simplicity, we define $A_{k,t}=\sum_{m\in\mathcal{B}_t}\mathbf{h}_{m,k}^\mathrm{H}\mathbf{p}_{m,t}$, and the $r_k$ can be rewritten as
$$
\check{r}_k=\sum_{t\in\mathcal{S}_U,t\neq k}A_{k,t}A_{k,t}^*+\sigma_z^2,
$$
and the WMMSE precoder is derived as
$$
\mathbf{p}_{l,k}^\text{WMMSE}=\frac{\omega_k\mathbf{h}_{l,k}u_k^*\mathrm{W}_k-\sum_{j\in\mathcal{S}_U}\omega_j\mathbf{h}_{l,j}u_j^*\sum_{m\neq l}\mathbf{h}_{m,j}^\mathrm{H}\mathbf{p}_{m,k}}{\sum_{j\in\mathcal{S}_U}\omega_j\mathbf{h}_{l,j}u_j^*\mathrm{W}_ju_j\mathbf{h}_{l,j}^\mathrm{H}+\lambda_k\mathbf{I}_{M_t}}
$$
where $u_k=A_{k,t}^*/\check{r}_k$ is the receive coefficient for $k$-th UT, $\mathrm{W}_k=1+A_{k,k}^*A_{k,k}/r_k$ is the MMSE weight, and $\lambda_k$ denotes the Lagrange multiplier, computed by bisection method。

​	As shown in WMMSE precoder equation, the WMMSE method involves the inversion of an $M_t\times M_t$ matrix, which resulted in high computational complexity, significantly limiting the practical application of massive MIMO-based UCN system. To avoid this, conventional gradient descent (GD) methods directly solve WSR problem. Here, we define $g\pqty{\hat{\mathbf{p}}}=-f\pqty{\hat{\mathbf{p}}}$. Under this method, the conventional update formula can be given by
$$
\mathbf{p}_{n+1}^\text{GD}=\mathbf{p}_n^\text{GD}-\alpha\grad g\pqty{\mathbf{p}_n^\text{GD}}
$$
where $\grad g\pqty{\mathbf{p}_n^\text{GD}}$ denotes the negative gradient of WSR-maximization objective function with variable $\mathbf{p}_n^\text{GD}$ and $\alpha$ is the step size. However, the convergence performance of this conventional GD method is unsatisfied, while a better method NAGD is widely used to enhance the performance. The NAGD introduce a momentum term, significantly speeding up the iterative process and improving the convergence, which has been used in satellite communication precoder design with user-centric rule. The NAGD is as follows
$$
\begin{align}
&\mathbf{q}_n=\mathbf{p}_n^\text{NAGD}+\mu\pqty{\mathbf{p}_n^\text{NAGD}-\mathbf{p}_{n-1}^\text{NAGD}}\\
&\mathbf{p}_{n+1}^\text{NAGD}=\mathbf{q}_n-\alpha\grad g\pqty{\mathbf{p}_n^\text{NAGD}}
\end{align}
$$
where $\mu$ is hyperparameter, $\alpha$ is step size in NAGD, and $\mathbf{q}$ denotes the novel momentum.

​	In the GD method, it's crucial to choose step size for guaranteeing convergence performance. Generally, line search methods are widely used to enhance convergence, while Armijo condition is a criterion to test whether the candidate step size is proper. However, the line search spends too much computational resources on each iteration of objective function to find a fit $\alpha$.

​	Symplectic optimization is an advanced novel method proposed to solve such a mathematical optimization problem, which relates the optimization problems with a dissipative dynamical system, where the potential energy is the objective function and a kinetic energy term is also included. Owing to energy dissipation, the continuous dissipative dynamical system always converges to a minimal potential energy, which is also a minimal value of the optimization problem. By using discretization that keeps the symplectic structure, and symplectic optimization method that keeps the properties of original continuous dynamical system is obtained, which means this symplectic optimization-based algorithm will run faster than conventional GD and more likely to escape local optimal point thanks to kinetic energy, then this algorithm might outperform the WMMSE who may be trapped in a local optimal point. We adopt such a symplectic optimization-based to precoder design for our considered massive MIMO-based UCN system in the following research.




#### · Semantic Communication

##### · *CSDN Blog Based on Paper -- Semantic Communications: overview, open issues, and Future Research Direction*

​	This notes is originate from a CSDN blog and the paper published in 2022 by Luo.

###### · Introduction

​	The semantic communication is generally divided into three parts

1. Technical layer: this layer is the first layer defined by Shannon Theory, concentrated on the transmission of signal
2. Semantic layer: this layer aims to detect the features of message, then extract and transmit the semantic information
3. Effective layer: this layer contains a Knowledge Base (KB) to ensure the communication efficiency

​	The semantic communication is proposed as a intelligent communication scheme concentrated on the transmission of semantic information, rather than security communication required accurate bit stream. The core of semantic communication is to extract the meaning of message at transmitter and translate it as a similar message at receiver, even definitely accurate (maybe some fault but understood).

​	some correlate keywords: Syntax, Polysemy, Synonym, Dialect

###### · Difference

​	The principle is to extract the characters or features from the information source and explain it at destination.

​	Compared to traditional communication, whose data will be source/channel/physical encoded at transmitter and the similar process at receiver to ensure the robustness. However, the semantic communication system is not just a system equipped with the above function but also a system contained AI agent at transmitters and receivers respectively to carry out intelligent algorithm for message translation. Besides, this agent is also required to percept the environment (in case to real-time update KB), and agent at transmitter will extract the feature, while that at receiver will understand and restore (even speculate) the meaning of message.

###### · Noise and Errors

​	**Noise:** Different with the traditional communication, the semantic communication system is intervened not just by imperfect noise channel but also by semantic noise due to ambiguity.

​	**Err:** The err usually caused by mismatching KB between source and destination (semantic aspect) or bit err due to physical and semantic noise (technical aspect).

###### · Encoding and Decoding





##### · *Deep Learning Enabled Semantic Communication Systems*

###### · Abstract

​	Inspired by NLP and application of Deep-Learning on communication, they propose semantic communication to provide a new aspect from **semantic level**, which can better analyze and understand the texts. This paper propose a Transformer based semantic communication system for transmission to improve the performance by transmitting at semantic level rather than bit or symbol. Besides, learning also assistant accelerate model training and guarantee generality. This paper define the metric to judge the sentences similarity with more considering semantics exchange in various scenario (low SINR).

###### · Introduction

​	Generally speaking, the communication can be categorized into 3 following layers:

	1. Symbols Transmission
	1. Semantic Exchange of Transmitted Symbol
	1. Effect of Semantic Information Exchange (Application)

​	The layer Ⅰ concerns about the symbol transmission, accuracy is measured at the bit level. The layer Ⅱ concentrated on semantic information transmission and interpret, e.g. semantic communication. The layer Ⅲ deals with the effects of communication in specific task environment.

​	From 1G to 5G, people pay their attention on the approach to improving the accuracy and efficiency of symbol transmission, measuring the performance by BER and SER. However, while the development of communication is approaching the Shannon Limit, there are still more various application from automobile to health-care demanding and generating more data and needing to be supported  massive connectivity over limited spectrum resources with lower latency, which challenged the traditional communication. Luckily, the semantic communication can process data in semantic domain by extracting meanings of original data, only proceeding and transmitting the useful, relevant and essential part, which hugely compressed the data while reserving the meanings. Moreover, considering the application, the semantic communication will be more robust in lower SINR to fit the demand on high reliability.

​	There were a few researchers already studied semantic communication since several decades ago. GMSC  (generic model of semantic communication) is applied in [8] to construct a lossless compression theory, opening the door for data compression at semantic level. The researchers of [9] proposed an E2E framework, integrally solving the semantic inference and physical layer problems, but only in the words meaning not sentences, while these works still provide us insights for future more efficient semantic communication systems. The development of DL and NLP inspires us its possible application on SC, when considered SC systems mainly focus on the joint semantic-channel coding and decoding. Conclusively, we are confront with the following questions:

1. *How to define the meaning behind the bits?*
2. *How to measure the semantic error of sentences?*
3. *How to jointly design the semantic and channel coding?*

​	In this paper, we proposed DeepSC based on machine translation in NLP for physical layer to address the above challenges, our contribution is as follows:

1. ***Proposed a framework*** can efficiently extract semantic information from texts, and ***designed a joint semantic-channel coding*** to cope with noise and distortion.
1. Our transceiver is composed by semantic encoder/decoder and channel encoder/decoder, while receiver is optimized with two loss function (cross-entropy and  mutual information). Besides, ***a new metric*** is also proposed to measure the performance of DeepSC at semantic level.
1. To the generality of DeepSC, we use the deep transfer learning to accelerate the model pre-training. With the pre-trained model, our DeepSC can recognize various knowledge and more easily recover the semantic information.
1. Our proposed DeepSC is much better than the conventional system, even robustness at SNR regime.

###### · Related  Work

​	This section briefly introduce the related work on the E2E physical layer system and the DNN adopted in NLP.

1. E2E Physical Layer Communication Systems

​	Some pioneering works show the potentiality of E2E system by adopting the auto-encoder in DL and removing the block structure, demonstrating the advantage than uncoded BPSK and Hamming coded BPSK in terms of BER. On this basis, more work researched on the missing channel gradient during training. Among, a DNN-based training processing two-phase are proposed, the transceiver is trained by stochastic channel model, while the receiver is fine-tuned under real channel, the RL also exploited in [19] to assistant to construct the channel gradient under unknown model, and also better performing than the DQPSK in real. A GAN is used in [20] to use DNN to represent the distortion so that adapting unknown channel model during training the E2E system. Conclusively, meta-learning with few pilots (like the model in transceiver, which is mentioned above and trained by stochastic channel)  has been developed for training fast with little data.

​	Besides, the old BER is unfitted to new type of data, like texts and images, so they adopt word-error rate and peak signal-to-noise-ratio (PSNR) for measuring the accuracy of semantic information recovery rather than digital bits.

2. Semantic Representation in NLP

​	Thanks to improvement of NLP, we have a strong tool to cope with the semantic information. Recalling the development of NLP, the **step Ⅰ** is **probability prediction**, predicting the next word by the above, with disadvantage on long sentence and syntax, **step Ⅱ** construct the **word2vec** model, trying to capture and describe the relationship among words, solving the long-term capture while still can't understand the syntax, **step Ⅲ** proposed a **underlying meaning** of texts by adopting various DL, able to extract the syntax and the semantic meaning of long sentence. In step Ⅲ, a deep contextualized word representation is proposed, modeling complex characteristics of word usage (syntax and semantics) and the difference across linguistic contexts. However, the model fits to specific tasks only, needing to redesign when task changes. In [27], a general word representation model, BERT, is developed to provide word vectors for various NLP without more redesign.

3. Comparison of State-of-Art NLP Techniques

​	There is three types of neural networks used for NLP tasks. Firstly, RNN, despite the RNN can effectively capture the syntax information and learn the whole sentence, the long sentence tasks are big trouble for RNN, also with the disadvantages on parallel computing. Secondly, CNN, thanks to the convolution kernel, the CNN can carry on parallel tasks, but the size of kernel is too small to extract the semantic information effectively, even a deeper network. However, there is also a third path, FCN, based on Transformer, the model pay more attention to the useful semantic information on NLP tasks, combining the virtue of CNN and RNN.

###### · System Model and Problem Formulation

​	This paper design a model consist of two levels. Among semantic level extract the semantic information for encoding and decoding, given the certain KB (trained under various data), while the transmission level is charge for the accuracy of exchange.

​	By the way, two kind of noise are defined. The Semantic Noise originates from the ambiguous interpretation in words. The Physical Channel Noise is caused by physical channel impairment.

​	As the Fig. 1, the transmitter intend to send a message $\mathbf{s}$, then maps it into a complex symbol stream $\mathbf{x}$ and passes it through the physical channel with noise and distortion. Next, the $\mathbf{y}$ is received and decoded at receiver to estimate the original sentence $\mathbf{s}$. This paper jointly design the transmitter and receiver with DNN to train model with variable-length sentence input.

​	Particularly, assume the sentence $\mathbf{s}=\bqty{w_1,w_2,\dots,w_L}$, where the $w_l$ means the $l$-th word in sentence $\mathbf{s}$. The transmitter is consist of semantic encoder and channel encoder to guarantee the semantic information extraction and transmission. The encoder can be expressed as:
$$
\mathbf{x}=C_\alpha\pqty{S_\beta\pqty{\mathbf{s}}}
$$
where $\mathbf{x}\in\mathbb{C}^{M\times 1}$, and $S_\beta\pqty{\cdot}$ is the semantic encoder network with the parameter set $\beta$ and $C_\alpha\pqty{\cdot}$ is the channel encoder with the parameter set $\alpha$. Then, if the $\mathbf{x}$ is sent, the signal received at receiver will be
$$
\mathbf{y}=h\mathbf{x}+\mathbf{n}
$$
where $\mathbf{y}\in C^{M\times 1}$, and $h$ represents the Rayleigh fading channel with $\mathcal{C}N\pqty{0,1}$ and $\mathbf{n}$ with $\mathcal{C}N\pqty{0,\sigma_n^2}$. For our E2E training of en/decoder, the channel is supposed to be allowed back-propagation, while physical channel is formulated by neural networks, i.e. the simple neural networks can assistant model the physical channel, like AWGN channel, as for fading channel, more complicated networks are used. In this paper, we only consider the AWGN channel and Rayleigh fading channel for simplicity, when we concentrate on coding and decoding.

​	Correspondingly, this paper design the receiver to recover the symbols and sentence, considering the process of encoding and the circumstance of transmission channel, finally as follows
$$
\hat{\mathbf{s}}=S_\chi^{-1}\pqty{C_\delta^{-1}\pqty{\mathbf{y}}}
$$
where the $\hat{\mathbf{s}}$ is the recovered sentence, and $C_\delta^{-1}\pqty{\cdot}$ and $S_\chi^{-1}\pqty{\cdot}$ are channel decoder and semantic decoder with parameters set $\delta$ and $\chi$, respectively.

​	This paper aims to minimize the semantic error while reduce the number of symbols. However, we face two challenges. Firstly, how to design the joint semantic-channel coding. Secondly, semantic transmission is not considered in the traditional communication. Especially to the second challenge, even if several bit error could lead to misunderstanding due to partial semantic information missed. Thus, to gain a successful recovery at semantic level, this paper design a joint semantic and channel coding to eliminate the difference between $\mathbf{s}$ and $\hat{\mathbf{s}}$, which is realized by a new DNN framework. The loss-function is designed as cross-entropy to measure the difference, and formulated as
$$
\mathcal{L}_{\mathrm{CE}}\pqty{\mathbf{s,\hat{s};\alpha,\beta,\gamma,\delta}}=-\sum_{l=1} q(w_l)\log\pqty{p\pqty{w_l}}+\pqty{1-q(w_l)}\log\pqty{1-p(w_l)}
$$
where $q(w_l)$ is the real probability that the $l$-th word in $\mathbf{s}$, while $p(w_l)$ is the predicted probability that the $l$-th word in predicted $\hat{\mathbf{s}}$. By reducing the value of loss-function, we steadily approach the original sentence $\mathbf{s}$ by learning the word distribution. Moreover, joint training and designing will assistant to pay more attention to protecting the relevant semantic information, neglecting the irrelevant ones, compared with the separate designing treating the information equally.

​	Maximize the capacity or data rate is the goal of any communication system, mutual information provides more information to train a receiver than BER, similarly defined the transmitted and received symbol as $\mathbf{x}$ and $\mathbf{y}$, respectively, thus, they can be computed by
$$
I\pqty{\mathbf{x;y}}=\int \mathcal{X}\times\mathcal{Y}p\pqty{x,y}\log\frac{p(x,y)}{p(x)p(y)}dxdy=\mathbb{E}_{p(x,y)}\bqty{\log\frac{p(x,y)}{p(y)p(x)}}
$$
where $\pqty{\mathbf{x,y}}$ is a pair of random variables with values over the space $\mathcal{X}\times\mathcal{Y}$, and $p(x)$, $p(y)$ are marginal probability of sending symbol $\mathbf{x}$ and received $\mathbf{y}$, while $p(x,y)$ is the joint probability of $\mathbf{x}$ and $\mathbf{y}$. The mutual information is equal to the Kullback-Leibler (KL) divergence between marginal probability and the joint probability, given by
$$
I(\mathbf{x};\mathbf{y})=D_{\mathrm{KL}}\pqty{p(x,y)|| p(x)p(y)}
$$
We also have the following theorem considering [31]

​	*Theorem 1*: The KL divergence admits the dual representation as follows
$$
D_{\mathrm{KL}}(P||Q)=\text{sup}_{T:\Omega\rightarrow R}E_P[T]-\log\pqty{E_Q[e^T]}
$$
where supremum is taken over all function $T$ so that the two expectations are finite.

​	According to the above theorem, the KL divergence can also be represented as
$$
D_{\mathrm{KL}}\pqty{p(x,y)||p(x)p(y)}\geq\mathbb{E}_{p(x,y)}\bqty{T}-\log\pqty{\mathbb{E}_{p(x)p(y)}\bqty{e^T}}
$$
​	Obviously, the lower bound of $I\pqty{\mathbf{x,y}}$ can be derived from
$$
I\pqty{\mathbf{x,y}}=D_{\mathrm{KL}}\pqty{p\pqty{x,y}||p(x)p(y)}\geq\mathbb{E}_{p(x,y)}\bqty{T}-\log\pqty{\mathbb{E}_{p(x)p(y)}\bqty{e^T}}
$$
​	However, to find a tight bound on the $I(\mathbf{x,y})$, we train the network $T$ unsupervised, it will be approximated by neural network. By the way, the expectation in above equation can be calculated by sampling, more samples, more converge. Then, we optimize the encoder by maximizing the mutual information and the related loss function can be given by
$$
\mathcal{L}_{\mathrm{MI}}\pqty{\mathbf{x,y};T}=\mathbb{E}_{p(x,y)}\bqty{f_T}-\log\pqty{\mathbb{E}_{p(x)p(y)}\bqty{e^{fT}}}
$$
where $f_T$ is composed by neural network, where the inputs are samples from $p(x,y)$, $p(x)$ and $p(y)$. In our scheme, $\mathbf{x}$ is generated by the $C_\alpha$ and $S_\beta$, the loss will be expressed as
$$
\mathcal{L}_\mathrm{MI}(\mathbf{x,y};T,\alpha,\beta)\leq I(\mathbf{x;y})
$$
thus, the loss can be computed and used to train the network weights to get $\alpha$, $\beta$ and $T$. Once the $\alpha$ and $\beta$ fixed, the mutual information can be estimated by training $T$. From the aspect of training, we can optimize the encoder by training $\alpha$ and $\beta$, when we obtain mutual information.

​	Performance criteria is another significant problem. The traditional system usually assess the performance through BER, which do not reflect the text information precisely. In the area of translation, BLEU is often used to evaluate the similarity between sentences, which is borrowed as one of our metric to compare the difference between words in two sentences in this paper, as well. Furthermore, to assess the similarity of semantic information, we also establish a new metric, to describe the similarity of two sentences from the aspects of semantics.

​	***BLEU-score***: By counting the difference of $n$-grams between transmitted and received texts. To illustrate the meaning of $n$-grams, 1-grams means a single words, while n-grams means block consist of n words.

​	Typically, we can simply define the BLEU as follows
$$
\log \text{BLEU}=\min\pqty{1-\frac{l_{\hat{\mathbf{s}}}}{l_{\mathbf{s}}},0}+\sum_{n=1}^N u_n\log p_n
$$
where $l_\mathbf{s}$ and $l_\hat{\mathbf{s}}$ are the length of the transmitted and received sentences $\mathbf{s}$ and $\hat{\mathbf{s}}$, respectively. and the $u_n$ is the weights of $n$-grams and $p_n$ is the $n$-grams score, which is expressed as
$$
p_n=\frac{\sum_k\min\pqty{C_k\pqty{\hat{\mathbf{s}}},C_k\pqty{s}}}{\sum_k\min\pqty{C_k\pqty{\hat{\mathbf{s}}}}}
$$
where $C_k\pqty{\cdot}$ is the frequency count function for $k$-th elements in $n$-th grams.

​	Finally, the range of the BLEU output is $\pqty{0,1}$, which indicates the degree of similarity between the decoded text and the original text. However, there is another trouble, even if the system restore the text concretely, the score is still not 1, sometimes, due to the existence of synonym. This paper beat the trouble by proposing a new metric, introducing the ***Sentence Similarity*** to the BLEU score system.

​	***Sentence Similarity***: Same words may have different meanings, traditional method, like *word2vec* cannot recognize the polysemy, but solve it by various vectors for same words. As our ***Semantic Similarity***, we intend to compute the sentence similarity between original $\mathbf{s}$ and recovered $\hat{\mathbf{s}}$, as follows
$$
\text{match}\pqty{\hat{\mathbf{s}},\mathbf{s}}=\frac{B_{\Phi}\pqty{\mathbf{s}}\cdot B_{\Phi}\pqty{\hat{\mathbf{s}}}^T}{\Vert B_{\Phi}\pqty{\mathbf{s}}\Vert\Vert B_{\Phi}\pqty{\hat{\mathbf{s}}}\Vert}
$$
where $B_{\Phi}$ means BERT, a pre-trained model with billions of parameters for extracting the semantic information, the similarity is the above equation with the value range of $\pqty{0,1}$, indicates the degree of similarity, the larger, the more similar.

​	Conclusively, BERT has been fed by billions of sentences, thus, this pre-trained model has learned semantics from countless sentence, able to generate vectors for different contexts. There, we express the transmitted semantics of sentence as $\mathbf{c}$, accordingly, the estimated semantics by receivers as $\hat{\mathbf{c}}$, which can make us assess the similarity by $\text{match}\pqty{\mathbf{c},\hat{\mathbf{c}}}$.

###### · Proposed Deep Semantic Communication Systems

​	In this section, we proposed a DNN to optimize the considered semantic Com system, named DeepSC. In our framework, Transformer is adopted for text understanding, while transfer learning for better ability of generality, in case to adapt dynamic and complicated communication environment.

​	Our framework is as shown in Fig. 2, consist of a transmitter, a channel interpreter and a receiver, among the transmitter is composed of a semantic encoder with a multiple Transformer for extracting semantics and a channel encoder with dense layer with different units for generating symbols, while the receiver is equipped with the corresponding part, i.e. channel decoder for symbol detection and semantic decoder for semantics estimation. Based on this framework, the loss function can be expressed as
$$
\mathcal{L}_\text{total}=\mathcal{L}_\mathrm{CE}\pqty{\mathbf{s},\hat{\mathbf{s}},\alpha,\beta,\chi,\delta}-\lambda\mathcal{L}_\text{MI}\pqty{\mathbf{x},\mathbf{y},T\alpha,\beta}
$$
where the first term is the loss representing the sentence similarity, aiming to minimize the difference between $\mathbf{s}$ and $\hat{\mathbf{s}}$, while the second term is the loss for mutual information, aiming to maximize the achieved data rate, and the parameter $\lambda$ is the weight for the latter.

​	What need to notice is the Transformer is well-known for its multi-head self-attention mechanism, which means it can predict the next words by viewing the predicted words, also means able to capture the distant dependency of the pronoun, demonstrating it can learn the semantics.

​	Owing to the different loss, we train our DeepSC with two phases. Having initialized the weight and bias $\mathbf{W,b}$ and represented the words by vectors, **Phase Ⅰ** is to train the mutual information model by unsupervised learning to estimate the data rate for next step. **Phase Ⅱ** is to train the whole system with $\mathcal{L}_\text{total}$ as loss. In this process, each phase aims to minimize the loss until criterion met, then the wanted max of iteration reached, or no longer decrease made by any term of $\mathcal{L}_\text{total}$. Our joint semantic-channel scheme will still preserve semantic information well, despite coder deal with bits rather than semantics like those separate schemes. More details about two-phases-training is as shown in Fig. 4 and following.

**Phase Ⅰ**: This phase estimate the mutual information by transmitted $\mathbf{X}$  and received $\mathbf{Y}$. Firstly, based on the KB, the set $\mathcal{K}$ generates a mini-batch of sentences $\mathbf{S}\in\mathfrak{R}^{B\times L \times 1}$, where $B$ means batch size and $L$ means the length of sentences. Secondly, passed the embedding layer the sentences are transformed into a dense word vector $\mathbf{E}\in\mathfrak{R}^{B\times L\times E}$, where $E$ is the dimension of vector. Thirdly, pass the semantic encoder to obtain $\mathbf{M}\in\mathfrak{R}^{B\times L \times V}$ as the extracted semantics, where the $V$ is the dimension of encoder output. Fourthly, the $\mathbf{M}$ is encoded into $\mathbf{X}$ to transmit in physical channel, where $\mathbf{X}\in\mathfrak{R}^{B\times NL\times 2}$. Fifthly, having passed physical channel, the receiver obtains distorted signal $\mathbf{Y}$. Sixthly, considering the equation of $\mathcal{L}_\text{MI}$, we can compute the mutual information with known $\mathbf{X}$ and $\mathbf{Y}$, under existent AWGN. Eventually, given the $\mathcal{L}_\text{MI}$, the stochastic gradient descent (SGD) is adopted to optimize the weights and bias of $f_T\pqty{\cdot}$.

**Phase Ⅱ**: This phase train the whole framework by cross-entropy and mutual information. Also like the above process of information stream transmission, generates sentences $\mathbf{S}$ from $\mathcal{K}$, from $\mathbf{S}$ to semantic $\mathbf{M}$, from $\mathbf{M}$ to transmitted symbols $\mathbf{X}$, from $\mathbf{X}$ to distorted received $\mathbf{Y}$ through channel, and restored the estimated semantic $\hat{\mathbf{M}}$. Finally, this whole network will be optimized by SGD with loss $\mathcal{L}_\text{total}$.

**Transfer Learning for Dynamic Environment**: Different communication scenario means different channel and data, re-train for each scenario is valueless, so the transfer learning is adopted to solve these different but related scenario. Similarly to above training process, the training modules contains mutual information estimation training and whole network training. However, in practice, we only need to download the pre-trained model and redesign and train semantic encoder/decoder with fixed channel encoder/decoder layer for different KB, while redesign and train the channel encoder/decoder with fixed semantic encoder/decoder. Even if confront the totally different scenario, our pre-trained model will assistant reduce the time of training. With the other modules trained, few epochs training will make us reach the global optimum.

###### · Numerical Results

###### · Citation List

**Deep learning based communication over the air**

Model-free training of end-to-end communication systems

**Deep learning based end-to-end wireless communication systems with conditional GAN as unknown channel**

**Model-driven deep learning for MIMO detection**

Model-free training of end-to-end communication systems

The semantic communication game

Deep joint source-channel coding for wireless image transmission



##### · *Semantic MIMO Systems for  Speech-to-Text Transmission*

###### · Abstract

​	This paper propose a semantic-aware speech-to-text transmission system for MIMO system, named SAC-ST. Particularly, we first design a speech-to-text semantic communication system for transmitter devices, which compress semantic information and generates the low-dimensional semantic features by Transformer. Moreover, by identifying the critical semantics and guaranteeing the accurate recovery, high semantic fidelity realized in our framework. Meanwhile, a DNN-based channel estimation model reduce the dependence on an accurate CSI, enforcing the feasibility and robustness in various circumstance. Over MIMO system, our proposed scheme is better than the ones without semantic-aware network, especially in low SINR. DNN-based channel estimation is also comparable to the perfect CSI.

###### · Introducion

​	Semantic Com is the second level of communication, delivering the semantics and optimizing the system in semantic level. To tackle with less bandwidth and higher demand of data, Semantic Com is regarded as the potential next communication technology. The conventional communication, the first level of communication, convert input message into bits, aiming to achieve a low BER by proper coding algorithm. Due to the bit era thrived past, the semantic level has been ignored. The aims of Semantic Com is to convey the message with inherent semantic information of itself. However, without unified theory to represent semantics from different resources and accurate mathematical formula to quantify the capacity, the research made less progress. Thanks to AI and its appliance on other part of communication system, it breaks the constraints of mathematical model.

​	DL-based Semantic Com maps the source message into a low-dimensional semantic features, reducing the data volume and mitigating the traffic without performance degradation. Thus, a Transformer-based joint semantic-channel codecs mechanism proposed in [9], named DeepSC, minimizing the semantic error by measuring the difference between input and recovery. Then, Peng proposed R-DeepSC to improve the fidelity and against the impairments inherent in the sentence by a error corrector in [10]. Jiang proposed a text Semantic Com system to ensure the reliability by hybrid automatic repeat request (HARQ) in [11]. Liang studied a reasoning-based paradigm to represent semantics by a graph-based knowledge architecture and semantic interpretation in [12]. Xie researched on task-oriented Semantic Com for text transmission in [13] and introduced a memory component to tackle with problems in dynamic environment. Nam developed a novel sequential system for text-to-image generation through multi-model generation in [15]. Weng designed a DeepSC-S by extracting and transmitting global semantic information for speech signal in [16]. Wei devised a real-time audio Semantic Com system to strengthen audio quality and alleviate physical channel attenuation by capturing the long-distance dependency in [17]. Han proposed a speech-oriented system with a redundancy removal module in [18]. Xiao designed a framework to reconstruct speech by a flexible rate-distortion trade-off in [19]. Weng proposed a task-oriented system, named DeepSC-ST in [20]. Zhang historically concluded a unified system for dimension controlling and data reducing with shared trainable parameters among various tasks in [21].





