#### Notes

$$
\gamma_{k}=\frac{\left|\int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{g}(\mathbf{s}) d \mathbf{s}\right|^{2}}{\sum_{i=1, i \neq g}^{G}\left|\int_{\mathcal{S}_{\mathrm{T}}} H_{k}(\mathbf{s}) J_{i}(\mathbf{s}) d \mathbf{s}\right|^{2}+\sigma_{k}^{2}}
$$

##### ·Without RIS:

$$
{H}_k(s)=\hat{\pmb{u}_k^T}\pmb{G}(r_k,s)\hat{\pmb{u}_y}
$$
​	Where $\mathbf{G}(r,s)$ satisfy:
$$
\mathbf{G}(\mathbf{r}, \mathbf{s})=-\frac{j \eta e^{-j \frac{2 \pi}{\lambda}\|\mathbf{r}-\mathbf{s}\|}}{2 \lambda\|\mathbf{r}-\mathbf{s}\|}\left(\mathbf{I}_{3}-\frac{(\mathbf{r}-\mathbf{s})(\mathbf{r}-\mathbf{s})^{T}}{\|\mathbf{r}-\mathbf{s}\|^{2}}\right)
$$

##### ·With RIS:

​	What we need to note is the multicast of RIS-aided MIMO is discreet so the transmission function can be derived into a matrix, while our CAPA multicast is continuous so we consider this condition according to the direct channel and split it into two parts . 


$$
H_k(s)=\hat{\pmb{u}_k^T}\pmb{G}_{d,k}(r_k,s)\hat{\pmb{u}_y}+\hat{\pmb{u}_k^T}\hat{\pmb{u}}_y
$$

$$
y_{k}=\left(\mathbf{h}_{\mathrm{r}, k}^{\mathrm{H}} \mathbf{E} \mathbf{H}_{\mathrm{dr}}+\mathbf{h}_{\mathrm{d}, k}^{\mathrm{H}}\right) \sum_{g=1}^{G} \mathbf{f}_{g} s_{g}+n_{k}
$$

##### ·System Model

​	We study a CAPA-based downlink multi-group multi-cast communication system, where a CAPA transmitter is equipped at BS to serve $K$ single antenna users. We model this system with reference to Wang's paper.

###### Transmit Signal

​	Simply know or define the following notation:

1. Surface area is $A_T=|\mathcal{S}_T|$
2. $\mathbf{J}(s, \omega)$ represents the Fourier Transform of the source current density vector at $\forall s\in \mathcal{S}_T$, where $\omega=2\pi/\lambda=2\pi f/c$, $f$ is the frequency, $\lambda$ is the wave length of signal.

​		ps: in our condition only a narrowband single-carrier and the current density vector of y-ax will be considered, simplified formula: $J(s)=\sum_{g=1}^GJ_g(s)x_g$ is the strength of source current density, where $J_g(s)$ is the source current density at point $s$ and $x_g$ is communication symbol of each multi-group, while the vector of current density is $\mathbf{J}(s)=J(s)\hat{\mathbf{u}}_y$

###### EM Channel and Receive Signal

3. $\mathbf{E}_k=\int_{\mathcal{S}_T}\mathbf{G}(\mathbf{r}_k,s)\mathbf{J}(s)ds$, represent the vector of the electric field at location $\mathbf{r}_k$ for each user $k$, contains both direction and strength of EF.
4. $\mathbf{G}(\mathbf{r}, \mathbf{s})=-\frac{j \eta e^{-j \frac{2 \pi}{\lambda}\|\mathbf{r}-\mathbf{s}\|}}{2 \lambda\|\mathbf{r}-\mathbf{s}\|}\left(\mathbf{I}_{3}-\frac{(\mathbf{r}-\mathbf{s})(\mathbf{r}-\mathbf{s})^{T}}{\|\mathbf{r}-\mathbf{s}\|^{2}}\right)$ is perceived as the transmission function for each user.
5. $E_k=\hat{\mathbf{u}}_k^T\mathbf{E}_k+n_k$ represents the eventual EF strength that user $k$ received in noisy condition, with considering the single-antenna vector direction $\hat{\mathbf{u}_k^T}$
6. Eventually, the transmission function for each user can be expressed as $H_k(s)=\hat{\mathbf{u}}_k^T\mathbf{G(r_k,s)}\hat{\mathbf{u}}_y$。
7. Rewrite $E_k$ as $E_k=\int_{\mathcal{S}_T}H_k(s)\sum_{g=1}^GJ_g(s)x_gds+n_k$
8. Mentioned that three parts
   1. Signal user $k$ in group $g$ want $\int_{\mathcal{S}_T}H_k(s)J_g(s)x_gds$
   2. Signal unexpected $\sum_{i=1,i\neq g}^G\int_{\mathcal{S}_T}H_kJ_g(s)x_gds$
   3. Noise $n_k$

###### Achievable Communication Rate

9. $P_k=\mathbb{E}\{\varepsilon_k \times \frac{|E_k|^2}{2} \times \frac{1}{\eta}\}=\mathbb{E}\{\frac{\varepsilon_k}{2\eta}|E_k|^2\}$ means the power user $k$ received, where $\varepsilon_k$ means absorption efficiency and $\eta$ as intrinsic impedance.
10. Also can be divided into three parts (due to satisfy $\mathbb{E}\{\mathbf{xx}^H\}=\mathbf{I}_K$, among $\mathbf{x}=[x_1,...,x_g]^T$):
    1. Signal user $k$ want $\frac{\varepsilon_k}{2\eta}|\int_{\mathcal{S}_T}H_k(s)J_g(s)ds|^2$
    1. Signal unexpected $\frac{\varepsilon_k}{2\eta}\sum_{i=1,i\neq g}^G|\int_{\mathcal{S}_T}H_k(s)J_i(s)ds|^2$
    1. Noise $\frac{\varepsilon_k}{2\eta}\sigma_k^2$
11. Thus, the SINR of each user $k$ in group $g$ will be expressed as: $\gamma_k=\frac{Signal_{expected}}{Signal_{unexpected}+Noise}$, also directly expressed as $\gamma_k=\frac{|\int_{\mathcal{S}_T}H_k(s)J_g(s)ds|^2}{\sum_{i=1,i\neq g}^G|\int_{\mathcal{S}_T}H_k(s)J_i(s)ds|^2+\sigma_k^2}$.

12. Power limit: $\sum_{g=1}^G\int_{\mathcal{S}_T}|J_g(s)|^2ds=P_T$.

#### ·Optimization Problem Analysis

$$
\begin{flalign}
\mathrm{J}_g(\mathbf{r})&=\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})\\
\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})&=\mathrm{H}_k(\mathbf{r})\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})\\
\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}&=\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})d\mathbf{r}\\
&=\sum_{i=1}^K\alpha_{i,g}\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}\\
&=\sum_{i=1}^K\alpha_{i,g}h_{k,i}
\end{flalign}
$$

​	We denote $h_{k,i}$ as $\int_\mathcal{A}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}$, and according to Gauss-Legendre also can be calculated by
$$
\begin{align}
h_{k,i/i,k}&=\int_\mathcal{S}\mathrm{H}^*_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}\\
&=\int_{-\frac{L_y}{2}}^{\frac{L_y}{2}}\int_{-\frac{L_x}{2}}^{\frac{L_x}{2}}\mathrm{H}^*_k(x)\mathrm{H}_i(y)dxdy\\
&\approx\frac{L_xL_y}{4}\sum_{m_y=1}^M\sum_{m_x=1}^M\omega_{m_x}\omega_{m_y}\mathrm{H}^*_k(\frac{\theta_{m_x}L_x}{2},\frac{\theta_{m_y}L_y}{2})\times\mathrm{H}_i(\frac{\theta_{m_x}L_x}{2},\frac{\theta_{m_y}L_y}{2})
\end{align}
$$
​	Where $\theta_m$ are roots of the $M$-th Legendre polynomial and $\omega_m$ are the quadrature weights. And when $M=10$ is the optimal choice to balance the complexity and accuracy.

​	Thus, excited electric field is given by
$$
\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}=\sum_{i=1}^K\alpha_{ig}h_{k,i}
$$
​	The $\gamma_k$ is thus given by
$$
\begin{flalign}
\gamma_k
&\approx\frac{|\mathcal{S}_k|\cdot|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'=1,g'\neq g}^G|\mathcal{S}_k|\cdot|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_{g'}(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}\\
&=\frac{|\mathcal{S}_k|\cdot|\sum_{i=1}^K\alpha_{ig}\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'=1,g'\neq g}^G|\mathcal{S}_k|\cdot|\sum_{i=1}^K\alpha_{ig'}\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}\\
&=\frac{|\mathcal{S}_k|\cdot|\sum_{i=1}^K\alpha_{ig}h_{k,i}|^2}{\sum_{g'=1,g'\neq g}^G|\mathcal{S}_k|\cdot|\sum_{i=1}^K\alpha_{ig'}h_{k,i}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}h^*_{i,j/j,i}}
\end{flalign}
$$
​	Consequently, we define $\mathbf{H}$ is a matrix constructed by $h_{ki}=\int_\mathcal{A}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}$, and $\mathbf{A}$ is a matrix constructed by $\alpha_{ik}$
$$
\left(
\begin{array}{ccc}
h_{11} & h_{12}&...&h_{1k}\\

\end{array}
\right)
$$

$$
\begin{pmatrix}
h_{11} & h_{12} & \cdots & h_{1k} \\
h_{21} & h_{22} & \cdots & h_{2k} \\
\vdots & \vdots & & \vdots \\
h_{K1} & h_{K2} & \cdots & h_{KK}
\end{pmatrix}
$$

$$
\begin{pmatrix}
\alpha_{11} & \alpha_{12} & \cdots &\alpha_{1G}\\
\alpha_{21} & \alpha_{22} & \cdots &\alpha_{2G}\\
\vdots &\vdots & & \vdots\\
\alpha_{K1}&\alpha_{K2}&\cdots&\alpha_{KG}
\end{pmatrix}
$$

$$
\begin{pmatrix}
e_{11}&e_{12}&\cdots&e_{1G}\\
e_{21}&e_{22}&\cdots&e_{2G}\\
\vdots&\vdots&&\vdots\\
e_{K1}&e_{K2}&\cdots&e_{KG}
\end{pmatrix}
$$

​	To better optimize the process of calculating integrals, simplify the complexity, we try to eliminate the cost of calculating the integrals of power.
$$
\begin{align}
\mathrm{J}_k(\mathbf{r})&=A'_k\mathrm{H}_k^*(\mathbf{r})-\sum_{i=1}^Kw_{k,i}B_i'\mathrm{H}_i^*(\mathbf{r})\\
\\
P_T
&=\sum_{k=1}^K\int_\mathcal{S}|\mathrm{J}_k(\mathbf{r})|^2d\mathbf{r}=\sum_{k=1}^K\int_\mathcal{S}|A'_k\mathrm{H}_k^*(\mathbf{r})-\sum_{i=1}^Kw_{k,i}B_i'\mathrm{H}_i^*(\mathbf{r})|^2d\mathbf{r}\\
&=\sum_{k=1}^K|A'_k|^2h_{k,k}-\sum_{k=1}^K\sum_{k=i}^K2\Re{w_{k,i}^*A'_kB'^*_ih_{k,i}}
+\sum_{k=1}^K\sum_{i=1}^K\sum_{j=1}^Kw_{k,i}w_{k,j}^*B'_iB'^*_jh_{i,j}
\end{align}
$$
​	Consider our $\mathrm{J}_k(\mathbf{r})$ and $\mathrm{H}_k(\mathbf{r})$ are both a complex scalar, so the power constraint can be expressed as:
$$
\begin{align}
P_T
&=\int_\mathcal{S}\sum_{g=1}^G|\mathrm{J}_g(\mathbf{r})|^2d\mathbf{r}=\sum_{g=1}^G\int_\mathcal{S}|\mathrm{J}_g(\mathbf{r})|^2d\mathbf{r}\\
&=\sum_{g=1}^G\int_\mathcal{S}\mathrm{J}_g(\mathbf{r})\mathrm{J}^*_g(\mathbf{r})d\mathbf{r}\\
&=\sum_{g=1}^G\int_\mathcal{S}\left(\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})\right)\left(\sum_{j=1}^K\alpha^*_{j,g}\mathrm{H}^*_j(\mathbf{r})\right)\\
&=\sum_{g=1}^G\int_\mathcal{S}\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}\mathrm{H}_i(\mathbf{r})\mathrm{H}^*_j(\mathbf{r})d\mathbf{r}\\
&=\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}\int_\mathcal{S}\mathrm{H}_i(\mathbf{r})\mathrm{H}^*_j(\mathbf{r})d\mathbf{r}\\
&=\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}h^*_{i,j/j,i}
\end{align}
$$
​	To better deal with those massive integrals in optimization, we try to simplify the calculating process, thus, the concrete approach are divided into two parts:
$$
\begin{flalign}
\mathrm{J}_g(\mathbf{r})&=\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})\tag{1a}\\
\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})&=\mathrm{H}_k(\mathbf{r})\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})\\
\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}&=\sum_{i=1}^K\alpha_{i,g}\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}\\
&=\sum_{i=1}^K\alpha_{i,g}h_{k,i}\label{EF}\tag{1b}
\end{flalign}
$$
​	Where $h_{k,i}\triangleq\int_\mathcal{A}\mathrm{H}_k(\mathbf{r})\mathrm{H}_i(\mathbf{r})d\mathbf{r}$ to simplify the integrals, we can easily gain the $E_k$ after knowing coefficient $\alpha_{i,g}$, then we continue to cope with the power constraint
$$
\begin{align}
P_T
&=\int_\mathcal{S}\sum_{g=1}^G|\mathrm{J}_g(\mathbf{r})|^2d\mathbf{r}=\sum_{g=1}^G\int_\mathcal{S}|\mathrm{J}_g(\mathbf{r})|^2d\mathbf{r}\tag{2a}\\
&=\sum_{g=1}^G\int_\mathcal{S}\mathrm{J}_g(\mathbf{r})\mathrm{J}^*_g(\mathbf{r})d\mathbf{r}\\
&=\sum_{g=1}^G\int_\mathcal{S}\left(\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i(\mathbf{r})\right)\left(\sum_{j=1}^K\alpha^*_{j,g}\mathrm{H}^*_j(\mathbf{r})\right)\\
&=\sum_{g=1}^G\int_\mathcal{S}\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}\mathrm{H}_i(\mathbf{r})\mathrm{H}^*_j(\mathbf{r})d\mathbf{r}\\
&=\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}\int_\mathcal{S}\mathrm{H}_i(\mathbf{r})\mathrm{H}^*_j(\mathbf{r})d\mathbf{r}\\
&=\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}h^*_{i,j/j,i}\label{power_constraint}\tag{2b}
\end{align}
$$
​	Considering the SINR  $\gamma_k$ and combining the equations $(\ref{EF})$ and $(\ref{power_constraint})$, we can eventually reorganize the SINR as following:
$$
\begin{flalign}
\gamma_k
&=\frac{|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'=1,g'\neq g}^G|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_{g'}(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}\\

&=\frac{|\sum_{i=1}^K\alpha_{ig}h_{k,i}|^2}{\sum_{g'=1,g'\neq g}^G|\sum_{i=1}^K\alpha_{ig'}h_{k,i}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}h^*_{i,j/j,i}}
\end{flalign}
$$
​	Where $h_{i,k}$ or $h^*_{i,k}$ can be calculated by Gauss-Legendre at first
$$
h_{k,i/i,k}\approx\frac{L_xL_y}{4}\sum_{m_y=1}^M\sum_{m_x=1}^M\omega_{m_x}\omega_{m_y}\mathrm{H}_k(\frac{\theta_{m_x}L_x}{2},\frac{\theta_{m_y}L_y}{2})\times\mathrm{H}_i(\frac{\theta_{m_x}L_x}{2},\frac{\theta_{m_y}L_y}{2})
$$
​	Now, we can optimize the $\alpha$ matrix without caution, $\gamma_k$ is just a function of $\alpha$, so is $R_k$, while our objective function is 
$$
\begin{align}
R_k
&=\log_2(1+\gamma_k)\\
&=\log_2\left(1+\frac{|\sum_{i=1}^K\alpha_{ig}h_{k,i}|^2}{\sum_{g'=1,g'\neq g}^G|\sum_{i=1}^K\alpha_{ig'}h_{k,i}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}h^*_{i,j/j,i}}\right)\\
\max_{\{\mathrm{J}_g(\mathbf{r})\}_{g=1}^G}\ &\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\{\log_2(1+\gamma_k)\}\\
\end{align}
$$
​	Thus, the problem is transformed into following structure
$$
\max_{\{\alpha_{k,g}\}_{k=1,g=1}^{K,G}}\ \sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\{\log_2\left(1+\frac{|\sum_{i=1}^K\alpha_{ig}h_{k,i}|^2}{\sum_{g'=1,g'\neq g}^G|\sum_{i=1}^K\alpha_{ig'}h_{k,i}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha^*_{j,g}h^*_{i,j/j,i}}\right)\}
$$

#### ·Matrix $\mathbf{A}$ Optimization Schemes

##### ·SDR

**Theory**

The problem satisfies the following form
$$
\begin{align}
&\min_{\mathbf{x}\in\mathbb{R}^n}\mathbf{x^T Cx}\\
&\text{s.t.}\mathbf{x^T A_i x}\unrhd_\mathrm{i}\mathrm{b_i}, \mathrm{i=1,\cdots,m}
\end{align}
$$
Where matrix $\mathbf{A_i,C}$ are both symmetric matrix, and $\unrhd_\mathrm{i}$ means one of $\geq,\leq,=$ satisfying QCQP

*Note*: this problem is not a convex problem, due to s.t. is not a convex set

Then, we introduce a new variable $\mathbf{X=xx}^\mathrm{T}$, and transform the original problem into
$$
\begin{align}
&\min_{\mathbf{X}\in\mathbb{S}^n} \Tr{(\mathbf{CX})}\\
&\text{s.t.}\Tr{(\mathbf{A_i X})}\unrhd_\mathrm{i}\mathrm{b_i,\ \ i=1,\cdots,m}\\
&\mathbf{X}\succeq0,\rank{(\mathbf{X})}=1
\end{align}
$$
Obviously, $\mathbf{x^T Cx}$ can be referred to as a $1\times 1$ matrix, i.e., $\mathbf{x^T Cx}=\tr{(\mathbf{x^T Cx})}$, and thanks to the property $\tr(\mathbf{AB})=\tr(\mathbf{BA})$ of trace, we can rewrite it as $\mathbf{x^T Cx}=\tr(\mathbf{x^T Cx})=\tr(\mathbf{Cxx^T})=\tr(\mathbf{CX})$, similarly, the constraint. Then, due to the $\mathbf{X}$ is linear and originated from $\mathbf{xx^T}$, the $\rank(\mathbf{X})$ must be 1 and non-negative.

**This problem thus transferred into a convex problem**, just adding another constraint $\rank(\mathbf{X})=1$, which can be relaxed alternatively.

Practically, we need to find a corresponding $\mathbf{X}$ and calculate the $\mathbf{x}$ from $\mathbf{X}$, but due to the relaxation of $\rank(\mathbf{X})=1$, there is not always be a correct $\mathbf{X}$ and $\mathbf{x}$. 

**Practice**

Thus, we need to approximate it  by $\min_{\mathbf{x}}\norm{\mathbf{X-xx^T}}_\mathrm{F}^2$, i.e. try to find the closest $\mathbf{x}$ satisfying the $\mathbf{X}$ we get. We gain the $\tilde{\mathbf{x}}=\sqrt{\lambda}\mathbf{q}$, where $\lambda,\mathbf{q}$ are the max eigenvector and eigenvalue respectively. 

However, the $\tilde{\mathbf{x}}$ we gained also might be not the feasible solution, and we need to find the closest feasible solution.

**Application**

We apply SDR method in a HBF problem.

The original problem is
$$
\begin{align}
&\min_{\mathbf{F}_\text{BB}}\ \ \norm{\mathbf{F}_\text{opt}-\mathbf{F}_\text{RF}\mathbf{F}_\text{BB}}_\mathbf{F}^2\\
&\text{s.t.} \norm{\mathbf{F}_\text{BB}}_\mathbf{F}=\frac{\mathrm{N_{RF}^t N_s}}{\mathrm{N_t}}
\end{align}
$$
Then, the authors proposed a SDR method to solve this problem as following steps.
$$
\begin{align}
\norm{\mathbf{F}_\text{opt}-\mathbf{F}_\text{RF}\mathbf{F}_\text{BB}}_\mathbf{F}^2
&=\norm{\text{vec}(\mathbf{F}_\text{opt}-\mathbf{F}_\text{RF}\mathbf{F}_\text{BB})}_2^2\\
&=\norm{\text{vec}(\mathbf{F}_\text{opt})-\text{vec}(\mathbf{F}_\text{RF}\mathbf{F}_\text{BB})}_2^2\\
&=\norm{\text{vec}(\mathbf{F}_\text{opt})-(\mathbf{I}_\mathrm{N_\mathrm{s}}\otimes\mathbf{F}_\text{RF})\text{vec}(\mathbf{F}_\text{BB})}_2^2\\
\end{align}
$$
And let
$$
\begin{align}
\mathbf{f}&=\text{vec}(\mathbf{F}_\mathrm{opt})\\
\mathbf{b}&=\text{vec}(\mathbf{F}_\mathrm{BB})\\
\mathbf{E}&=\mathbf{I}_\mathrm{N_s}\otimes\mathbf{F}_\mathrm{RF}
\end{align}
$$
The problem is transformed into
$$
\min_\mathbf{b}\norm{\mathrm{t}\mathbf{f-Eb}}_2^2\\
\text{s.t.}
\begin{cases}
\norm{\mathbf{b}}_2^2=\frac{\mathrm{N_RF^t N_s}}{\mathrm{N_t}}\\
\mathrm{t^2}=1
\end{cases}
$$
where $\mathrm{t}$ is an auxiliary variable to better solve our question. If $\mathrm{t}\pm1$ then the $\mathbf{f}$ is the optimal solution.

Considering the $\mathbf{x^T Cx}$ form, the objective function can be derived as
$$
\norm{\mathrm{t}\mathbf{f-Eb}}_2^2=\left[\mathbf{b}^\mathrm{H}\ \ \mathrm{t}\right]
\begin{bmatrix}
\mathbf{E}^\mathrm{H}\mathbf{E}&-\mathbf{E}^\mathrm{H}\mathbf{f}\\
\mathbf{f}^\mathrm{H}\mathbf{E}&\mathbf{f}^\mathrm{H}\mathbf{f}
\end{bmatrix}
\begin{bmatrix}
\mathbf{b}\\
\mathrm{t}
\end{bmatrix}
$$
which satisfied the form, let $\mathbf{y}=\begin{bmatrix}\mathbf{b}\\\mathrm{t}\end{bmatrix}$ and $\mathbf{C}=\begin{bmatrix}\mathbf{E}^\mathrm{H}\mathbf{E}&-\mathbf{E}^\mathrm{H}\mathbf{f}\\
\mathbf{f}^\mathrm{H}\mathbf{E}&\mathbf{f}^\mathrm{H}\mathbf{f}\end{bmatrix}$, the above formulation then become $\mathbf{y}^\mathrm{H}\mathbf{Cy}$, among $\mathbf{C}$ is obviously a semidefinite matrix. Similarly, the constraints can be
$$
\begin{align}
\norm{\mathbf{b}}_2^2&=\begin{bmatrix}\mathbf{b}^\mathrm{H}&\mathrm{t}\end{bmatrix}\begin{bmatrix}\mathbf{I}_\mathrm{N_RF^t}&0\\\mathbf{0}&0\end{bmatrix}\begin{bmatrix}\mathbf{b}\\\mathrm{t}\end{bmatrix}=\frac{\mathrm{N_RF^t N_s}}{\mathrm{N_t}}\\
\mathrm{t^2}&=\begin{bmatrix}\mathbf{b}^\mathrm{H}&\mathrm{t}\end{bmatrix}\begin{bmatrix}0_\mathrm{N_RF^t N_s}&0\\0&1\end{bmatrix}\begin{bmatrix}\mathbf{b}\\\mathrm{t}\end{bmatrix}=1
\end{align}
$$
Then, $\norm{\mathbf{b}}_2^2=\mathbf{y}^\mathrm{H}\mathbf{A_1y}$ and $\mathrm{t^2}=\mathbf{y}^\mathrm{H}\mathbf{A_2y}$, thus the problem can be
$$
\begin{align}
&\min_{\mathbf{Y}}\ \ \Tr(\mathbf{CY})\\
&\text{s.t.}
\begin{cases}
\Tr(\mathbf{A_1Y})=\frac{\mathrm{N_RF^t N_s}}{\mathrm{N_t}}\\
\Tr(\mathbf{A_2Y})=1\\
\mathbf{Y}\succeq0,\rank{(\mathbf{Y})}=1
\end{cases}
\end{align}
$$
Obviously, we relax the last constraint and solve this convex problem by CVX in matlab.

*Note*: whether the $\mathrm{t^2}=1$ and $\norm{\mathbf{b}}_2^2=\frac{\mathrm{N_RF^t N_s}}{\mathrm{N_t}}$. If not just enforce the $\mathrm{t^2}=1$ and $\norm{\mathbf{b}}_2^2=\frac{\mathrm{N_RF^t N_s}}{\mathrm{N_t}}$ by normalization.

**Reference List**

Paper list from ***Multi-Group Multicast Beamforming: Optimal  Structure and Efficient Algorithms***. The SDR is used in [2,4,6,10,11], and SCA is used in [14,15(uni),16,17], and a method combined ZF with SCA is used in [18]. Among SDR, the [2,4] is randomization methods, and [29] is rank-reducing methods.

**How They Transform the SINR to a SDR problem?**

​	Actually, they do not transform the $R_k$ into a SDR but transform the SINR to a constraint, and the original power constraint become the objective function. This problem is a totally brand-new problem as following
$$
\begin{align}
&\min\ \sum_{g=1}^G|P_g|\\
\text{s.t.}&\gamma_k=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma_k^2}\geq\gamma_\text{threshold}
\end{align}
$$
​	In this problem, they try to minimize the power under the constraint of threshold SINR $\gamma_\text{threshold}$, reducing the cost within a appropriate communication rate/SINR.

​	Then, we continue our analysis on this problem and derive the core steps of SDR scheme, and try to rewrite our original problem in such a approach and respective.

​	Naturally, from aspect of objective function, the power of each group $g$ can be expressed as $\mathbf{a}_g^\mathrm{H}\mathbf{Ha}_g$, where $\mathbf{H}=\mathbf{H}^\mathrm{H}$, and the objective function will be $\min\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{Ha}$, and is thus derived as $\mathbf{a}_g^\mathrm{H}\mathbf{Ha}=\tr(\mathbf{a}_g^\mathrm{H}\mathbf{Ha})=\tr(\mathbf{Ha}\mathbf{a}_g^\mathrm{H})=\tr(\mathbf{HA}_g)$, among $\mathbf{A}_g\mathbf{=aa}^\mathrm{H}_g$, and objective function is $\min\sum_{g=1}^G\tr(\mathbf{HA}_g)$.

​	Retrospect the SINR constraint, we reorganize the inequation as following steps.
$$
\begin{align}
\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma_k^2}&\geq\gamma_\text{threshold}\\
\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2-\gamma_\text{threshold}\pqty{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma_k^2}&\geq0\\
\mathbf{a}_g^\mathrm{H}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g-\gamma_\text{threshold}\pqty{\sum_{g'\neq g}^G\mathbf{a}_{g'}^\mathrm{H}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}+\sigma_k^2}\geq0\\
\tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H})-\gamma_\text{threshold}\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}\mathbf{a}_{g'}^\mathrm{H}}+\sigma_k^2}\geq0
\end{align}
$$
​	If we define the $\mathbf{h}_k\mathbf{h}_k^\mathrm{H}$ as $\mathbf{Q}_k$, and apparently satisfying $\mathbf{Q}_k^\mathrm{H}=\mathbf{Q}_k$ and semidefinite, and we already define the $\mathbf{a}_g\mathbf{a}_g^\mathrm{H}$ as $\mathbf{A}_g$, the constraint is given by
$$
\tr(\mathbf{Q}_k\mathbf{A}_g)-\gamma_\text{threshold}\pqty{\sum_{g'\neq g}^G\tr(\mathbf{Q}_k\mathbf{A}_{g'})+\sigma_k^2}\geq0
$$
​	Reorganize the whole problem, we finally gain the following concave problem as 
$$
\begin{align}
&\min\sum_{g=1}^G\tr(\mathbf{HA}_g)\\
\text{s.t.}\ \ &\tr(\mathbf{Q}_k\mathbf{A}_g)-\gamma_\text{threshold}\pqty{\sum_{g'\neq g}^G\tr(\mathbf{Q}_k\mathbf{A}_{g'})+\sigma_k^2}\geq0,\forall k\in\mathcal{K}\\
&\mathbf{A}_g\succeq0,\ \ \rank(\mathbf{A}_g)=1,\ \ \forall g\in\mathcal{G}\\
\end{align}
$$

​	Which can be easily solve by the solver like CVX in matlab.

***Note***: The above scheme is only used in the problem within the SINR constraint and the objective function of minimizing the power cost, aiming to reduce the cost under a given rate or SINR. And we just transform this problem with the conventional technique. I do not even know why the problem is non-concave, while the transformed one is concave.

**How we cope with the $\log$ in our problem and formulate the SDR problem?**

​	In early research, find a trace representing the rate with $\log$ form baffles us hardly, so I have to transform my original problem to another similar problem, but now thanks to the technique of handling $\log$, I can transform my original problem to the one of trace form easier and gain the optimal solution with conventional SDR solver.

​	We can easily convert the power constraint of original problem into a trace form, but baffled by the rate $R_k$, then we recall the rate $R_k$ and try to complete the transforming.
$$
\begin{align}
R_k
&=\log_2\left(1+\frac{|\sum_{i=1}^K\alpha_{ig}h_{k,i}|^2}{\sum_{g'=1,g'\neq g}^G|\sum_{i=1}^K\alpha_{ig'}h_{k,i}|^2+\sigma_k^2}\right)\\
&=\log_2\pqty{1+\frac{\vqty{\mathbf{a}_g\mathbf{h}_k}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{a}_{g'}\mathbf{h}_k}^2+\sigma_k^2}}\\
&=\log_2\pqty{\frac{{\sum_{g=1}^G\vqty{\mathbf{a}_{g}\mathbf{h}_k}^2+\sigma_k^2}}{{\sum_{g'\neq g}^G\vqty{\mathbf{a}_{g'}\mathbf{h}_k}^2+\sigma_k^2}}}\\
&=\log_2\pqty{{\sum_{g=1}^G\vqty{\mathbf{a}_{g}\mathbf{h}_k}^2+\sigma_k^2}}-\log_2\pqty{{\sum_{g'\neq g}^G\vqty{\mathbf{a}_{g'}\mathbf{h}_k}^2+\sigma_k^2}}\\
&=\log_2\pqty{\sum_{g=1}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H})+\sigma_k^2}-\log_2\pqty{\sum_{g'\neq g}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}\mathbf{a}_{g'}^\mathrm{H})+\sigma_k^2}
\end{align}
$$
​	If we replace the $\log$ with two auxiliary variables $\tau_k,\varepsilon_k$, then the rate will be $R_k=\tau_k-\varepsilon_k$, which is obviously a convex objective function, and the variables is derived as
$$
\begin{align}
e^{\tau_k}&=\sum_{g=1}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H})+\sigma_k^2\\
e^{\varepsilon_k}&=\sum_{g'\neq g}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}\mathbf{a}_{g'}^\mathrm{H})+\sigma_k^2
\end{align}
$$
​	Reorganize the problem, we can reformulate it as
$$
\begin{align}
\max_{\qty{\tau_k,\varepsilon_k}_{k=1}^K,\qty{\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}_{g=1}^G}\ \ \ &\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\qty{\tau_k-\varepsilon_k}\\
\text{s.t.}\ \ \ &\sum_{g=1}^G\Tr(\mathbf{Ha}_g\mathbf{a}_g^\mathrm{H})\leq P_T\\
&e^{\tau_k}\leq\sum_{g=1}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H})+\sigma_k^2\\
&e^{\varepsilon_k^n}\pqty{\varepsilon_k-\varepsilon_k^n+1}\geq\sum_{g'\neq g}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}\mathbf{a}_{g'}^\mathrm{H})+\sigma_k^2\\
&\mathbf{a}_g\mathbf{a}_g^\mathrm{H}\succeq0, \ \ \text{Rank}(\mathbf{a}_g\mathbf{a}_g^\mathrm{H})=1\\
\end{align}
$$
​	Where $\varepsilon_k^n$ is the $n$ iteration and as the parameter of the first-order Taylor expansion of original $e^{\varepsilon_k}$.

​	Combined with my analysis, I intend to interpret the process of this algorithm in a more understandable way.

**Algorithm**:

1. The problem is the above problem we reorganize as the paper, which is regard as a convex problem
2. We have some fixed parameters like $P_T$ and calculation accuracy $\Delta_\text{set}$, and collaboratively optimize the variables $\qty{\tau_k^n,\varepsilon_k^n}_{k=1}^K$ and $\qty{\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}_{g=1}^G$, then finally derive the WSR $\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\qty{\tau_k-\varepsilon_k}$ to compare with the one of last time until satisfying the difference is less than we set $\Delta_\text{set}$.
3. Concretely, I think the process is a circle as following
   1. We set a original matrix $\mathbf{a}_g\mathbf{a}_g^\mathrm{H}$/$\mathbf{X}_g$ or gain from last iteration, then calculate the limit for $\tau_k,\varepsilon_k$;
   2. Having known the limit of $\tau_k,\varepsilon_k$, the main objective function can be easily optimize as a convex;
   3. Then, we gain again the $\tau_k$ and $\varepsilon_k$ from above convex formulation, but this time we calculated $\mathbf{X}_g$ preparing for next iteration, these $\mathbf{X}_g$ is the one mentioned in step 3.1
   4. What need to note that we still need to compare the WSR of this time with the one of last time, deciding whether or not continue our circle by the $\Delta_\text{set}$ and $\Delta^n$.
4. When the circle is ended, we can gain the $\tau_k^n,\varepsilon_k^n$ as the optimal solution for $\tau_k,\varepsilon_k$, we can calculate the $\mathbf{X}_g$ by 

$$
\begin{align}
e^{\tau_k}&=\sum_{g=1}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H})+\sigma_k^2\\
e^{\varepsilon_k}&=\sum_{g'\neq g}^G\Tr(\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}\mathbf{a}_{g'}^\mathrm{H})+\sigma_k^2
\end{align}
$$

5. Gain $\mathbf{a}_g$ from $\mathbf{X}_g$ or approximate it as the optimal $\mathbf{a}_g$.

**Our Problem:**

​	To further simplify the problem we reformulate, we introduce the SDR and Taylor expansion to convert our problem into a convex one, and eventually solve them through using SCA scheme.

​	Firstly, let rewrite the problem
$$
\begin{align}
\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{1+\gamma'_k}}
&=\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{1+\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}}}\\
&=\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{\frac{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}}}\\
&=\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}-\log_2\pqty{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}}\\
&=\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{\sum_{g=1}^G\Tr\pqty{\mathbf{a}_g^\mathrm{H}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}}-\log_2\pqty{\sum_{g'\neq g}^G\Tr\pqty{\mathbf{a}_g^\mathrm{H}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}}}\\
&=\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{\sum_{g=1}^G\Tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}}-\log_2\pqty{\sum_{g'\neq g}^G\Tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}}}\\
&=\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}}-\log_2\pqty{\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}}}\\
\end{align}
$$
​	Thus, the problem is converted into a convex one
$$
\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}}-\log_2\pqty{\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}}}\\
$$
​	If we denote the $\log$ form as following
$$
\begin{align}
\tau_k&=\log_2\pqty{\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}}\\
\varepsilon_k&=\log_2\pqty{\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}}
\end{align}
$$
​	The problem can be more simple as
$$
\max_{\qty{\mathbf{A}}}\sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\tau_k-\varepsilon_k}
$$
​	Then the problem combined with new constraint of $\tau_k$ and $\varepsilon_k$ will given by
$$
\begin{align}
\max_{\qty{\tau_k,\varepsilon_k}_{k=1}^K,\qty{\mathbf{X}_g}_{g=1}^G}\ \ \ &\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\qty{\tau_k-\varepsilon_k}\\
&e^{\tau_k}\leq\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}\\
&e^{\varepsilon_k^n}\pqty{\varepsilon_k-\varepsilon_k^n+1}\geq\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}\\
&\mathbf{X}_g\succeq0, \ \ \text{Rank}\pqty{\mathbf{X}_g}=1,\ \ \forall g\in\mathcal{G}\\
\end{align}
$$

**Multicast Communication**

​	However, the objective function $\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\qty{\tau_k-\varepsilon_k}$ is still non-convex, so we need to transform it in a similar way

​	First, we restate the original problem again
$$
\begin{flalign}
\max_{\{\mathrm{J}_g(\mathbf{r})\}_{g=1}^G}\ \ \ &\sum_{g=1}^G\alpha_g\min_{k\in\mathcal{K}_g}\{\log_2(1+\gamma'_k)\}
\end{flalign}
$$
​	Where 
$$
\gamma'_k=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}
$$
​	Then, the $R'_k$ can be derived as
$$
R'_k=\log_2\pqty{1+\gamma'_k}=\log_2\pqty{\frac{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}}
$$
​	To simplify the problem, we can express our problem as
$$
\max_{\qty{\varepsilon_g}_{g=1}^G}\sum_{g=1}^G\omega_g\varepsilon_g
$$
​	Where
$$
\begin{align}
e^{\varepsilon_g}&\leq \frac{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g},\forall k\in\mathcal{K}_g\\
0&\geq e^{\varepsilon_g}\pqty{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}-\pqty{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}\\
0&\geq e^{\varepsilon_g}\pqty{\sum_{g'\neq g}^G\mathbf{a}_{g'}^{\mathrm{H}}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}-\pqty{\sum_{g=1}^G\mathbf{a}_g^{\mathrm{H}}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}\\
0&\geq e^{\varepsilon_g}\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}\mathbf{a}_{g'}^{\mathrm{H}}}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\tr\pqty{\mathbf{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}}-\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^{\mathrm{H}}}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\tr\pqty{\mathbf{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}}\\
0&\geq e^{\varepsilon_g}\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\tr\pqty{\mathbf{H}\mathbf{X}_g}}-\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\tr\pqty{\mathbf{H}\mathbf{X}_g}}\\
\end{align}
$$
​	Thus, the problem can be transformed as
$$
\begin{align}
\max_{\qty{\varepsilon_g}_{g=1}^G}\ \ \  \ &\sum_{g=1}^G\omega_g\varepsilon_g\\
\text{s.t.}\ \  \ \ \ &\mathbf{X}_g\succeq0,\text{Rank}\pqty{\mathbf{X}_g}=1,\forall g\in\mathcal{G}\\
&0\geq\pqty{e^{\varepsilon_g}-1}\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\tr\pqty{\mathbf{HX}_g}}+\tr\pqty{\mathbf{Q}_k\mathbf{X}_g},\forall g\in\mathcal{G}
\end{align}
$$

$$
\varepsilon_g=\min_{k\in\mathcal{K}_g}\qty{\log_2\pqty{1+\gamma'_k}},\forall k\in\mathcal{K}_g
$$



Then
$$
\begin{align}
\gamma_k
&=\frac{\vqty{\int_\mathcal{S}\mathrm{H}_k(\mathbf{s})\sum_{i=1}^Ka_{i,g}\mathrm{H}_i(\mathbf{s})d\mathbf{s}}^2}{\sum_{g'\neq g}^G\vqty{\int_\mathcal{S}\mathrm{H}_k(\mathbf{s})\sum_{i=1}^Ka_{i,g'}\mathrm{H}_i(\mathbf{s})d\mathbf{s}}^2+\sigma_k^2}\\
&=\frac{\vqty{\sum_{i=1}^Ka_{i,g}h_{k,i}}^2}{\sum_{g'\neq g}^G\vqty{\sum_{i=1}^Ka_{i,g'}h_{k,i}}^2+\sigma_k^2}\\
&=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma^2},\\
\end{align}
$$
Then
$$
\begin{align}
\max_{\mathbf{A}}\ \ \ &\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\{\log_2\pqty{\frac{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2+\sigma^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma^2}}\}\label{problem_a},\\
\mathrm{s.t.}\ \ \ &\sum_{g=1}^G\tr\pqty{\mathbf{H}\mathbf{X}_g}=P_T,\label{problem_b}
\end{align}
$$
Then
$$
\begin{align}
e^{R_g}\leq \frac{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2},\forall k\in\mathcal{K}_g\\
\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}e^{R_g}\leq\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}
\end{align}
$$

However, the formulation is still non-concave, so we further exchange to another method to find a more appropriate transformation to gain a complete concave formulation as follows
$$
\begin{align}
R_g
&\leq\log_2\pqty{\frac{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}}\\
&=\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}\\
\end{align}
$$
where the $\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}$ is concave, but $-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}$ ($\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}$ is concave but which is not fit to its negative), so we try to find its first-order Taylor Expansion to replace it for a totally concave constraint for $R_g$.

Considering the property of trace, the derivative of $\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}$ will given by
$$
\frac{\mathbf{Q}_k}{\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}\ln2}
$$

Thereby, the first-order Taylor Expansion will be expressed as
$$
\begin{align}
f\pqty{\mathbf{X}_1,\mathbf{X}_2,\cdots,\mathbf{X}_G}&\approx f\pqty{\mathbf{X}^n_1,\mathbf{X}^n_2,\cdots,\mathbf{X}^n_G}+\sum_{j\neq g}^G\frac{\partial f}{\partial\mathbf{X}_j}\vert_{\mathbf{X}_{j}=\mathbf{X}_{j}^n}\pqty{\mathbf{X}_j-\mathbf{X}_j^n}\\
&=\log_2\pqty{\sum_{j\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_j^n}+\sigma_k^2}+\sum_{j\neq g}^G\frac{\mathbf{Q}_k:\pqty{\mathbf{X}_j-\mathbf{X}_j^n}}{\pqty{\sum_{j\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_j^n}+\sigma_k^2}\ln2}
\end{align}
$$
where $:$ means Frobenius inner product, and the $\mathbf{A}:\mathbf{B}$ with same dimension is calculated by $\mathbf{A}:\mathbf{B}=\sum_{i}^M\sum_{j}^Na_{ij}b_{ij}$.

Consequently, the problem is eventually transformed into 
$$
\begin{align}
\max_{\qty{\mathbf{X}_g,R_g}_{g=1}^G}\ \ \ \  \ &\sum_{g=1}^G\omega_gR_g\\
\text{s.t.}\ \ \ \ \ &R_g\leq\log_2\pqty{\sum_{j=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_j}+\sigma_k^2}-\log_2\pqty{\sum_{j\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_j^n}+\sigma_k^2}-\frac{\sum_{j\neq g}^G\tr\pqty{\mathbf{Q}_k\pqty{\mathbf{X}_j-\mathbf{X}_j^n}}}{\pqty{\sum_{j\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_j^n}+\sigma_k^2}\ln2},\forall k\in\mathcal{K}_g, g\in\mathcal{G}\\
&\sum_{g=1}^G\tr\pqty{\mathbf{HX}_g}=P_T\\
&\mathbf{X}_g\succeq0,\text{rank}\pqty{\mathbf{X}_g}=1,\forall g\in\mathcal{G}
\end{align}
$$


Having relax the rank-1 constraint the problem will become a totally convex problem, which is easily solved by CVX solvers.



Or we try another method as follows
$$
\begin{align}
R_g&=\tau_g-\varepsilon_g\\
\tau_g&\leq\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}\\
\varepsilon_g&\geq\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}
\end{align}
$$
Then
$$
\begin{align}
e^{\tau_g}&\leq\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2\\
e^{\varepsilon_g}&\geq\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2
\end{align}
$$
Introduce the first-order Taylor Expansion to expentional variable $\varepsilon_g$, we can eventually derive the constraint as
$$
\begin{align}
e^{\tau_g}&\leq\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2\\
e^{\varepsilon_g^n}\pqty{x-e^{\varepsilon_g^n}+1}&\geq\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2
\end{align}
$$




##### ·SOCP-Based Majorization-Minimization

##### ·DeepLearning

##### ·Wang's scheme







$$
\mathrm{Y}_k
=E_k+\mathrm{N}_k=\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{g=1}^G\mathrm{J}_g(\mathbf{r})s_gd\mathbf{r}+\mathrm{N}_k(\mathbf{s})
$$

#### Formula Reorganize

$$
\begin{align}
\max_{\qty{\mathbf{X}_g}_{g=1}^G}\ \ \ \ &\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}}\qty{\log_2\pqty{\frac{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}}},\\
\text{s.t.}\ \ \ \ 
&\sum_{g=1}^G\tr\pqty{\mathbf{HX}_g}=P_T,\\
&\mathbf{X}_g\succeq0,\rank\pqty{\mathbf{X}_g}=1,\forall g\in\mathcal{G}.
\end{align}
$$



#### Broadcast


$$
\begin{align}
R_g&\leq\log_2\pqty{\frac{P_k}{P_k-\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}}},\forall k\in \mathcal{K}_g\\
&\leq\log_2\pqty{P_k}-\log_2\pqty{P_k-\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}}
\end{align}
$$
where $P_k=\sum_{j=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_j}+\sigma_k^2$, is the total electric power at point $\mathbf{r}_k$, derived from total observed electric field $\mathrm{Y}_k$ at point $\mathbf{r}_k$.
$$
\begin{align}
\log_2\pqty{P_k-\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}}
\end{align}
$$



$$
\begin{align}
\gamma_k
&=\frac{\vqty{\int_\mathcal{S}\mathrm{J}\pqty{\mathbf{s}}\mathrm{H}_k\pqty{\mathbf{s}}d\mathbf{s}}^2}{\sigma_k^2}\\
&=\frac{\vqty{\int_\mathcal{S}\sum_{i=1}^K\alpha_i\mathrm{H}^*_i\pqty{\mathbf{s}}\mathrm{H}_k\pqty{\mathbf{s}}d\mathbf{s}}^2}{\sigma_k^2}\\
&=\frac{\vqty{\sum_{i=1}^K\alpha_i\int_\mathcal{S}\mathrm{H}^*_i\pqty{\mathbf{s}}\mathrm{H}_k\pqty{\mathbf{s}}d\mathbf{s}}^2}{\sigma_k^2}\\
&=\frac{\vqty{\sum_{i=1}^K\alpha_ih_{i,k}}^2}{\sigma_k^2}=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}}^2}{\sigma_k^2}\\
&=\frac{\mathbf{h}_k^\mathrm{H}\mathbf{a}\pqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}}^\mathrm{H}}{\sigma_k^2}=\frac{\mathbf{h}_k^\mathrm{H}\mathbf{aa}^\mathrm{H}\mathbf{h}_k}{\sigma_k^2}\\
&=\frac{\tr\pqty{\mathbf{h}_k^\mathrm{H}\mathbf{aa}^\mathrm{H}\mathbf{h}_k}}{\sigma_k^2}=\frac{\tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{aa}^\mathrm{H}}}{\sigma_k^2}\\
&=\frac{\tr\pqty{\mathbf{Q}_k\mathbf{X}}}{\sigma_k^2}
\end{align}
$$
where $\mathbf{Q}_k=\mathbf{h}_k\mathbf{h}_k^\mathrm{H}$ and $\mathbf{X}=\mathbf{aa}^\mathrm{H}$.

Thus, the ratio of communication will be given by
$$
R=\min\qty{\log_2\pqty{1+\gamma_k}}
$$
which can be further represented by
$$
R\leq\log_2\pqty{1+\gamma_k}=\log_2\pqty{1+\frac{\tr\pqty{\mathbf{Q}_k\mathbf{X}}}{\sigma_k^2}}
$$
and similarly, the constraint of power can also be derived from following steps
$$
\begin{align}
&\int_\mathcal{S}\vqty{\mathrm{J}\pqty{\mathbf{s}}}^2d\mathbf{s}\\
=&\int_\mathcal{S}\sum_{i=1}^K\sum_{j=1}^K\alpha_i\mathrm{H}_i\pqty{\mathbf{s}}\alpha_j\mathrm{H}_j^*\pqty{\mathbf{s}}d\mathbf{s}\\
=&\sum_{i=1}^K\sum_{j=1}^K\alpha_i\alpha_j\int_\mathcal{S}\mathrm{H}_i\pqty{\mathbf{s}}\mathrm{H}_j^*\pqty{\mathbf{s}}d\mathbf{s}\\
=&\sum_{i=1}^K\sum_{j=1}^K\alpha_i\alpha_jh_{i,j}=\mathbf{a}^\mathrm{H}\mathbf{H}\mathbf{a}
\end{align}
$$
Eventually, the optimization problem can be formulated as
$$
\begin{align}
\max_{\mathbf{X}}\ \ \ \ \ &R\\
\text{s.t}\ \ \ \ \ &R\leq\log_2\pqty{\frac{\tr\pqty{\mathbf{Q}_k\mathbf{X}}+\sigma_k^2}{\sigma_k^2}}\\
&\mathbf{a}^\mathrm{H}\mathbf{Ha}=P_T\\
&\mathbf{X}\succeq0,\rank\pqty{\mathbf{X}}=1
\end{align}
$$


Multicast
$$
\begin{align}
\gamma_k&=\frac{|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'\neq g}^G|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_{g'}(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}\\
&=\frac{\vqty{\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{i=1}^K\alpha_{i,g}\mathrm{H}_i\pqty{\mathbf{r}}d\mathbf{r}}^2}{\sum_{g'=1}^G\vqty{\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{i=1}^K\alpha_{i,g'}\mathrm{H}_i\pqty{\mathbf{r}}d\mathbf{r}}^2+\sigma_k^2}\\
&=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma_k^2}\\
R_k&=\log_2\pqty{1+\gamma_k}=\log_2\pqty{\frac{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2+\sigma_k^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma_k^2}}\\
R_g&=\min_{k\in\mathcal{K}_g}\qty{R_k}\\
\text{SumRate}&=\text{weights}*R_g*\text{Users/PerGroup}
\end{align}
$$
Unicast
$$
\begin{align}
\gamma_k&=\frac{|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_k(\mathbf{r})d\mathbf{r}|^2}{\sum_{k'\neq k}^K|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_{k'}(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}\\
&=\frac{\vqty{\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{i=1}^K\alpha_{i,k}\mathrm{H}_i\pqty{\mathbf{r}}d\mathbf{r}}^2}{\sum_{k'=1}^K\vqty{\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{i=1}^K\alpha_{i,k'}\mathrm{H}_i\pqty{\mathbf{r}}d\mathbf{r}}^2+\sigma_k^2}\\
&=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_k}^2}{\sum_{k'\neq k}^K\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{k'}}^2+\sigma_k^2}\\
R_k&=\log_2\pqty{1+\gamma_k}=\log_2\pqty{\frac{\sum_{k=1}^K\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_k}^2+\sigma_k^2}{\sum_{k'\neq k}^K\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{k'}}^2+\sigma_k^2}}\\
\text{SumRate}&=\text{weights}*R_k
\end{align}
$$
Broadcast
$$
\begin{align}
\gamma_k&=\frac{|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}(\mathbf{r})d\mathbf{r}|^2}{\sigma_k^2}\\
&=\frac{\vqty{\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{i=1}^K\alpha_{i}\mathrm{H}_i\pqty{\mathbf{r}}d\mathbf{r}}^2}{\sigma_k^2}\\
&=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}}^2}{\sigma_k^2}\\
R_k&=\log_2\pqty{1+\gamma_k}=\log_2\pqty{\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}}^2+\sigma_k^2}{\sigma_k^2}}\\
R_g&=\min_{k\in\mathcal{K}_g}\qty{R_k}\\
\text{SumRate}&=\text{weights}*R_g*\text{Users/PerGroup}
\end{align}
$$


Broadcast
$$
\begin{align}
R_k&=\log_2\pqty{\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}}^2+\sigma_k^2}{\sigma_k^2}}\\
&=\log_2\pqty{\frac{\tr\pqty{\mathbf{Q}_k\mathbf{X}}+\sigma_k^2}{\sigma_k^2}}\\
\int_\mathcal{S}\vqty{\mathrm{J}\pqty{\mathbf{s}}}^2d\mathbf{s}&=\int_\mathcal{S}\sum_{i=1}^K\sum_{j=1}^K\alpha_i\alpha_j\mathrm{H}_i(\mathbf{s})\mathrm{H}_j(\mathbf{s})d\mathbf{s}\\
&=\sum_{i=1}^K\sum_{j=1}^K\alpha_i\alpha_j\int_\mathcal{S}\mathrm{H}_i(\mathbf{s})\mathrm{H}_j(\mathbf{s})d\mathbf{s}\\
&=\mathbf{a}^\mathrm{H}\mathbf{Ha}\leq P_T
\end{align}
$$
Multicast
$$
\begin{align}
R_k&=\log_2\pqty{\frac{\sum_{g=1}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2+\sigma_k^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\sigma_k^2}}\\
&=\log_2\pqty{\frac{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}}\\
\sum_{g=1}^G\int_\mathcal{S}\vqty{\mathrm{J}_g\pqty{\mathbf{s}}}^2d\mathbf{s}&=\int_\mathcal{S}\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha_{j,g}\mathrm{H}_i(\mathbf{s})\mathrm{H}_j(\mathbf{s})d\mathbf{s}\\
&=\sum_{i=1}^K\sum_{j=1}^K\alpha_{i,g}\alpha_{j,g}\int_\mathcal{S}\mathrm{H}_i(\mathbf{s})\mathrm{H}_j(\mathbf{s})d\mathbf{s}\\
&=\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{Ha}_g\leq P_T
\end{align}
$$

$$
\begin{align}
h_{i,j}&=\int_\mathcal{S}\mathrm{H}_i(\mathbf{s})\mathrm{H}_j(\mathbf{s})d\mathbf{s}\\
\mathbf{h}_k&=\bqty{h_{k,1},h_{k,2},\cdots,h_{k,K}}^\mathrm{H}\\
\mathbf{H}&=\bqty{\mathbf{h}_1,\mathbf{h}_2,\cdots,\mathbf{h}_K}^\mathrm{H}\\

\mathbf{a}_g&=\bqty{a_{1,g},a_{2,g},\cdots,a_{K,g}}^\mathrm{H}
\end{align}
$$




Taylor
$$
\begin{align}
R_k&=\log_2\pqty{\frac{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}}\\
&=\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}
\end{align}
$$



$$
\begin{align}
&\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}\\
\approx&\frac{\sum_{g'\neq g}^G\mathbf{Q}_k:\pqty{\mathbf{X}_j-\mathbf{X}_j^\pqty{n}}}{\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}\ln2}+\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}
\end{align}
$$

$$
R_k\approx\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}-\frac{\sum_{g'\neq g}^G\mathbf{Q}_k:\pqty{\mathbf{X}_j-\mathbf{X}_j^\pqty{n}}}{\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}\ln2}-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}
$$

$$
R_g\leq\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}-\frac{\sum_{g'\neq g}^G\mathbf{Q}_k:\pqty{\mathbf{X}_j-\mathbf{X}_j^\pqty{n}}}{\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}\ln2}-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}
$$






$$
\begin{align}
\max_{\qty{\mathbf{X}_g}_{g=1}^G}\ \ \ \ &\sum_{g=1}^G\omega_gR_g\\
\text{s.t.}\ \ \ \ &\mathbf{X}_g\succeq0,\rank\pqty{\mathbf{X}_g}=1,\forall g\in\mathcal{G},\\
&\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{Ha}_g\leq P_T\\
&R_g\leq\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}\\
&\ \ \ \ \ \ \ \ \ \ -\frac{\sum_{g'\neq g}^G\mathbf{Q}_k:\pqty{\mathbf{X}_j-\mathbf{X}_j^\pqty{n}}}{\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}^\pqty{n}_{g'}}+\sigma_k^2}\ln2}
\end{align}
$$
