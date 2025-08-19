### ·System Model

1. Surface area is $\mathcal{A}$ and the transmitted signal at $\mathbf{r}\in\mathcal{A}$ can be expressed as $x(\mathbf{r})=\sum_{k=1}^K\mathrm{V}_k(\mathbf{r})s_k$, where $\mathrm{V}_k$ is distribution on BS surface to convey $s_k$, satisfying $\mathbb{E}\{|s_k|^2\}=1$.
2. $E_k=\int_\mathcal{A}\mathrm{H}_k(\mathbf{r})\sum_{i=1}^K\mathrm{V}_i(\mathbf{r})s_id\mathbf{r}+n_k$ denote the practical EF $k$-th user received and eventually derived the SINR as:

$$
\begin{flalign}
&\gamma_k\approx\frac{|\mathcal{A}_k|\cdot|\int_\mathcal{A}\mathrm{H}_k(\mathbf{r})\mathrm{V}_k(\mathbf{r})d\mathbf{r}|^2}{\sum_{i=1,i\neq k}^K|\mathcal{A}_i|\cdot|\int_\mathcal{A}\mathrm{H}_k(\mathbf{r})\mathrm{V}_i(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}\\
&\mathrm{H}_k(\mathbf{r})=\sqrt{\frac{\mathbf{e}_r^T(\mathbf{s}_k-\mathbf{r})}{\norm{\mathbf{r-s}_k}}}\cdot\frac{jk_0\eta e^{-jk_0\norm{\mathbf{r-s}_k}}}{4\pi\norm{\mathbf{r-s}_k}}(1+\frac{j/k_0}{\norm{\mathbf{r-s}_k}}-\frac{1/k_0^2}{\norm{\mathbf{r-s}_k}^2})
\end{flalign}
$$

​	Where $\mathbf{e}_r$ is the normal vector of CAPA surface at BS, $\eta=120\pi$ is impedance of free space and $k_0$ is wavenumber.

3. Objective function and power limit as following:

$$
\max_{\{\mathrm{V}_k(\mathbf{r})\}_{k=1}^K}\sum_{k=1}^Klog_2(1+\gamma_k)\\
\mathrm{s.t.}\ \ \ \sum_{k=1}^K\int_\mathcal{A}|\mathrm{V}_k(\mathbf{r})|^2d\mathbf{r}=P_{max}
$$

4. We can prove that optimal distribution can be expressed as: $\mathrm{V}^*_k(\mathbf{r})=\sum_{j=1}^Ka_{jk}\mathrm{H}_j(\mathbf{r})$, which can also be expressed by matrices $\mathbf{V=HA}$, where $\mathbf{V}=[\mathrm{V}_1,...,\mathrm{V}_K]\in\mathbb{C}^{1\times K}$, $\mathbf{H}=[\mathrm{H}_1,...,\mathrm{H}_K]\in\mathbb{C}^{1\times K}$ and $\mathbf{A}=[\mathbf{a}_1,...,\mathbf{a}_K]\in\mathbb{C}^{K\times K}$, $\mathbf{a}_k=[a_{1k},...,a_{Kk}]^T\in\mathbb{C}^{K\times 1}$.



#### System Model and Problem Formulation

​	We study a CAPA-based downlink multigroup multicast communication system, where a CAPA transmitter is equipped at both the BS and user to serve $G$ multicasting groups. We assume that there are a total of $K(K\ge G)$ users in all multicast groups. The users in the same group share the same information data, and the data of different groups is different and independent. We define the set of all multicast groups by $\mathcal{G}=\{1,2,...,G\}$, the set of all users by $\mathcal{K}=\{1,2,...,K\}$ and the user set belonging to group  $g\in\mathcal{G}$ by $\mathcal{K}_g$. Each user belongs to only one group, and each group contains at least one user.

###### ·Signal Model

​	We consider a CAPA transmitter with a continuous surface ST with an area of AT = |ST|, which contains sinusoidal source currents to emit EM waves for wireless communications.

​	We denote the aperture spaces of CAPAs at the BS and the $k$-th user as $\mathcal{S}$ and $\mathcal{S}_k\subseteq\mathbb{R}^{3\times 1}$, respectively, and the CAPA of $k$-th user is centered at $\mathbf{s}_k\in\mathbb{R}^{3\times 1}$. The transmit signal at point $\mathbf{r}\in\mathcal{S}$ can be expressed as
$$
x(\mathbf{r})=\sum_{g=1}^G\mathrm{J}_g(\mathbf{r})s_g,
$$
where $s_g$ is the symbol to be transmitted to the users of $g$-th group with $\mathbb{E}\{|s_g|^2\}=1$, and $\mathrm{J}_g(\mathbf{r})$ is the continuous associated current distribution to convey $s_g$. We denote $\mathbf{H}_k(\mathbf{r,s})\in\mathbb{C}$ as the channel response of point $\mathbf{s}\in\mathcal{S}_k$ from point $\mathbf{r}$ of BS's aperture and $E_k(\mathbf{s})\in\mathbb{C}$ as the excited electric field at point $\mathbf{s}\in\mathcal{S}_k$. Given both the propagation distance and BS aperture size, the aperture size of each user is generally negligible, and it is acceptable for each user to represent the channel response by that at point $\mathbf{s}_k$ from the point $\mathbf{r}$ of BS's aperture. Consequently,  we have $\mathrm{H}_k(\mathbf{r,s})\approx\mathrm{H}_k(\mathbf{r,s}_k)\triangleq\mathrm{H}_k(\mathbf{r})$, $E_k(\mathbf{s})\approx E_k(\mathbf{s}_k)\triangleq E_k$, the excited electric field at $\mathbf{s}\in\mathcal{S}_k$ within user $k$'s aperture can be approximately written as
$$
E_k\approx\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})x(\mathbf{r})d\mathbf{r}=\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\left(\sum_{g=1}^G\mathrm{J}_g(\mathbf{r})s_g\right)d\mathbf{r}
$$

If there are only line-of-sight (LoS) channels between the BS and the users, $\mathrm{H}_k(\mathbf{r})$ can be modeled as
$$
\mathrm{H}_k(\mathbf{r})=\sqrt{\frac{\mathbf{e}_r^T(\mathbf{s}_k-\mathbf{r})}{\norm{\mathbf{r-s}_k}}}\cdot\frac{jk_0\eta e^{-jk_0\norm{\mathbf{r-s}_k}}}{4\pi\norm{\mathbf{r-s}_k}}(1+\frac{j/k_0}{\norm{\mathbf{r-s}_k}}-\frac{1/k_0^2}{\norm{\mathbf{r-s}_k}^2}),
$$
where $\mathbf{e}_r\in\mathbb{R}^{3\times 1}$ is the normal vector of the BS's CAPA surface, $\eta=120\pi$ is the impedance of free space, $k_0=2\pi/\lambda$ with wavelength $\lambda$ denotes the wavenumber, and $(\cdot)^T$ denotes matrix transpose.

The total observed electric field $\mathrm{Y}_k(\mathbf{s})$ at point $\mathbf{s}\in\mathcal{S}_k$ is the sum of the information-carrying electric fields $E_k$ as well as a random noise field $\mathrm{N}_k(\mathbf{s})$, i.e.,
$$
\begin{flalign}
\mathrm{Y}_k
&=E_k+\mathrm{N}_k\\
&=\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\mathrm{J}_g(\mathbf{r})s_gd\mathbf{r}+\int_\mathcal{S}\mathrm{H}_k\pqty{\mathbf{r}}\sum_{g'\neq g}^G\mathrm{J}_{g'}(\mathbf{r})s_{g'}d\mathbf{r}+\mathrm{N}_k(\mathbf{s})
\end{flalign}
$$
where $\mathrm{N}_k(\mathbf{s})$ denotes the thermal noise and the noise field is modeled as a zero-mean complex Gaussian random process satisfying $\mathbb{E}\{\mathrm{N}(\mathbf{s})\mathrm{N}^*(\mathbf{s'})\}=\sigma^2\delta(\mathbf{s-s'})$, where $\delta(\cdot)$ represents Dirac delta function and $\sigma^2$ describes the noise intensity.

Therefore, the signal-to-noise-plus-interference ratio for decoding the desired signal at $k$-th user is given by
$$
\gamma_k=\frac{|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'\neq g}^K|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}_g(\mathbf{r})d\mathbf{r}|^2+\sigma_k^2}
$$
and the achievable rate of $k$-th user is thus given by $R_k=log_2(1+\gamma_k)$, according to Shannon's theorem.

However, due to the nature of the multicast, the practical achievable rate of $g$-th group is limited by the minimum user rate or worst channel condition in this group and defined as
$$
R_g=\min_{k\in\mathcal{K}_g}\{R_k\}
$$
###### ·Problem Formulation

In this paper, we aim to optimize the current distribution $\{\mathrm{J}_g(\mathbf{s})\}_{g=1}^G$ to maximize the sum-rate of the multicast system under the transmit power constrain. Consequently, the resulting optimization problem is reorganized as
$$
\begin{flalign}
\max_{\{\mathrm{J}_g(\mathbf{r})\}_{g=1}^G}\ \ \ &\sum_{g=1}^G\alpha_g\min_{k\in\mathcal{K}_g}\{\log_2(1+\gamma_k)\}\\
\mathrm{s.t.}\ \ \ &\sum_{g=1}^G\int_{\mathcal{S}_T}|\mathrm{J_g(\mathbf{r})|^2}d\mathbf{r}=P_T
\end{flalign}
$$

Where $\omega_g$ is the weight specified for $g$-th multicast group, given the requirements of fairness and quality of service (QoS). The total transmit power of the BS equipped with CAPA is limited by the constraint. Apparently, this problem is a non-convex problem, which is tricky to solve due to the following two reasons. First, in our problem, the optimization variable $\mathrm{J}_g(\mathbf{r})$ is a continuous function thanks to the new property of CAPA's ability to optimize the current distribution flexibly, thus, the conventional multicast beamforming schemes applied in MIMO cannot directly adapt to the optimization of the infinite dimensional current distribution in CAPA. Second, it is the continuous current distribution that brings us the objective function and constraints containing integrals, significantly improving the computational complexity. To tackle with above challenges, we proposed XXX method to optimize the current distribution of BS equipped with CAPA in multicast with a low computational cost.

### ·Approach

###### ·Problem Reformulation

In this section, we propose a SDR-based current distribution optimization scheme to maximize the sum-rate of multicast communication.

To optimize the objective function by adjusting the multigroup power allocation rather than increasing the transmit power, we introduce unconstrained current distribution $\mathrm{J}'_g(\mathbf{r})$ to define a new form objective function as
$$
\gamma'_k=\frac{|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}'_g(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'\neq g}^K|\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{J}'_g(\mathbf{r})d\mathbf{r}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\int_\mathcal{S}|\mathrm{J}'_g(\mathbf{r})|^2d\mathbf{r}}
$$
the corresponding objective function is
$$
\max_{\{\mathrm{J}'_g(\mathbf{r})\}_{g=1}^G}\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\{\log_2(1+\gamma'_k)\}
$$
If denote $\bar{\mathrm{J}}'_g(\mathbf{r})$ as the optimal solution to the problem (), then the optimal solution $\bar{\mathrm{J}}_g(\mathbf{r})$ to the original problem () can be derived by $\bar{\mathrm{J}}_g(\mathbf{r})=C\bar{\mathrm{J}}'_g(\mathbf{r})$, where $C=\sqrt{\frac{P_T}{\sum_{g=1}^G\int_\mathcal{S}|\bar{\mathrm{J}}'_g(\mathbf{r})|^2d\mathbf{r}}}$.

Considering that the optimal current distribution is in a function subspace spanned by $\{\mathrm{H}_k(\mathbf{r})\}_{k=1}^K$, the $\mathrm{J}'_g(\mathbf{r})$ can thus be expressed as $\mathrm{J}'_g(\mathbf{r})=\sum_{i=1}^Ka_{i,g}\mathrm{H}_i(\mathbf{r})$. Consequently, the SINR can be rewritten as
$$
\begin{align}
\gamma'_k&=\frac{|\sum_{i=1}^Ka_{i,g}\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{H}^*_i(\mathbf{r})d\mathbf{r}|^2}{\sum_{g'\neq g}^K|\sum_{i=1}^Ka_{i,g'}\int_\mathcal{S}\mathrm{H}_k(\mathbf{r})\mathrm{H}^*_i(\mathbf{r})d\mathbf{r}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\int_\mathcal{S}\left(\sum_{i=1}^Ka_{i,g}\mathrm{H}_i(\mathbf{r})\right)\left(\sum_{i=1}^Ka^*_{i,g}\mathrm{H}^*_i(\mathbf{r})\right)d\mathbf{r}}\\
&=\frac{|\sum_{i=1}^Ka_{ig}h_{k,i}|^2}{\sum_{g'=1,g'\neq g}^G|\sum_{i=1}^Ka_{ig'}h_{k,i}|^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\sum_{i=1}^K\sum_{j=1}^Ka_{i,g}a^*_{j,g}h_{i,j}}\\
&=\frac{\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2}{\sum_{g'\neq g}^G\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_{g'}}^2+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\mathbf{a}_g^\mathrm{H}\mathbf{H}\mathbf{a}_g}
\end{align}
$$
where $\mathbf{H}=[\mathbf{h}_1,\cdots,\mathbf{h}_K]$, $\mathbf{h}_k=[h_{k,1},\cdots,h_{k,K}]^T$, $\mathbf{a}_g=[a_{g,1},\cdots,a_{g,K}]^T$ and $h_{i,j/j,i}\triangleq\int_\mathcal{S}\mathrm{H}_i(\mathbf{r})\mathrm{H}_j^*(\mathbf{r})d\mathbf{r}$. Obviously, since the current distribution $\{\mathrm{J}_g(\mathbf{r})\}_{g=1}^G$ can be represented by the matrix $\mathbf{A}=[\mathbf{a}_1,\cdots,\mathbf{a}_G]$, problem () can be converted to
$$
\max_{\mathbf{A}}\sum_{g=1}^G\omega_g\min_{k\in\mathcal{K}_g}\{\log_2(1+\gamma'_k)\}
$$
​	Assuming that the optimal solution to problem () is $\bar{\mathbf{A}}=[\bar{\mathbf{a}}_1,\cdots,\bar{\mathbf{a}}_G]$, where $\bar{\mathbf{a}}_g=[\bar{a}_{g,1},\cdots,\bar{a}_{g,K}]^T$, then that of problem () is $\bar{\mathrm{J}}'_g(\mathbf{r})=\sum_{i=1}^G\bar{a}_{g,i}\mathrm{H}_i(\mathbf{r})$, and that of original problem () is thus derived as
$$
\bar{\mathrm{J}}_g(\mathbf{r})=\sqrt{\frac{P_T}{\sum_{g=1}^G\int_\mathcal{S}|\bar{\mathrm{J}}'_g(\mathbf{r})|^2d\mathbf{r}}}\bar{\mathrm{J}}'_g(\mathbf{r})
$$
Based on the above analysis, we transform the infinite dimensional current distribution optimization problem into a finite dimensional optimization problem of matrix $\mathbf{A}$.

###### ·SDR-based Scheme

To solve the above matrix $\mathbf{A}$ optimization problem by adopting SDR, we convert the above objective function into that composed of matrix trace by the property of trace. Concretely, we express the received signal as $\vqty{\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}^2=\mathbf{a}_g^\mathrm{H}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g=\tr\pqty{\mathbf{a}_g^\mathrm{H}\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g}=\tr\pqty{\mathbf{h}_k\mathbf{h}_k^\mathrm{H}\mathbf{a}_g\mathbf{a}_g^\mathrm{H}}$, given the property of  $\tr\pqty{\mathbf{AB}}=\tr\pqty{\mathbf{BA}}$. Subsequently, we define $\mathbf{Q}_k=\mathbf{h}_k\mathbf{h}_k^\mathrm{H}$ and $\mathbf{X}_g=\mathbf{a}_g\mathbf{a}_g^\mathrm{H}$, and note that $\mathbf{X}_g=\mathbf{a}_g\mathbf{a}_g^\mathrm{H}$ holds if and only if $\mathbf{X}_g\succeq0$ and $\text{Rank}\pqty{\mathbf{X}_g}=1$, then the received signal is given by $\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}$, the transmit power is similarly given by $\tr\pqty{\mathbf{HX}_g}$, thus the objective function and the transmit power constraint can be derived as
$$
\begin{align}
\max_{\qty{R_g}_{g=1}^G} &\sum_{g=1}^G\omega_gR_g\label{trace_sum_rate},\\
\mathrm{s.t.} &\sum_{g=1}^G\tr\pqty{\mathbf{H}\mathbf{X}_g}=P_T\label{trace_power},
\end{align}
$$
​	Where $R_g=\min_{k\in\mathcal{K}_g}\qty{R_k},\forall g\in\mathcal{G}$ is the achievable rate of $g$-th multicast group, satisfying
$$
\begin{align}
R_g\leq&\log_2\pqty{\frac{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}}\\
\leq&\log_2\pqty{\sum_{g=1}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\sigma_k^2}-\log_2\pqty{\sum_{g'\neq g}^G\tr\pqty{\mathbf{Q}_k\mathbf{X}_{g'}}+\sigma_k^2}
\end{align}
$$
​	However, to further address the non-convexity of $R_g$, the Taylor Expansion is adapted. Specifically, the first-order Taylor expansion of $\varepsilon_k$ is introduced to transform () into
$$
R_g\leq
$$






 $\tau_k$ and $\varepsilon_k$ are auxiliary variables defined as follow
$$
\begin{align}
e^{\tau_k}&=\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}\\
e^{\varepsilon_k}&=\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}\\
\end{align}
$$



​	Obviously, () is a convex function with the constraints of $\tau_k$ and $\varepsilon_k$ given as
$$
\begin{align}
e^{\tau_k}&\leq\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g},\\
e^{\varepsilon_k}&\geq\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g},
\end{align}
$$
​	where inequalities hold equality at optimum point. With respect to the non-convexity of (), we further substitute the first-order Taylor expansion $e^{\varepsilon_k^n}\pqty{\varepsilon_k-\varepsilon_k^n+1}$ for $e^{\varepsilon_k}$ to ensure the convexity of constraints. Eventually, the optimization problem can be reformulated as
$$
\begin{align}
&\max_{\qty{\tau_k,\varepsilon_k}_{k=1}^K,\qty{\mathbf{X}_g}_{g=1}^G}\ \sum_{g=1}^G\min_{k\in\mathcal{K}_g}\qty{\tau_k-\varepsilon_k}\\
&\text{s.t.}\ \mathbf{X}_g\succeq0,\ \text{Rank}\pqty{\mathbf{X}_g}=1, \forall g\in\mathcal{G},\\
&e^{\tau_k}\leq\sum_{g=1}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g},\\
&e^{\varepsilon_k^n}\pqty{\varepsilon_k-\varepsilon_k^n+1}\geq\sum_{g'\neq g}^G\Tr\pqty{\mathbf{Q}_k\mathbf{X}_g}+\frac{\sigma_k^2}{P_T}\sum_{g=1}^G\Tr\pqty{\mathbf{H}\mathbf{X}_g}.
\end{align}
$$
​	To obtain a complete convex problem, we temporarily relax the non-convex constraint of $\text{Rank}\pqty{\mathbf{X}_g}=1,\forall g\in\mathcal{G}$ in () aided with SDR. Thereby, the problem () becomes a convex optimization problem being solved by CVX solver. Furthermore, since the optimal solution $\qty{\bar{\mathbf{X}}_g}_{g=1}^G$ is not always content with the $\text{rank-1}$, the approximation approaches can be exploited to calculate a feasible solution $\qty{\bar{\mathbf{a}}_g}_{g=1}^G$ as the optimal solution to problem (), and that to problem () will be thus derived from ().
