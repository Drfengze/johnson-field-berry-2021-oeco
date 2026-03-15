# Python Variable Name Lookup Table v2

## 1. Naming Rule

- `sif.ipynb` 旧表里已经有的标准名，直接沿用。
- 新增变量尽量贴近旧风格，不用超长全拼。
- 这一版以“最终标准名”作为后续 Python 代码规范化的候选版本。

## 2. Core Inputs And Shared Variables

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `pathway_opt` | `pathway_option` | 光合途径选项，当前支持 `C3`、`Type-I-C3-C4`、`NADP-ME-C4` |
| `Q` | `PPFD` | 光合有效辐射通量密度 |
| `T` | `Temp` | 叶温 |
| `P` | `Pressure` | 总压 |
| `C_m` | `CO2_m` | 叶肉细胞 CO2 分压 |
| `O_m` | `O2_m` | 叶肉细胞 O2 分压 |
| `Abs` | `Absorptance` | 叶片总吸收率 |
| `abs_frac` | `Abs_fraction_s` | 维管束鞘吸收占总吸收的比例 |
| `a2_m_frac` | `PSII_fraction` | 叶肉 PSII 吸收占比 |
| `a2_s_frac` | `PSII_fraction_s` | 维管束鞘 PSII 吸收占比 |
| `eps1` | `epsilon_PSI` | PSI 到探测器的传递函数 |
| `eps2` | `epsilon_PSII` | PSII 到探测器的传递函数 |
| `CB6F` | `Cytbf_density` | Cyt b6f 密度 |
| `vq_frac` | `Cytbf_fraction_s` | 维管束鞘 Cyt b6f 占比 |
| `RUB` | `Rubisco_density` | Rubisco 位点密度 |
| `Rdsc` | `Resp_scalar` | 线粒体暗呼吸标量 |
| `vc_frac` | `Rubisco_fraction_s` | 维管束鞘 Rubisco 占比 |
| `Vpmax` | `V_p_max` | PEP carboxylase 最大活性 |
| `gbs` | `g_bs_CO2` | 维管束鞘对 CO2 的传导度 |
| `gbso` | `g_bs_O2` | 维管束鞘对 O2 的传导度 |
| `Kf` | `k_F` | 荧光速率常数 |
| `Kd` | `k_D` | 基础热耗散速率常数 |
| `Kp1` | `k_P1` | PSI 光化学速率常数 |
| `Kn1` | `k_N1` | PSI 调节性热耗散速率常数 |
| `Kp2` | `k_P2` | PSII 光化学速率常数 |
| `Ku2` | `k_U2` | PSII 激子共享速率常数 |
| `kq` | `k_q` | Cyt b6f 对 PQH2 的催化常数 |
| `nl` | `ATP_e_ratio_linear` | 线性电子流每电子 ATP 产额 |
| `nc` | `ATP_e_ratio_cyclic` | 环式电子流每电子 ATP 产额 |
| `kc` | `k_cat_CO2` | Rubisco 对 CO2 的 kcat |
| `ko` | `k_cat_O2` | Rubisco 对 O2 的 kcat |
| `Kc` | `K_m_CO2` | Rubisco 对 CO2 的 Km |
| `Ko` | `K_m_O2` | Rubisco 对 O2 的 Km |
| `Kp` | `K_m_PEPC_CO2` | PEPC 对 CO2 的 Km |
| `Abs_m` | `Absorptance_m` | 叶肉总吸收率 |
| `Abs_s` | `Absorptance_s` | 维管束鞘总吸收率 |
| `a1_m` | `a_PSI_m` | 叶肉 PSI 吸收截面 |
| `a2_m` | `a_PSII_m` | 叶肉 PSII 吸收截面 |
| `a1_s` | `a_PSI_s` | 维管束鞘 PSI 吸收截面 |
| `a2_s` | `a_PSII_s` | 维管束鞘 PSII 吸收截面 |
| `T_K` | `Temp_K` | 开尔文温度 |
| `Tref_K` | `Temp_ref_K` | 参考温度（K） |
| `Vqmax` | `V_q_max` | 全叶 Cyt b6f 最大活性 |
| `Vqmax_m` | `V_q_max_m` | 叶肉 Cyt b6f 最大活性 |
| `Vqmax_s` | `V_q_max_s` | 维管束鞘 Cyt b6f 最大活性 |
| `Rd` | `R_d` | 全叶暗呼吸 |
| `Rd_m` | `R_d_m` | 叶肉暗呼吸 |
| `Rd_s` | `R_d_s` | 维管束鞘暗呼吸 |
| `S` | `Specificity` | Rubisco 对 CO2/O2 的特异性 |
| `Vcmax` | `V_c_max` | 全叶 Rubisco 最大羧化能力 |
| `Vcmax_m` | `V_c_max_m` | 叶肉 Rubisco 最大羧化能力 |
| `Vcmax_s` | `V_c_max_s` | 维管束鞘 Rubisco 最大羧化能力 |
| `Vpmax_m` | `V_p_max_m` | 叶肉 PEPC 最大活性 |
| `CB6F_m` | `Cytbf_density_m` | 叶肉 Cyt b6f 密度 |
| `CB6F_s` | `Cytbf_density_s` | 维管束鞘 Cyt b6f 密度 |
| `RUB_m` | `Rubisco_density_m` | 叶肉 Rubisco 密度 |
| `RUB_s` | `Rubisco_density_s` | 维管束鞘 Rubisco 密度 |
| `CB6F_ma` | `Cytbf_active_m_actual` | 叶肉实际活化的 Cyt b6f |
| `CB6F_sa` | `Cytbf_active_s_actual` | 维管束鞘实际活化的 Cyt b6f |

## 3. Symbolic Solver And State Selector Variables

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `c3c4_solve_jj` | `solve_C3C4_jj` | Type-I-C3-C4 情况下，`jj` 状态的符号求解器 |
| `c3c4_solve_jc` | `solve_C3C4_jc` | Type-I-C3-C4 情况下，`jc` 状态的符号求解器 |
| `c3c4_solve_cj` | `solve_C3C4_cj` | Type-I-C3-C4 情况下，`cj` 状态的符号求解器 |
| `c3c4_solve_cc` | `solve_C3C4_cc` | Type-I-C3-C4 情况下，`cc` 状态的符号求解器 |
| `c4_solve_jj` | `solve_C4_jj` | NADP-ME-C4 情况下，`jj` 状态的符号求解器 |
| `c4_solve_jc` | `solve_C4_jc` | NADP-ME-C4 情况下，`jc` 状态的符号求解器 |
| `c4_solve_cj` | `solve_C4_cj` | NADP-ME-C4 情况下，`cj` 状态的符号求解器 |
| `c4_solve_cc` | `solve_C4_cc` | NADP-ME-C4 情况下，`cc` 状态的符号求解器 |
| `which_JP700_ma` | `which_J_PSI_m` | 叶肉实际态选择器，区分 `m_j` 和 `m_c` |
| `which_JP700_sj` | `which_J_PSI_s_j` | 当叶肉为 `j` 时，维管束鞘在 `s_jj` 和 `s_jc` 间的选择器 |
| `which_JP700_sc` | `which_J_PSI_s_c` | 当叶肉为 `c` 时，维管束鞘在 `s_cj` 和 `s_cc` 间的选择器 |

## 4. Mesophyll State Variables

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `JP700_mj` | `J_PSI_m_j` | 叶肉 `j` 限制下的 PSI 电子传递速率 |
| `JP680_mj` | `J_PSII_m_j` | 叶肉 `j` 限制下的 PSII 电子传递速率 |
| `Vc_mj` | `V_c_m_j` | 叶肉 `j` 限制下的 Rubisco 羧化速率 |
| `Vo_mj` | `V_o_m_j` | 叶肉 `j` 限制下的 Rubisco 氧化速率 |
| `Ag_mj` | `A_gross_m_j` | 叶肉 `j` 限制下的总同化速率 |
| `An_mj` | `A_net_m_j` | 叶肉 `j` 限制下的净同化速率 |
| `Vp_mj` | `V_p_m_j` | 叶肉 `j` 限制下的 PEPC / C4 泵速率 |
| `Vc_mc` | `V_c_m_c` | 叶肉 `c` 限制下的 Rubisco 羧化速率 |
| `Vo_mc` | `V_o_m_c` | 叶肉 `c` 限制下的 Rubisco 氧化速率 |
| `Ag_mc` | `A_gross_m_c` | 叶肉 `c` 限制下的总同化速率 |
| `An_mc` | `A_net_m_c` | 叶肉 `c` 限制下的净同化速率 |
| `JP680_mc` | `J_PSII_m_c` | 叶肉 `c` 限制下的 PSII 电子传递速率 |
| `JP700_mc` | `J_PSI_m_c` | 叶肉 `c` 限制下的 PSI 电子传递速率 |
| `Vp_mc` | `V_p_m_c` | 叶肉 `c` 限制下的 PEPC / C4 泵速率 |
| `Vg_mj` | `V_g_m_j` | 叶肉 `j` 限制下的 glycine shuttling 速率 |
| `Vg_mc` | `V_g_m_c` | 叶肉 `c` 限制下的 glycine shuttling 速率 |
| `JP700_ma` | `J_PSI_m_actual` | 叶肉实际 PSI 电子传递速率 |
| `JP680_ma` | `J_PSII_m_actual` | 叶肉实际 PSII 电子传递速率 |
| `Vg_ma` | `V_g_m_actual` | 叶肉实际 glycine shuttling 速率 |
| `Vp_ma` | `V_p_m_actual` | 叶肉实际 PEPC / C4 泵速率 |
| `An_ma` | `A_net_m_actual` | 叶肉实际净同化速率 |
| `Ag_ma` | `A_gross_m_actual` | 叶肉实际总同化速率 |

## 5. Bundle-Sheath State Variables

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `JP700_sjj` | `J_PSI_s_jj` | 维管束鞘在 `jj` 状态下的 PSI 电子传递速率 |
| `JP680_sjj` | `J_PSII_s_jj` | 维管束鞘在 `jj` 状态下的 PSII 电子传递速率 |
| `An_sjj` | `A_net_s_jj` | 维管束鞘在 `jj` 状态下的净同化速率 |
| `C_sjj` | `CO2_s_jj` | 维管束鞘在 `jj` 状态下的 CO2 分压 |
| `O_sjj` | `O2_s_jj` | 维管束鞘在 `jj` 状态下的 O2 分压 |
| `Ag_sjj` | `A_gross_s_jj` | 维管束鞘在 `jj` 状态下的总同化速率 |
| `Vc_sjj` | `V_c_s_jj` | 维管束鞘在 `jj` 状态下的 Rubisco 羧化速率 |
| `Vo_sjj` | `V_o_s_jj` | 维管束鞘在 `jj` 状态下的 Rubisco 氧化速率 |
| `JP700_sjc` | `J_PSI_s_jc` | 维管束鞘在 `jc` 状态下的 PSI 电子传递速率 |
| `JP680_sjc` | `J_PSII_s_jc` | 维管束鞘在 `jc` 状态下的 PSII 电子传递速率 |
| `An_sjc` | `A_net_s_jc` | 维管束鞘在 `jc` 状态下的净同化速率 |
| `C_sjc` | `CO2_s_jc` | 维管束鞘在 `jc` 状态下的 CO2 分压 |
| `O_sjc` | `O2_s_jc` | 维管束鞘在 `jc` 状态下的 O2 分压 |
| `Ag_sjc` | `A_gross_s_jc` | 维管束鞘在 `jc` 状态下的总同化速率 |
| `Vc_sjc` | `V_c_s_jc` | 维管束鞘在 `jc` 状态下的 Rubisco 羧化速率 |
| `Vo_sjc` | `V_o_s_jc` | 维管束鞘在 `jc` 状态下的 Rubisco 氧化速率 |
| `JP700_scj` | `J_PSI_s_cj` | 维管束鞘在 `cj` 状态下的 PSI 电子传递速率 |
| `JP680_scj` | `J_PSII_s_cj` | 维管束鞘在 `cj` 状态下的 PSII 电子传递速率 |
| `An_scj` | `A_net_s_cj` | 维管束鞘在 `cj` 状态下的净同化速率 |
| `C_scj` | `CO2_s_cj` | 维管束鞘在 `cj` 状态下的 CO2 分压 |
| `O_scj` | `O2_s_cj` | 维管束鞘在 `cj` 状态下的 O2 分压 |
| `Ag_scj` | `A_gross_s_cj` | 维管束鞘在 `cj` 状态下的总同化速率 |
| `Vc_scj` | `V_c_s_cj` | 维管束鞘在 `cj` 状态下的 Rubisco 羧化速率 |
| `Vo_scj` | `V_o_s_cj` | 维管束鞘在 `cj` 状态下的 Rubisco 氧化速率 |
| `JP700_scc` | `J_PSI_s_cc` | 维管束鞘在 `cc` 状态下的 PSI 电子传递速率 |
| `JP680_scc` | `J_PSII_s_cc` | 维管束鞘在 `cc` 状态下的 PSII 电子传递速率 |
| `An_scc` | `A_net_s_cc` | 维管束鞘在 `cc` 状态下的净同化速率 |
| `C_scc` | `CO2_s_cc` | 维管束鞘在 `cc` 状态下的 CO2 分压 |
| `O_scc` | `O2_s_cc` | 维管束鞘在 `cc` 状态下的 O2 分压 |
| `Ag_scc` | `A_gross_s_cc` | 维管束鞘在 `cc` 状态下的总同化速率 |
| `Vc_scc` | `V_c_s_cc` | 维管束鞘在 `cc` 状态下的 Rubisco 羧化速率 |
| `Vo_scc` | `V_o_s_cc` | 维管束鞘在 `cc` 状态下的 Rubisco 氧化速率 |
| `JP700_sj` | `J_PSI_s_j` | 当叶肉为 `j` 时选中的维管束鞘 PSI 电子传递速率 |
| `JP680_sj` | `J_PSII_s_j` | 当叶肉为 `j` 时选中的维管束鞘 PSII 电子传递速率 |
| `C_sj` | `CO2_s_j` | 当叶肉为 `j` 时选中的维管束鞘 CO2 分压 |
| `O_sj` | `O2_s_j` | 当叶肉为 `j` 时选中的维管束鞘 O2 分压 |
| `An_sj` | `A_net_s_j` | 当叶肉为 `j` 时选中的维管束鞘净同化速率 |
| `JP700_sc` | `J_PSI_s_c` | 当叶肉为 `c` 时选中的维管束鞘 PSI 电子传递速率 |
| `JP680_sc` | `J_PSII_s_c` | 当叶肉为 `c` 时选中的维管束鞘 PSII 电子传递速率 |
| `C_sc` | `CO2_s_c` | 当叶肉为 `c` 时选中的维管束鞘 CO2 分压 |
| `O_sc` | `O2_s_c` | 当叶肉为 `c` 时选中的维管束鞘 O2 分压 |
| `An_sc` | `A_net_s_c` | 当叶肉为 `c` 时选中的维管束鞘净同化速率 |
| `JP700_sa` | `J_PSI_s_actual` | 维管束鞘实际 PSI 电子传递速率 |
| `JP680_sa` | `J_PSII_s_actual` | 维管束鞘实际 PSII 电子传递速率 |
| `Ag_sa` | `A_gross_s_actual` | 维管束鞘实际总同化速率 |
| `An_sa` | `A_net_s_actual` | 维管束鞘实际净同化速率 |
| `C_sa` | `CO2_s_actual` | 维管束鞘实际 CO2 分压 |
| `O_sa` | `O2_s_actual` | 维管束鞘实际 O2 分压 |
| `L_C_sa` | `Leak_CO2_s_actual` | 维管束鞘实际 CO2 泄漏速率 |

## 6. Whole-Leaf Electron Transport And Assimilation

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `JP700_a` | `J_PSI_actual` | 全叶实际 PSI 电子传递速率 |
| `JP680_a` | `J_PSII_actual` | 全叶实际 PSII 电子传递速率 |
| `An_a` | `A_net_actual` | 全叶实际净同化速率 |
| `Ag_a` | `A_gross_actual` | 全叶实际总同化速率 |

## 7. Fluorescence Yield And Quenching Variables

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `phi1P_ma` | `phi_P1_m_actual` | 叶肉实际 PSI 光化学量子产额 |
| `phi1N_ma` | `phi_N1_m_actual` | 叶肉实际 PSI 调节性热耗散量子产额 |
| `phi1D_ma` | `phi_D1_m_actual` | 叶肉实际 PSI 基础热耗散量子产额 |
| `phi1F_ma` | `phi_F1_m_actual` | 叶肉实际 PSI 荧光量子产额 |
| `phi1P_sa` | `phi_P1_s_actual` | 维管束鞘实际 PSI 光化学量子产额 |
| `phi1N_sa` | `phi_N1_s_actual` | 维管束鞘实际 PSI 调节性热耗散量子产额 |
| `phi1D_sa` | `phi_D1_s_actual` | 维管束鞘实际 PSI 基础热耗散量子产额 |
| `phi1F_sa` | `phi_F1_s_actual` | 维管束鞘实际 PSI 荧光量子产额 |
| `phi2p_ma` | `phi_p2_m_actual` | 叶肉实际 PSII 基础光化学量子产额 |
| `phi2n_ma` | `phi_n2_m_actual` | 叶肉实际 PSII 基础调节性热耗散量子产额 |
| `phi2d_ma` | `phi_d2_m_actual` | 叶肉实际 PSII 基础热耗散量子产额 |
| `phi2f_ma` | `phi_f2_m_actual` | 叶肉实际 PSII 基础荧光量子产额 |
| `phi2u_ma` | `phi_u2_m_actual` | 叶肉实际 PSII 基础激子共享量子产额 |
| `phi2P_ma` | `phi_P2_m_actual` | 叶肉实际 PSII 真实光化学量子产额 |
| `phi2N_ma` | `phi_N2_m_actual` | 叶肉实际 PSII 真实调节性热耗散量子产额 |
| `phi2D_ma` | `phi_D2_m_actual` | 叶肉实际 PSII 真实热耗散量子产额 |
| `phi2F_ma` | `phi_F2_m_actual` | 叶肉实际 PSII 真实荧光量子产额 |
| `phi2p_sa` | `phi_p2_s_actual` | 维管束鞘实际 PSII 基础光化学量子产额 |
| `phi2n_sa` | `phi_n2_s_actual` | 维管束鞘实际 PSII 基础调节性热耗散量子产额 |
| `phi2d_sa` | `phi_d2_s_actual` | 维管束鞘实际 PSII 基础热耗散量子产额 |
| `phi2f_sa` | `phi_f2_s_actual` | 维管束鞘实际 PSII 基础荧光量子产额 |
| `phi2u_sa` | `phi_u2_s_actual` | 维管束鞘实际 PSII 基础激子共享量子产额 |
| `phi2P_sa` | `phi_P2_s_actual` | 维管束鞘实际 PSII 真实光化学量子产额 |
| `phi2N_sa` | `phi_N2_s_actual` | 维管束鞘实际 PSII 真实调节性热耗散量子产额 |
| `phi2D_sa` | `phi_D2_s_actual` | 维管束鞘实际 PSII 真实热耗散量子产额 |
| `phi2F_sa` | `phi_F2_s_actual` | 维管束鞘实际 PSII 真实荧光量子产额 |
| `q1_ma` | `q_P1_m_actual` | 叶肉实际 PSI 光化学淬灭 |
| `q2_ma` | `q_P2_m_actual` | 叶肉实际 PSII 光化学淬灭 |
| `Kn2_ma` | `k_N2_m_actual` | 叶肉实际 PSII 调节性热耗散速率常数 |
| `q1_sa` | `q_P1_s_actual` | 维管束鞘实际 PSI 光化学淬灭 |
| `q2_sa` | `q_P2_s_actual` | 维管束鞘实际 PSII 光化学淬灭 |
| `Kn2_sa` | `k_N2_s_actual` | 维管束鞘实际 PSII 调节性热耗散速率常数 |

## 8. Fluorescence Level And PAM Variables

| Current Var | Final Standard Name | Description |
| --- | --- | --- |
| `Fs_ma` | `F_s_m_actual` | 叶肉实际稳态荧光 |
| `Fm_ma` | `F_m_m_actual` | 叶肉实际暗适应最大荧光 |
| `Fo_ma` | `F_o_m_actual` | 叶肉实际暗适应最小荧光 |
| `Fmp_ma` | `F_m_prime_m_actual` | 叶肉实际光适应最大荧光 |
| `Fop_ma` | `F_o_prime_m_actual` | 叶肉实际光适应最小荧光 |
| `Fs_sa` | `F_s_s_actual` | 维管束鞘实际稳态荧光 |
| `Fm_sa` | `F_m_s_actual` | 维管束鞘实际暗适应最大荧光 |
| `Fo_sa` | `F_o_s_actual` | 维管束鞘实际暗适应最小荧光 |
| `Fmp_sa` | `F_m_prime_s_actual` | 维管束鞘实际光适应最大荧光 |
| `Fop_sa` | `F_o_prime_s_actual` | 维管束鞘实际光适应最小荧光 |
| `Fs_a` | `F_s` | 全叶稳态荧光 |
| `Fm_a` | `F_m` | 全叶暗适应最大荧光 |
| `Fo_a` | `F_o` | 全叶暗适应最小荧光 |
| `Fmp_a` | `F_m_prime` | 全叶光适应最大荧光 |
| `Fop_a` | `F_o_prime` | 全叶光适应最小荧光 |
| `PAM1_a` | `qP` | 由当前代码公式得到的 qP |
| `PAM2_a` | `qL` | 由当前代码公式得到的 qL |
| `PAM3_a` | `NPQ` | 由当前代码公式得到的 NPQ |
| `PAM4_a` | `Phi_PSII` | 由当前代码公式得到的 PSII 有效量子产额 |
| `PAM5_a` | `Phi_NPQ` | 由当前代码公式得到的非光化学耗散量子产额 |
| `PAM6_a` | `Phi_f_d` | 由当前代码公式得到的荧光与基础耗散项 |

## 9. Notes

- 这版已经把新增变量显式展开了，不再依赖 `Suffix Mapping` 去猜。
- 旧表中已有的名字，比如 `PPFD`、`Temp`、`V_c_max`、`A_net_actual`、`F_s`，这版都保留。
- 新加的命名尽量向旧风格靠拢，主要补的是：
  - `Pressure`
  - `CO2_m` / `O2_m`
  - `CO2_s_*` / `O2_s_*`
  - `J_PSI_*` / `J_PSII_*`
  - `A_net_s_*` / `A_gross_s_*`
  - `solve_C3C4_*` / `solve_C4_*`
- `PAM1_a` 到 `PAM6_a` 的含义，这里按当前 [model_fun_c3c4.py](/Users/huaizefeng/Documents/GitHub/johnson-field-berry-2021-oeco/src/scripts/model_fun_c3c4.py) 和 [model_fun_c3c4.m](/Users/huaizefeng/Documents/GitHub/johnson-field-berry-2021-oeco/scripts/model_fun_c3c4.m) 的公式解释，不按旧 notebook 的编号解释。
