f_pi [GeV]: $\frac{93}{1000}$
eta generator: $\left[\begin{matrix}- \frac{\sqrt{3} \sin{\left(\theta_{P} \right)}}{3} + \frac{\sqrt{6} \cos{\left(\theta_{P} \right)}}{6} & 0 & 0\\0 & - \frac{\sqrt{3} \sin{\left(\theta_{P} \right)}}{3} + \frac{\sqrt{6} \cos{\left(\theta_{P} \right)}}{6} & 0\\0 & 0 & - \frac{\sqrt{3} \sin{\left(\theta_{P} \right)}}{3} - \frac{\sqrt{6} \cos{\left(\theta_{P} \right)}}{3}\end{matrix}\right]$
eta-prime generator: $\left[\begin{matrix}\frac{\sqrt{6} \sin{\left(\theta_{P} \right)}}{6} + \frac{\sqrt{3} \cos{\left(\theta_{P} \right)}}{3} & 0 & 0\\0 & \frac{\sqrt{6} \sin{\left(\theta_{P} \right)}}{6} + \frac{\sqrt{3} \cos{\left(\theta_{P} \right)}}{3} & 0\\0 & 0 & - \frac{\sqrt{6} \sin{\left(\theta_{P} \right)}}{3} + \frac{\sqrt{3} \cos{\left(\theta_{P} \right)}}{3}\end{matrix}\right]$
generic diagonal C_A: $\left[\begin{matrix}c_{u} & 0 & 0\\0 & c_{d} & 0\\0 & 0 & c_{s}\end{matrix}\right]$
Delta-S spurion C_A=c_sd E_sd: $\left[\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & c_{sd} & 0\end{matrix}\right]$
Vector convention: ideal omega=diag(1,1,0)/sqrt(2), phi=diag(0,0,1).
Pair convention: pair1 precedes pair2 in MESON_ORDER.

Enumerated channels: 33 (charge Q=0, total S=0 or +1).

[diag C_A]  required nonzero final S=0
channel | nonzero spectator -> Vbar-pair : trace product
-----------|------------
$\pi^+ \pi^- \pi^0$ | $a \to \pi^+ \rho^-[\to \pi^-\,\pi^0]$ : $\sqrt{2} \left(c_{d} - c_{u}\right)$; $a \to \pi^- \rho^+[\to \pi^+\,\pi^0]$ : $\sqrt{2} \left(c_{d} - c_{u}\right)$
$\pi^+ \pi^- \eta$ | --
$\pi^+ \pi^- \eta'$ | --
$\pi^+ K^- K^0$ | $a \to \pi^+ \rho^-[\to K^-\,K^0]$ : $- c_{d} + c_{u}$; $a \to K^- K^{*+}[\to \pi^+\,K^0]$ : $- c_{s} + c_{u}$; $a \to K^0 \bar{K}^{*0}[\to \pi^+\,K^-]$ : $c_{d} - c_{s}$
$\pi^- K^+ \bar{K}^0$ | $a \to \pi^- \rho^+[\to K^+\,\bar{K}^0]$ : $- c_{d} + c_{u}$; $a \to K^+ K^{*-}[\to \pi^-\,\bar{K}^0]$ : $- c_{s} + c_{u}$; $a \to \bar{K}^0 K^{*0}[\to \pi^-\,K^+]$ : $c_{d} - c_{s}$
$\pi^0 \pi^0 \pi^0$ | --
$\pi^0 \pi^0 \eta$ | --
$\pi^0 \pi^0 \eta'$ | --
$\pi^0 K^+ K^-$ | $a \to K^+ K^{*-}[\to \pi^0\,K^-]$ : $- \frac{\sqrt{2} \left(c_{s} - c_{u}\right)}{2}$; $a \to K^- K^{*+}[\to \pi^0\,K^+]$ : $- \frac{\sqrt{2} \left(c_{s} - c_{u}\right)}{2}$
$\pi^0 K^0 \bar{K}^0$ | $a \to K^0 \bar{K}^{*0}[\to \pi^0\,\bar{K}^0]$ : $- \frac{\sqrt{2} \left(c_{d} - c_{s}\right)}{2}$; $a \to \bar{K}^0 K^{*0}[\to \pi^0\,K^0]$ : $- \frac{\sqrt{2} \left(c_{d} - c_{s}\right)}{2}$
$\pi^0 \eta \eta$ | --
$\pi^0 \eta \eta'$ | --
$\pi^0 \eta' \eta'$ | --
$K^+ K^- \eta$ | $a \to K^+ K^{*-}[\to K^-\,\eta]$ : $\frac{\sqrt{6} \left(c_{s} - c_{u}\right) \cos{\left(\theta_{P} \right)}}{2}$; $a \to K^- K^{*+}[\to K^+\,\eta]$ : $\frac{\sqrt{6} \left(c_{s} - c_{u}\right) \cos{\left(\theta_{P} \right)}}{2}$
$K^+ K^- \eta'$ | $a \to K^+ K^{*-}[\to K^-\,\eta']$ : $\frac{\sqrt{6} \left(c_{s} - c_{u}\right) \sin{\left(\theta_{P} \right)}}{2}$; $a \to K^- K^{*+}[\to K^+\,\eta']$ : $\frac{\sqrt{6} \left(c_{s} - c_{u}\right) \sin{\left(\theta_{P} \right)}}{2}$
$K^0 \bar{K}^0 \eta$ | $a \to K^0 \bar{K}^{*0}[\to \bar{K}^0\,\eta]$ : $- \frac{\sqrt{6} \left(c_{d} - c_{s}\right) \cos{\left(\theta_{P} \right)}}{2}$; $a \to \bar{K}^0 K^{*0}[\to K^0\,\eta]$ : $- \frac{\sqrt{6} \left(c_{d} - c_{s}\right) \cos{\left(\theta_{P} \right)}}{2}$
$K^0 \bar{K}^0 \eta'$ | $a \to K^0 \bar{K}^{*0}[\to \bar{K}^0\,\eta']$ : $- \frac{\sqrt{6} \left(c_{d} - c_{s}\right) \sin{\left(\theta_{P} \right)}}{2}$; $a \to \bar{K}^0 K^{*0}[\to K^0\,\eta']$ : $- \frac{\sqrt{6} \left(c_{d} - c_{s}\right) \sin{\left(\theta_{P} \right)}}{2}$
$\eta \eta \eta$ | --
$\eta \eta \eta'$ | --
$\eta \eta' \eta'$ | --
$\eta' \eta' \eta'$ | --

[C_A=c_sd E_sd]  required nonzero final S=1
channel | nonzero spectator -> Vbar-pair : trace product
-----------|------------
$\pi^+ \pi^- K^0$ | $a \to \pi^- K^{*+}[\to \pi^+\,K^0]$ : $- c_{sd}$; $a \to K^0 \rho^0[\to \pi^+\,\pi^-]$ : $- c_{sd}$
$\pi^- \pi^0 K^+$ | $a \to \pi^- K^{*+}[\to \pi^0\,K^+]$ : $- \frac{\sqrt{2} c_{sd}}{2}$; $a \to \pi^0 K^{*0}[\to \pi^-\,K^+]$ : $\frac{\sqrt{2} c_{sd}}{2}$; $a \to K^+ \rho^-[\to \pi^-\,\pi^0]$ : $\sqrt{2} c_{sd}$
$\pi^- K^+ \eta$ | $a \to \pi^- K^{*+}[\to K^+\,\eta]$ : $\frac{\sqrt{6} c_{sd} \cos{\left(\theta_{P} \right)}}{2}$; $a \to \eta K^{*0}[\to \pi^-\,K^+]$ : $- \frac{\sqrt{6} c_{sd} \cos{\left(\theta_{P} \right)}}{2}$
$\pi^- K^+ \eta'$ | $a \to \pi^- K^{*+}[\to K^+\,\eta']$ : $\frac{\sqrt{6} c_{sd} \sin{\left(\theta_{P} \right)}}{2}$; $a \to \eta' K^{*0}[\to \pi^-\,K^+]$ : $- \frac{\sqrt{6} c_{sd} \sin{\left(\theta_{P} \right)}}{2}$
$\pi^0 \pi^0 K^0$ | $a \to \pi^0 (\times2)K^{*0}[\to \pi^0\,K^0]$ : $- \frac{c_{sd}}{2}$
$\pi^0 K^0 \eta$ | $a \to \pi^0 K^{*0}[\to K^0\,\eta]$ : $- \frac{\sqrt{3} c_{sd} \cos{\left(\theta_{P} \right)}}{2}$; $a \to \eta K^{*0}[\to \pi^0\,K^0]$ : $\frac{\sqrt{3} c_{sd} \cos{\left(\theta_{P} \right)}}{2}$
$\pi^0 K^0 \eta'$ | $a \to \pi^0 K^{*0}[\to K^0\,\eta']$ : $- \frac{\sqrt{3} c_{sd} \sin{\left(\theta_{P} \right)}}{2}$; $a \to \eta' K^{*0}[\to \pi^0\,K^0]$ : $\frac{\sqrt{3} c_{sd} \sin{\left(\theta_{P} \right)}}{2}$
$K^+ K^- K^0$ | $a \to K^+ \rho^-[\to K^-\,K^0]$ : $- c_{sd}$; $a \to K^0 \rho^0[\to K^+\,K^-]$ : $- \frac{c_{sd}}{2}$; $a \to K^0 \omega[\to K^+\,K^-]$ : $\frac{c_{sd}}{2}$; $a \to K^0 \phi[\to K^+\,K^-]$ : $c_{sd}$
$K^0 K^0 \bar{K}^0$ | $a \to K^0 (\times2)\rho^0[\to K^0\,\bar{K}^0]$ : $\frac{c_{sd}}{2}$; $a \to K^0 (\times2)\omega[\to K^0\,\bar{K}^0]$ : $\frac{c_{sd}}{2}$; $a \to K^0 (\times2)\phi[\to K^0\,\bar{K}^0]$ : $c_{sd}$
$K^0 \eta \eta$ | $a \to \eta (\times2)K^{*0}[\to K^0\,\eta]$ : $\frac{3 c_{sd} \cos^{2}{\left(\theta_{P} \right)}}{2}$
$K^0 \eta \eta'$ | $a \to \eta K^{*0}[\to K^0\,\eta']$ : $\frac{3 c_{sd} \sin{\left(2 \theta_{P} \right)}}{4}$; $a \to \eta' K^{*0}[\to K^0\,\eta]$ : $\frac{3 c_{sd} \sin{\left(2 \theta_{P} \right)}}{4}$
$K^0 \eta' \eta'$ | $a \to \eta' (\times2)K^{*0}[\to K^0\,\eta']$ : $\frac{3 c_{sd} \sin^{2}{\left(\theta_{P} \right)}}{2}$

Assertions passed: channel Q/S, both-vertex flavor flow, conjugate-vector maps, and all trace products.
Approximation: ideal omega-phi mixing. Exact symbolic inputs: cu, cd, cs, c_sd and theta_P.
S=-1 classification is obtained by Hermitian conjugating fields and E_sd -> E_ds.
