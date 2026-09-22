The scope of the project is to re-derive the expressions for ALP process interacting with light mesons by using a consistent framework that we have developed in a previous work.

Our previous work is not yet publicly available, but the aspects relevant for mesonic interactions are firm. The source is available in `/memory/sources/Alda:2026bar.pdf`. Use the temporary identifier `Alda:2026bar` for all references to it.

Other related works, that you might want to consult, but being aware of different assumptions, conventions, and in some cases even mistakes, are:

- arXiv:1811.03474
- arXiv:2506.15637
- arXiv:2501.04525
- arXiv:2609.16104
- arXiv:2505.24822

Some steps to implement:

- Systematic derivation of all a-> 3 mesons and meson-> a 2 mesons amplitudes in pure chiral perturbation theory
- Implementation of efficient routines for phase space numerical integration of decays into three bodies using Mandelstam variables and optimized for previously obtained amplitudes
- Add heavier mesons and resonances. In the case of vector and axial-vector mesons, evaluate different theoretical models (resonance chiPT, hidden symmetries, background fields for WZW) and make an informed decision.
- When needed, implement specific forms of meson propagators and S-wave amplitudes
- For resonant contributions to a-> 3 mesons and meson -> a 2 mesons, find systematic ways to introduce these contributions in the total amplitude
- Some unstable resonances will result in a->4 mesons. Implement efficient routines for phase space numerical integration using Mandelstam variables tailored for this case
- For processes into 5 or more final states, produce reasonable estimations that do not require numerical integration
- Calculate amplitudes for a-> 4 mesons resulting from the 5-dimensional topological term of WZW