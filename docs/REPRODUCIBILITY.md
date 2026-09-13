# Reproducibility map

The default branch contains two frozen layers. The root-level modules preserve the v1.2 public baseline; `current/` adds the calculations used by the PRD v2.2.0 manuscript.

## Quick regression check

```bash
pip install -r requirements.txt
pytest -q
```

The root pytest configuration collects the baseline suite plus current source-only deterministic suites. The v1.6 branch-parity and v1.7 causal-kinetics full regressions additionally require frozen descendant-history tables from the submission-matched reviewer archive; those multi-megabyte intermediates are intentionally not duplicated in the GitHub default branch.

## Current PRD add-on suites

- `current/reproducibility_followup/` — nonlinear PM, matched ZA–PM, descendant anchoring, selector ceiling and threshold inversion.
- `current/research_v1.4/` — geometry-state identities, leaf-global field and scale-closure diagnostics.
- `current/research_v1.5/` — control-basis and constitutive-driver audits.
- `current/research_v1.6/` — one-sector/branch-parity and signed-history diagnostics.
- `current/research_v1.7/` — causal phase kinetics and passive-memory audits.
- `current/research_v1.8/` — retarded/Kubo response identification.
- `current/research_v1.9/` — nonlinear 1PI, thin-wall nucleation and local-KMS implementation checks.
- `current/research_v2.0/` — reciprocal/nonreciprocal phase–geometry feedback benchmarks.
- `current/research_v2.2/` — shared-control/global-competition and global-budget identifiability audit.

These modules reproduce deterministic identities and stress tests used in the PRD manuscript. They do not turn the low-resolution PM pilot into a precision cosmological simulation, and they do not supply a microscopic derivation of the remaining response coefficients.

## Historical v1.2 baseline

The exact submission snapshot is preserved by tag `v1.2.0`. The root-level baseline remains executable for continuity. The old documentation below is retained as a historical description; its original test-count wording should be read as applying to that development stage, not the current default branch.

---

# Reproducibility map

`python run_all.py` reproduces all analyses for which required data are present:

1. `analysis_phase.py` — normalized sextic coexistence/spinodal constants and Miyamoto–Nagai geometry factor.
2. `analysis_regional_state.py` — exact EdS identities for the regional coordinate and finite standard dynamical landmarks.
3. `analysis_excursion_phase.py` — exact two-absorbing-barrier first-passage fractions and the illustrative Eulerian volume weighting.
4. `analysis_power_clock.py` — Planck-normalized variance-clock audit and inverse early-selector target.
5. `analysis_curvature_selector.py` — exact regional matter-curvature-equality landmark, operational density/velocity estimator, and its fixed-scale first-passage benchmark.
6. `analysis_basin_selector.py` — maximal-active basin antichain audit for non-recursive domain accounting.
7. `analysis_conversion_sign.py` — monotone-conversion residual-EOS sign audit and external-growth stress test.
8. `analysis_basin_dynamics.py` — history-dependent basin-occupation bookkeeping and absorbing-activation limiting audit.
9. `analysis_interface_coarsening.py` — exact interface-energy EOS identity, hydraulic-scale benchmark and coarsening crossing threshold.
10. `analysis_phase_interface.py` — normalized sextic planar-wall tension and thickness constants.
11. `analysis_phase_front.py` — Allen-Cahn sharp-interface velocity, critical radius, and capillary-barrier audit.
12. `analysis_emg.py` — one-parameter SPARC RAR consistency diagnostic; skipped when `data/RAR.mrt` is absent.

The frozen v1.2 suite contains 30 deterministic regression tests. Numerical outputs are written under `results/`; manuscript figures are written under `figures/`.

## Third-party SPARC table

For a GitHub-safe checkout, obtain the public table with `python scripts/fetch_sparc.py`. The helper verifies the manuscript-snapshot hash before the fit is run.

## Scope of the excursion calculation

The excursion-set calculation is a structural prototype, not a precision cosmology likelihood. It uses a sharp-`k` Markov Gaussian first-passage model, standard spherical collapse/void barriers, and the variance clock `tau=D^2 S_Omega`. No redshift is inferred until a physical basin variance and parent fluctuation spectrum are specified. The original two-barrier void abundance model is also known to require corrections for accurate simulation abundances; the manuscript uses it only to replace arbitrary conversion-history shapes by an explicit standard-dynamics kernel.


## Scope of the curvature-equality calculation

The identity `Upsilon^2=1/Omega_m,local` and the spherical matter-curvature-equality landmark are exact within the EdS separate-universe prototype. The subsequent `phi_V(z)` curves use the same Markov first-passage bookkeeping and an external Planck-normalized LambdaCDM growth benchmark. They are a scale/timing consistency test, not a self-consistent background solution. Curvature domination itself is not identified with vacuum energy; it is only tested as a non-fitted regional event that could bias a separate V-phase order parameter.

## Scope of the basin-hierarchy calculation

The basin selector code is a combinatorial regression audit, not a cosmic-web finder. It verifies the manuscript proposition that maximal active nodes of a laminar merge-tree hierarchy form an antichain and therefore remove nested-domain double counting. A real application must reconstruct a density/Morse/watershed hierarchy, evaluate the regional expansion-density balance on its nodes, and evolve the history-dependent order parameter. In the static synthetic example the eligibility and activity predicates are identified only for testing the non-recursion logic.

## Scope of the monotone-conversion sign audit

The analytic inequality assumes a fixed positive intrinsic V-phase energy density, negligible interface energy, pressureless D+C bookkeeping, and a nondecreasing active V filling fraction. The reported `w_res(z=0)` values use the same external LCDM growth benchmark as the curvature-selector timing audit and are not a self-consistent cosmological likelihood. Their purpose is to expose a sign/timing constraint: matching the present filling fraction is not enough if conversion remains too rapid today.

## Scope of the history-dependent basin audit

The v0.8 basin-dynamics script is a deterministic bookkeeping test on a fixed synthetic laminar hierarchy, not an N-body result. It demonstrates that an instantaneous `Gamma >= 1` classifier can chatter near threshold while a phase label with strong hysteresis retains history. In the absorbing limiting case the active-node set grows by inclusion; maximal active nodes give the same geometric union without nested double counting. Therefore the active filling fraction is monotone for a fixed hierarchy/measure, activating the analytic `w_res <= -1` sign theorem. Real basin geometry evolves and nodes merge or disappear, so a simulation can evade this limiting corollary only through explicitly modeled deactivation/restructuring, interface physics, or nonconstant intrinsic phase energy.


## Scope of the interface-coarsening audit

The interface calculation is an exact effective bookkeeping identity once `phi_V(a)`, boundary area density `s_Sigma(a)`, positive interface tension `sigma_Sigma`, and intrinsic V-phase density `epsilon_V` are specified. The NEXUS wall area and catalogued void radii are used only as external morphology proxies to check the hydraulic scale `L_Sigma=3 phi_V/s_Sigma`; NEXUS walls are not assumed to be the physical phase boundary. The interval-averaged crossing estimate is therefore a scale/sign audit, not a dark-energy likelihood or a measurement of `sigma_Sigma`. A real test requires one simulation/catalogue pipeline that measures the maximal active-basin union and its boundary consistently at multiple redshifts.


## Scope of the sextic phase-interface audit

The wall calculation uses the same normalized sextic potential as the three-phase toy model and adds one standard positive square-gradient stiffness `kappa`. The reported coefficients are exact numerical quadratures of that normalized potential at coexistence. They do **not** determine the physical value of `kappa`, the dimensional phase free-energy scale `E0`, or its ratio to the intrinsic V-phase density `epsilon_V`. The comparison to Mpc-scale cosmic-web sheet widths is therefore conditional and only tests whether the interface route requires an obviously absurd mesoscopic scale.
