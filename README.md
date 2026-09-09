# Projected Yukawa Line and Generation-Level Assignment

*Functorial Closure of $L_Y = \wedge^2(\mathcal{S}_\Pi)$ and Spectral Ordering on $C^3_{\mathrm{gen}}$*

J. Beau, Independent Researcher, France

## Status

Working paper, v2.0. DOI: [10.5281/zenodo.20767265](https://doi.org/10.5281/zenodo.20767265)

## Summary

Companion note of the Cosmochrony fermionic-matter sub-programme. It opens the mass-sector frontier and closes only
its first stage: the functorial Yukawa coupling line and the generation-level ordering. It does **not** derive masses.

## Result

Two structural propositions, both verified by exact symbolic computation (no sampling).

1. **Projected Yukawa line.** On the rank-two admissible spinor bundle `S_Π`, the determinant line `L_Y = ∧²(S_Π)`
   is the unique functorial line of the carrier, invisible to the `SU(2)_L` adjoint sector `Sym²(S_Π) ≃ ad_C(P_Π)`
   (action `det = 1` on `SL(2,C)`), and the carrier of the abelian `U(1)_Y` weight (action by the trace). The projected
   Yukawa sector is therefore the determinant-line sector of the admissible spinor functorial closure, not an external
   bundle datum.

2. **Generation-level assignment.** On the gauge-singlet triplet `C³_gen = span(e₀, e₊, e₋)`, the squared projective
   residue takes the closed normal form `E_Π²|gen = diag(1, ½+u, ½−u)`, whose exit deficits `{0, ½−u, ½+u}` order the
   three generation levels: `e₀` lightest (central, non-exiting), `e₋` heaviest (for `u > 0`). The even sector
   `diag(1, ½, ½)` is the algebraic value of `(C₂ − J₃²)/C₂` at `C₂ = 2`; reading it as the Born–Infeld even sector
   would require an identification with the conditional `3×3` model of O30 that is not available and is not used.

## Method lock

`L_Y` fixes the coupling line; `E_Π²|gen` fixes the levels; the mass comes after. The projected Yukawa operator
`Y_Π : S_{L,Π} ⊗ L_Y → S_{R,Π}`, its norm, the sign of `u`, and the complex mixing phase are downstream open data. No
mass value is claimed.

## Anchors

Q14 (Theorem A, exit-deficit dictionary), PRS (the `E_Π²` normal form), A4-note (radial factor). See the bibliography in `tex/cosmochrony-bibliography.bib`.

## Build

```bash
bash compile.sh   # -> out/ProjectedYukawaLine.pdf
```

## Audit

`code/front3_yukawa_line.py` (exact symbolic, 15 checks): `python3 code/front3_yukawa_line.py`.
