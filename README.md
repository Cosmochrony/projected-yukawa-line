# Toward Projected Yukawa Couplings: A Conditional Weak Determinant Line and Generation-Level Ordering

J. Beau, Independent Researcher, France

## Status

Working paper, v3.0. DOI: [10.5281/zenodo.20767265](https://doi.org/10.5281/zenodo.20767265)

## Summary

Companion note of the Cosmochrony fermionic-matter sub-programme. It opens the mass-sector frontier and records two
statements, each with the hypothesis it needs. It does **not** derive masses, Yukawa couplings, or mixing.

## Result

1. **Determinant line.** Algebra, proved on supplied model data: on a rank-two bundle with structure group in
   $\mathrm{SL}(2,\mathbb{C})$, $\wedge^2$ is the unique one-dimensional functor among $\mathrm{Sym}^k$ and $\wedge^k$
   ($k \ge 1$; the powers $\det^b$ are also one-dimensional and lie outside this family), and it is trivial, hence
   invisible to the adjoint sector $\mathrm{Sym}^2$. Under the spin solder
   [H-Spin] of Q14, $\wedge^2(S_L)$ is therefore a trivial line and carries no hypercharge. Under the distinct weak
   factor [H-Weak] of Q14, which is supplied by no source and not constructed, $L_Y := \wedge^2(E_{\mathrm{weak}})$ with
   structure group $\mathrm{U}(2)$ carries the abelian weight: within this functor route, a conditional weak
   determinant line, not a line of the spinor bundle.

2. **Generation-level ordering.** On the gauge-singlet triplet $\mathbb{C}^3_{\mathrm{gen}} = \mathrm{span}(e_0, e_+,
   e_-)$, the model operator $\mathrm{diag}(1, \tfrac12 + u, \tfrac12 - u)$ has exit deficits
   $\{0, \tfrac12 - u, \tfrac12 + u\}$ that order three levels for $0 < u < \tfrac12$, with a level ratio independent of
   the spectral scale.
   Reading the model operator as the restriction of $E_\Pi^2$ is a named hypothesis, [H-Res], which is not supplied;
   the map from levels to generations is not established. The even sector $\mathrm{diag}(1, \tfrac12, \tfrac12)$ is the
   algebraic value of $(C_2 - J_3^2)/C_2$ at $C_2 = 2$; reading it as the Born--Infeld even sector would require an
   identification with the conditional $3\times3$ model of O30 that is not available and is not used.

## Method

Under [H-Weak] $L_Y$ supplies the line from which the weights of the coupling are built (the exponents $k$, $m$ are not
fixed by [H-Weak] alone; Q14 selects them conditionally under both hypotheses, a supplied colour module, anomaly
constraints and a minimal normalisation of $L_Y$).
The model operator on $\mathbb{C}^3_{\mathrm{gen}}$ carries the levels (read as levels of $E_\Pi^2$ under [H-Res]);
the mass comes after.
The projected Yukawa operator needs a further carrier linking the weak doublet to the singlets, which is a missing
element, and a Lorentz-invariant coupling between the chiralities: a linear Lorentz-equivariant map is ruled out under
[H-Spin], a sesquilinear pairing or an antilinear parity lift is not, and none is supplied here.
Its norm, the sign of $u$, and the mixing are downstream open data.
A complex metaplectic phase is ruled out as a source of $R_{\mathrm{mix}}$ within the model of Q14;
the internal block $e_0 \leftrightarrow e_\pm$ is not covered by that remark and is treated in the companion note PYO.
No mass value is claimed.

## Interpretive outlook (a reading, not a result)

Where the Yukawa coupling attaches to hypercharge becomes the question of where a distinct weak factor comes from; the
generation levels are a dimensionless pattern whose reading as a mass hierarchy awaits the missing identifications.

## Anchors

Q14 (Theorem 2.1, [H-Spin], [H-Weak], exit-deficit dictionary), PRS (normal form on the supplied triplet), BIM
(radial factor of the split). See
the bibliography in `tex/cosmochrony-bibliography.bib`.

## Build

```bash
bash compile.sh   # -> out/ProjectedYukawaLine.pdf
```

## Audit

`code/front3_yukawa_line.py` (exact symbolic, 23 checks on the algebra and the model operator, including the region
$0 < u < \tfrac12$; it does not test [H-Spin], [H-Weak], [H-Res], or the level-to-generation map):
`python3 code/front3_yukawa_line.py`.

`code/` also holds four scripts that audit statements of the companion note PYO (not used by a statement of this note):
`front3c_polar_class_audit.py` (polar class under the rephasings that commute with each carrier's $J_3$),
`front3c_polar_class_nontriviality.py` (transverse part of the generator), `front3d_transverse_route.py` and
`front3e_reality_structure.py` (the $\mathfrak{sl}_2(\mathbb{C})$ step under the antilinear parity $J_\Pi$ of Q14).
Dependency: `sympy` (`pip install -r code/requirements.txt`).
