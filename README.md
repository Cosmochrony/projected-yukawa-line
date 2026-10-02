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
   [H-Spin] of Q14, $\wedge^2(S_L)$ is therefore a trivial line and carries no abelian weight of the spin group.
   Under the distinct weak factor [H-Weak] of Q14, which is supplied by no source and not constructed,
   $L_Y := \wedge^2(E_{\mathrm{weak}})$ with structure group $\mathrm{U}(2)$ carries the abelian weight: within this
   functor route, a conditional weak determinant line, not a line of the spinor bundle.

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
The projected Yukawa operator needs a weak linking carrier $K$ (a $\mathrm{U}(2)$-module linking the weak doublet to the
right sector, supplied by no source; the right fermion $P_R S \otimes L_Y^m$ has no factor $E_{\mathrm{weak}}$ and is a
$\mathrm{U}(2)$-character, a weak singlet carrying the hypercharge twist $L_Y^m$), which is the only missing element for
the existence of an invariant coupling (given [H-Spin] and [H-Weak]): the
Lorentz-invariant sesquilinear pairing exists (it contains the Lorentz scalar once, PYO Section 2), while a linear
Lorentz-equivariant map is ruled out under
[H-Spin] (non-zero, by Schur). [H-Fac] of PYO states two independent conditions: (F1) such a $K$ for which the space
$I_K$ of $K$-valued Lorentz-invariant sesquilinear $\mathrm{U}(2)$-invariant couplings is non-zero, and (F2) $\dim I_K
\le 1$; it
is not used here.
Its norm, the sign of $u$, and the mixing are downstream open data.
A complex metaplectic phase is ruled out as a source of the external block $R_{\mathrm{mix}}$ within the model of Q14,
which the image of $\mathfrak{sl}_2(\mathbb{C})$ does not reach (Q14 Remark 6.4).
For the internal block $e_0 \leftrightarrow e_\pm$: $A_\Pi$, the anti-Hermitian part of the $J_\Pi$-odd part of the
$\mathfrak{sl}_2$ lift of the step generator, vanishes identically for all complex coefficients, because with the
internal antilinear parity $J_\Pi$ of Q14 Section 6 (acting on the generation copy, not a spinor-level object) the
$J_\Pi$-odd part of the lift is its Hermitian part; this is a statement about $A_\Pi$ as defined ($A_\Pi$ is the
generator defined in PYO, not the anomaly density of Q14). The $J_\Pi$-even anti-Hermitian part has the non-zero
internal entry $\tfrac{\sqrt2}{2}(q - \bar p)$ (non-zero for real $p \ne q$), so the image of $\mathfrak{sl}_2$ does
reach the internal block through it. No physical exclusion is claimed: whether $A_\Pi$ is the right object is a
modelling choice of the companion note PYO that no source justifies, and Q14 excludes the internal block only by
hypothesis (Proposition 6.3 (i)-(ii)).
No mass value is claimed.

## Interpretive outlook (a reading, not a result)

Where the Yukawa coupling attaches to hypercharge becomes the question of where a distinct weak factor comes from; the
generation levels are a dimensionless pattern whose reading as a mass hierarchy awaits the missing identifications.

## Anchors

Q14 (Theorem 2.1, [H-Spin], [H-Weak], exit-deficit dictionary), PRS (normal form on the supplied triplet), BIM
(which does not predict $|u|$). See
the bibliography in `tex/cosmochrony-bibliography.bib`.

## Build

```bash
bash compile.sh   # -> out/ProjectedYukawaLine.pdf
```

## Audit

`code/front3_yukawa_line.py` (exact symbolic, 23 checks on the algebra and the model operator, including the region
$0 < u < \tfrac12$; it does not test [H-Spin], [H-Weak], [H-Res], or the level-to-generation map):
`python3 code/front3_yukawa_line.py`.

`code/` also holds four audit scripts of statements of the companion note PYO (not used by a statement of this
note; PYO's own scripts, which reproduce PYO, are in the `code/` directory of the PYO repository):
`front3c_polar_class_audit.py` (polar class under the rephasings that commute with each carrier's $J_3$),
`front3c_polar_class_nontriviality.py` (transverse part of the generator), `front3d_transverse_route.py` and
`front3e_reality_structure.py` (the $\mathfrak{sl}_2(\mathbb{C})$ step under the internal antilinear parity $J_\Pi$ of
Q14 Section 6).
Dependency: `sympy` (`pip install -r code/requirements.txt`).
