"""Front 3: determinant-line algebra and the level ordering of a model operator.

Exact symbolic verification (no sampling, except the labelled finite grid in check F) of the algebraic content of the
two statements opening the mass-sector frontier. Nothing here is new representation theory. The tensor identities are
those of Q14 Theorem 2.1(a) (Sym^2(C^2) = sl_2(C), wedge^2(C^2) = determinant line), the model operator is that of PRS
(diag(1, 1/2+u, 1/2-u) on C^3_gen = Sym^2(V_gen)), and the exit deficits are those of Q14 Def. 6.6 and Prop. 6.7.
Method: under [H-Weak] the determinant line L_Y := wedge^2(E_weak) supplies the line from which the weights of the
coupling are built (the exponents k, m are not fixed by the hypothesis); the model operator carries the levels, read as
levels of E_Pi^2 only under the named identification [H-Res]; the mass comes after (Yukawa norm + sign of u + mixing).

Typing (Q14): the rank-two fibre is either the spinor factor S_L, with structure group SL(2,C) under [H-Spin], on which
wedge^2 is a trivial line with no abelian weight of the spin group, or the weak factor E_weak of [H-Weak], with
structure group U(2), on
which L_Y = wedge^2(E_weak) carries the abelian weight. E_weak is not constructed in Q14. The script works on C^2 and
does not construct either bundle.

Seven groups of checks are run.
  (A) Dimensions, computed from the defining subspaces of C^2 (x) C^2: symmetric tensors (dim 3 = dim sl_2, computed
      as the traceless 2x2 matrices) and antisymmetric tensors (dim 1). The uniqueness of wedge^2 among Sym^k,
      wedge^k (k >= 1) rests on the dimension count in the paper's proof; the script does not enumerate the
      functors. The powers det^b are also one-dimensional and are not in that family.
  (B) Action on the line: the action of g (x) g on the antisymmetric tensor e1^e2 is computed explicitly and equals
      det(g); on the SL(2,C) slice and on full SU(2) (general element, reduced by x^2+y^2+z^2+t^2 = 1) it is 1, so
      wedge^2 is trivial there (any subgroup of SL(2,C)). Negative controls: g = 2 * 1 acts by 4, and the SU(2)
      remainder test fails when the unit-norm relation is not imposed.
  (C) Weight of a gl_2 generator: the derivative action Y (x) 1 + 1 (x) Y on the antisymmetric tensor is the trace
      y1 + y2. This generator lies in gl_2 and not in sl_2; it belongs to the structure algebra only if the structure
      group is U(2), which is the choice of [H-Weak] for E_weak. Checks B and C are separate and the script does NOT
      test that a generator belongs to the structure group of S_L or of E_weak: that membership is an input.
  (D) Normal form of the model operator: (C2 - J3^2)/C2 + u J3 (v = 0, CP-even) with C2 = 2, J3 = diag(0,1,-1) gives
      exactly diag(1, 1/2+u, 1/2-u) on (e_0, e_+, e_-).
  (E) Even sector: at u = 0 this reduces to diag(1, 1/2, 1/2), the algebraic value of (C2 - J3^2)/C2 at C2 = 2
      (PRS eq:even); its Born-Infeld reading through O30 is not used.
  (F) Exit deficits d(e_i) = s_0 - s_i, read from the diagonal entries of the model operator of (D) (s_0 is the e_0
      level, d(e_0) = 0 by definition): d(e_+) = 1/2 - u, d(e_-) = 1/2 + u. The script solves the system
      {0 < d(e_+) < d(e_-)} and checks that its solution set is exactly 0 < u < 1/2; that positivity of the two outer
      levels adds no condition (it is implied: s_- = d(e_+), s_+ = d(e_-)); that e_0 has the smallest deficit
      exactly for |u| <= 1/2; a finite exact grid of rational u inside and outside the interval (u = 7/10 and
      u = 1/2 violate it); and the identity d(e_-) - d(e_+) = 2u, kept as a labelled identity. The labelling
      e_+/e_- (J_3 = +1/-1) is the convention of Q14 Sec. 6.
  (G) Scale independence of the level RATIO: (1/2 + u)/(1/2 - u), computed from the rescaled operator lambda * E, is
      invariant under rescaling (an identity, labelled; the level difference is not, as a negative control), so the
      dimensionless level structure is fixed by u alone, while the absolute mass normalisation is NOT fixed here.

Not tested: the identification [H-Res] of the model operator with the restriction of E_Pi^2 to C^3_gen, the map from
levels to generations, [H-Spin], [H-Weak], and the Yukawa operator. No number is produced for any mass. No figures.
English.
"""

import sympy as sp


def main():
    checks = {}

    # ---- (A) dimensions from the defining subspaces of C^2 (x) C^2 ----------------------------------------
    T = sp.Matrix(2, 2, sp.symbols("t0:4"))           # a 2-tensor on C^2 (x) C^2 as a 2x2 matrix
    t_syms = list(T)
    sym_eqs = list(T - T.T)                            # symmetric tensors: T = T^T
    asym_eqs = list(T + T.T)                           # antisymmetric tensors: T = -T^T
    dim_sym2 = 4 - sp.Matrix([[sp.diff(e, v) for v in t_syms] for e in sym_eqs]).rank()
    dim_wedge2 = 4 - sp.Matrix([[sp.diff(e, v) for v in t_syms] for e in asym_eqs]).rank()
    dim_sl2 = 4 - sp.Matrix([[sp.diff(T.trace(), v) for v in t_syms]]).rank()   # traceless 2x2 matrices
    checks["A_dim_sym2_is_3"] = (dim_sym2 == 3)
    checks["A_dim_wedge2_is_1"] = (dim_wedge2 == 1)
    checks["A_sym2_dim_equals_dim_sl2_computed"] = (dim_sym2 == dim_sl2)   # dim sl_2 computed as traceless matrices

    # ---- (B) action of g (x) g on the antisymmetric tensor e1^e2 is det g ---------------------------------------
    a, b, c, d = sp.symbols("a b c d")
    g = sp.Matrix([[a, b], [c, d]])
    e1, e2 = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    w = sp.kronecker_product(e1, e2) - sp.kronecker_product(e2, e1)       # e1 ^ e2 in C^2 (x) C^2
    gw = sp.kronecker_product(g, g) * w
    checks["B_wedge_action_is_det"] = sp.simplify(gw - g.det() * w) == sp.zeros(4, 1)
    # SL(2,C) slice: d = (1 + b c)/a gives det g = 1, hence g (x) g acts trivially on w
    g_sl = g.subs(d, (1 + b * c) / a)
    gw_sl = sp.simplify(sp.kronecker_product(g_sl, g_sl) * w - w)
    checks["B_trivial_on_SL2"] = gw_sl == sp.zeros(4, 1)
    # full SU(2): g = [[x + i y, z + i t], [-(z - i t), x - i y]] with x^2 + y^2 + z^2 + t^2 = 1 (real x, y, z, t);
    # reduce g (x) g w - w by the relation (polynomial remainder in x)
    x1, y1_, z1, t1 = sp.symbols("x1 y1 z1 t1", real=True)
    gsu2 = sp.Matrix([[x1 + sp.I * y1_, z1 + sp.I * t1], [-(z1 - sp.I * t1), x1 - sp.I * y1_]])
    diff_su2 = sp.expand(sp.kronecker_product(gsu2, gsu2) * w - w)
    reduced = [sp.rem(sp.expand(e), x1**2 - (1 - y1_**2 - z1**2 - t1**2), x1) for e in diff_su2]
    checks["B_trivial_on_full_su2"] = all(sp.simplify(e) == 0 for e in reduced)
    # the same remainder test applied to a non-unit-determinant element (relation violated) does not vanish
    bad = sp.expand(sp.kronecker_product(gsu2, gsu2) * w - w)
    checks["B_su2_relation_needed"] = not all(sp.simplify(e) == 0 for e in bad)
    # contrast: off SL(2), e.g. g = 2 * identity acts by det = 4 != 1 (the check can fail)
    g2 = 2 * sp.eye(2)
    checks["B_nontrivial_off_SL2"] = sp.simplify(sp.kronecker_product(g2, g2) * w - 4 * w) == sp.zeros(4, 1)

    # ---- (C) weight of a gl_2 generator: derivative action on e1^e2 is the trace --------------------------------
    y1, y2 = sp.symbols("y1 y2")
    Y = sp.diag(y1, y2)                             # diagonal gl_2 generator on C^2 (not in sl_2 unless y1 + y2 = 0)
    Yw = (sp.kronecker_product(Y, sp.eye(2)) + sp.kronecker_product(sp.eye(2), Y)) * w
    checks["C_gl2_weight_is_trace"] = sp.simplify(Yw - (y1 + y2) * w) == sp.zeros(4, 1)
    checks["C_traceless_generator_weight_zero"] = sp.simplify(Yw.subs(y2, -y1)) == sp.zeros(4, 1)

    # ---- (D) normal form of the model operator diag(1, 1/2+u, 1/2-u) on C^3_gen ----------------------------
    u, v = sp.symbols("u v", real=True)
    C2 = sp.Integer(2)
    J3 = sp.diag(0, 1, -1)
    Rmix = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])  # R_mix of Q14 Sec. 6 (antisymmetric on e_+, e_-)
    Epi2 = (C2 * sp.eye(3) - J3**2) / C2 + u * J3 + v * Rmix
    Epi2_even = Epi2.subs({u: 0, v: 0})
    checks["D_normalform_diag"] = sp.simplify(Epi2.subs(v, 0) - sp.diag(1, sp.Rational(1, 2) + u,
                                                                       sp.Rational(1, 2) - u)) == sp.zeros(3)
    # ---- (E) algebraic even sector ---------------------------------------------------------------
    checks["E_even_algebraic"] = sp.simplify(Epi2_even - sp.diag(1, sp.Rational(1, 2), sp.Rational(1, 2))) \
        == sp.zeros(3)

    # ---- (F) exit deficits computed from the model operator, and their ordering ---------------------------------
    Emod = Epi2.subs(v, 0)                                   # the model operator of (D), read entry by entry
    s0, s_plus, s_minus = Emod[0, 0], Emod[1, 1], Emod[2, 2]  # levels on (e_0, e_+, e_-); s_0 is the e_0 level
    dplus, dminus = s0 - s_plus, s0 - s_minus                 # exit deficits d(e_i) = s_0 - s_i  (Q14 Def. 6.6)
    checks["F_deficit_eplus_from_operator"] = sp.simplify(dplus - (sp.Rational(1, 2) - u)) == 0
    checks["F_deficit_eminus_from_operator"] = sp.simplify(dminus - (sp.Rational(1, 2) + u)) == 0
    checks["F_deficit_difference_identity"] = sp.simplify(dminus - dplus - 2 * u) == 0   # an identity, kept as such
    ur = sp.symbols("ur", real=True)
    order_conds = [sp.Rational(0) < dplus.subs(u, ur), dplus.subs(u, ur) < dminus.subs(u, ur)]
    outer_conds = [s_plus.subs(u, ur) > 0, s_minus.subs(u, ur) > 0]
    region = sp.reduce_inequalities(order_conds, ur)
    open_interval = sp.Interval.open(0, sp.Rational(1, 2))
    # the solution set of {0 < d(e_+) < d(e_-)} is exactly the open interval 0 < u < 1/2
    checks["F_region_is_0_lt_u_lt_half"] = sp.simplify(region.as_set() - open_interval) == sp.EmptySet \
        and sp.simplify(open_interval - region.as_set()) == sp.EmptySet
    # positivity of the two outer levels is IMPLIED by the ordering (s_- = d(e_+) and s_+ = d(e_-)), not an extra test
    region_full = sp.reduce_inequalities(order_conds + outer_conds, ur)
    checks["F_outer_levels_positive_implied_by_ordering"] = sp.simplify(
        region_full.as_set() - region.as_set()) == sp.EmptySet and sp.simplify(
        region.as_set() - region_full.as_set()) == sp.EmptySet
    # e_0 has zero deficit by definition; it is the smallest deficit exactly when both deficits are >= 0,
    # that is |u| <= 1/2
    region_e0 = sp.reduce_inequalities([dplus.subs(u, ur) >= 0, dminus.subs(u, ur) >= 0], ur)
    checks["F_e0_smallest_deficit_iff_abs_u_le_half"] = sp.simplify(
        region_e0.as_set() - sp.Interval(-sp.Rational(1, 2), sp.Rational(1, 2))) == sp.EmptySet and sp.simplify(
        sp.Interval(-sp.Rational(1, 2), sp.Rational(1, 2)) - region_e0.as_set()) == sp.EmptySet

    def ordered(uv):
        """Exact test of 0 < d(e_+) < d(e_-) at a rational u."""
        return all(bool(c_.subs(ur, uv)) for c_ in order_conds)

    inside = [sp.Rational(k, 100) for k in range(1, 50)]
    checks["F_grid_inside_all_ordered"] = all(ordered(uv) for uv in inside)
    outside = [sp.Rational(-1, 10), sp.Integer(0), sp.Rational(1, 2), sp.Rational(7, 10), sp.Integer(1)]
    checks["F_grid_outside_none_ordered"] = not any(ordered(uv) for uv in outside)
    # at u = 7/10 one has d(e_+) < 0 = d(e_0): e_0 is not the smallest deficit
    checks["F_bound_needed_example"] = bool(dplus.subs(u, sp.Rational(7, 10)) < 0)

    # ---- (G) scale-independence of the level ratio, computed from the scaled operator --------------------------
    lam = sp.symbols("lambda", positive=True)
    Escaled = lam * Emod
    ratio = (Emod[1, 1]) / (Emod[2, 2])
    ratio_scaled = (Escaled[1, 1]) / (Escaled[2, 2])
    checks["G_ratio_scale_invariant_identity"] = sp.simplify(ratio_scaled - ratio) == 0   # an identity, labelled
    # negative control: the level DIFFERENCE is not scale invariant
    checks["G_level_difference_not_scale_invariant"] = sp.simplify(
        (Escaled[1, 1] - Escaled[2, 2]) - (Emod[1, 1] - Emod[2, 2])) != 0

    # ---------------------------------------------------------------------------------------------
    print("Front 3 - determinant-line algebra and level ordering of a model operator (exact symbolic)")
    print("=" * 90)
    print("  (A) dim Sym^2(C^2) = 3 (= dim sl_2);  dim wedge^2(C^2) = 1 (determinant line)")
    print("  (B) action of g (x) g on wedge^2 = det g, computed; = 1 on SL(2,C) and on full SU(2) => trivial line")
    print("  (C) gl_2 generator diag(y1,y2) acts on wedge^2 by trace y1+y2 (needs structure group U(2), [H-Weak])")
    print("  (D) model operator (C2-J3^2)/C2 + u J3 = diag(1, 1/2+u, 1/2-u)   (v=0, CP-even)")
    print("  (E) u=0 reduces to diag(1, 1/2, 1/2)  (algebraic value of (C2-J3^2)/C2 at C2=2; BI reading not used)")
    print("  (F) exit deficits d(e_0)=0, d(e_+)=1/2-u, d(e_-)=1/2+u (from the operator): ordered iff 0 < u < 1/2")
    print("  (G) level ratio (1/2+u)/(1/2-u) scale-invariant (identity); absolute mass NOT fixed here")
    print("-" * 90)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 90)
    print(f"  checks run: {len(checks)}")
    print("METHOD: under [H-Weak] L_Y supplies the line from which weights are built (exponents k, m not fixed by")
    print("        [H-Weak] alone); the model operator carries the levels (read as E_Pi^2 only under [H-Res]); the")
    print("        mass comes after (Yukawa operator Y_Pi open: norm + sign(u) + mixing).")
    print("NOT TESTED: [H-Res], the level-to-generation map, [H-Spin], [H-Weak], the uniqueness enumeration of")
    print("            functors, membership of a generator in a structure group.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
