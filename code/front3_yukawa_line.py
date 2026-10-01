"""Front 3: determinant-line algebra and the level ordering of a model operator.

Exact symbolic verification (no sampling, except the labelled finite grid in check F) of the algebraic content of the
two statements opening the mass-sector frontier. Nothing here is new representation theory. The tensor identities are
those of Q14 Theorem 2.1(a) (Sym^2(C^2) = sl_2(C), wedge^2(C^2) = determinant line), the model operator is that of PRS
(diag(1, 1/2+u, 1/2-u) on C^3_gen = Sym^2(V_gen)), and the exit deficits are those of Q14 Prop. 6.7.
Method: under [H-Weak] the determinant line L_Y := wedge^2(E_weak) supplies the line from which the weights of the
coupling are built (the exponents k, m are not fixed by the hypothesis); the model operator carries the levels, read as
levels of E_Pi^2 only under the named identification [H-Res]; the mass comes after (Yukawa norm + sign of u + mixing).

Typing (Q14): the rank-two fibre is either the spinor factor S_L, with structure group SL(2,C) under [H-Spin], on which
wedge^2 is a trivial line with no hypercharge, or the weak factor E_weak of [H-Weak], with structure group U(2), on
which L_Y = wedge^2(E_weak) carries the abelian weight. E_weak is not constructed in Q14. The script works on C^2 and
does not construct either bundle.

Seven groups of checks are run.
  (A) Dimensions, computed from the defining subspaces of C^2 (x) C^2: symmetric tensors (dim 3 = dim sl_2, computed
      as the traceless 2x2 matrices) and antisymmetric tensors (dim 1). The uniqueness of wedge^2 among Sym^k,
      wedge^k (k >= 1) rests on the dimension count in the paper's proof; the script does not enumerate the
      functors. The powers det^b are also one-dimensional and are not in that family.
  (B) Action on the line: the action of g (x) g on the antisymmetric tensor e1^e2 is computed explicitly and equals
      det(g); on the sl_2 and su_2 slices it is 1, so wedge^2 is trivial there (any subgroup of SL(2,C)).
  (C) Weight of a gl_2 generator: the derivative action Y (x) 1 + 1 (x) Y on the antisymmetric tensor is the trace
      y1 + y2. This generator lies in gl_2 and not in sl_2; it belongs to the structure algebra only if the structure
      group is U(2), which is the choice of [H-Weak] for E_weak. Checks B and C are separate and the script does NOT
      test that a generator belongs to the structure group of S_L or of E_weak: that membership is an input.
  (D) Normal form of the model operator: (C2 - J3^2)/C2 + u J3 (v = 0, CP-even) with C2 = 2, J3 = diag(0,1,-1) gives
      exactly diag(1, 1/2+u, 1/2-u) on (e_0, e_+, e_-).
  (E) Even sector: at u = 0 this reduces to diag(1, 1/2, 1/2), the algebraic value of (C2 - J3^2)/C2 at C2 = 2
      (PRS eq:even); its Born-Infeld reading through O30 is not used.
  (F) Exit deficits d(e_i) = s_0 - s_i with s_0 = 1: d(e_0) = 0, d(e_+) = 1/2 - u, d(e_-) = 1/2 + u. The script solves
      the system {0 < d(e_+), d(e_+) < d(e_-), both outer levels > 0} and checks that its solution set is exactly
      0 < u < 1/2; it checks a finite exact grid of rational u inside and outside that interval (u = 7/10 and u = 1/2
      violate it); it also keeps the identity d(e_-) - d(e_+) = 2u. The labelling e_+/e_- (J_3 = +1/-1) is the
      convention of Q14 Sec. 6.
  (G) Scale independence of the level RATIO: (1/2 + u)/(1/2 - u) is invariant under rescaling, so the dimensionless
      level structure is fixed by u alone, while the absolute mass normalisation is NOT fixed here.

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
    # real-section SU(2) slice [[al, be], [-be, al]] with al^2 + be^2 = 1: reduce by the relation (polynomial remainder)
    al, be = sp.symbols("alpha beta", real=True)
    gsu2 = sp.Matrix([[al, be], [-be, al]])
    diff_su2 = sp.kronecker_product(gsu2, gsu2) * w - w
    reduced = [sp.rem(sp.expand(e), be**2 - (1 - al**2), be) for e in diff_su2]
    checks["B_trivial_on_su2_slice"] = all(sp.simplify(e) == 0 for e in reduced)
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
    Rmix = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])   # off-diagonal generation mixing
    Epi2 = (C2 * sp.eye(3) - J3**2) / C2 + u * J3 + v * Rmix
    Epi2_even = Epi2.subs({u: 0, v: 0})
    checks["D_normalform_diag"] = sp.simplify(Epi2.subs(v, 0) - sp.diag(1, sp.Rational(1, 2) + u,
                                                                       sp.Rational(1, 2) - u)) == sp.zeros(3)
    # ---- (E) algebraic even sector ---------------------------------------------------------------
    checks["E_even_algebraic"] = sp.simplify(Epi2_even - sp.diag(1, sp.Rational(1, 2), sp.Rational(1, 2))) \
        == sp.zeros(3)

    # ---- (F) exit-deficit ordering and its region of validity -------------------------------------
    s0, sp_, sm = 1, sp.Rational(1, 2) + u, sp.Rational(1, 2) - u
    d0, dplus, dminus = s0 - s0, s0 - sp_, s0 - sm
    checks["F_deficit_e0_zero"] = (d0 == 0)
    checks["F_deficit_eplus"] = sp.simplify(dplus - (sp.Rational(1, 2) - u)) == 0
    checks["F_deficit_eminus"] = sp.simplify(dminus - (sp.Rational(1, 2) + u)) == 0
    checks["F_deficit_difference_identity"] = sp.simplify(dminus - dplus - 2 * u) == 0   # an identity, kept as such
    ur = sp.symbols("ur", real=True)
    conds = [sp.Rational(0) < dplus.subs(u, ur), dplus.subs(u, ur) < dminus.subs(u, ur),
             sp_.subs(u, ur) > 0, sm.subs(u, ur) > 0]
    region = sp.reduce_inequalities(conds, ur)
    # the solution set of {0 < d(e_+) < d(e_-), both outer levels > 0} is exactly the open interval 0 < u < 1/2
    checks["F_region_is_0_lt_u_lt_half"] = sp.simplify(region.as_set() - sp.Interval.open(0, sp.Rational(1, 2))) \
        == sp.EmptySet and sp.simplify(sp.Interval.open(0, sp.Rational(1, 2)) - region.as_set()) == sp.EmptySet

    def ordered(uv):
        """Exact test of 0 < d(e_+) < d(e_-) with both outer levels strictly positive at a rational u."""
        return all(bool(c_.subs(ur, uv)) for c_ in conds)

    inside = [sp.Rational(k, 100) for k in range(1, 50)]
    checks["F_grid_inside_all_ordered"] = all(ordered(uv) for uv in inside)
    outside = [sp.Rational(-1, 10), sp.Integer(0), sp.Rational(1, 2), sp.Rational(7, 10), sp.Integer(1)]
    checks["F_grid_outside_none_ordered"] = not any(ordered(uv) for uv in outside)
    # at u = 7/10 the central level is no longer the smallest deficit: d(e_+) < 0 = d(e_0)
    checks["F_bound_needed_example"] = bool(dplus.subs(u, sp.Rational(7, 10)) < 0)

    # ---- (G) scale-independence of the level ratio ----------------------------------------------
    lam = sp.symbols("lambda", positive=True)
    ratio = (sp.Rational(1, 2) + u) / (sp.Rational(1, 2) - u)
    ratio_scaled = (lam * (sp.Rational(1, 2) + u)) / (lam * (sp.Rational(1, 2) - u))
    checks["G_ratio_scale_invariant"] = sp.simplify(ratio_scaled - ratio) == 0

    # ---------------------------------------------------------------------------------------------
    print("Front 3 - determinant-line algebra and level ordering of a model operator (exact symbolic)")
    print("=" * 90)
    print("  (A) dim Sym^2(C^2) = 3 (= dim sl_2);  dim wedge^2(C^2) = 1 (determinant line)")
    print("  (B) action of g (x) g on wedge^2 = det g, computed; = 1 on SL(2,C) and su(2) slices => trivial line")
    print("  (C) gl_2 generator diag(y1,y2) acts on wedge^2 by trace y1+y2 (needs structure group U(2), [H-Weak])")
    print("  (D) model operator (C2-J3^2)/C2 + u J3 = diag(1, 1/2+u, 1/2-u)   (v=0, CP-even)")
    print("  (E) u=0 reduces to diag(1, 1/2, 1/2)  (algebraic value of (C2-J3^2)/C2 at C2=2; BI reading not used)")
    print("  (F) exit deficits d(e_0)=0, d(e_+)=1/2-u, d(e_-)=1/2+u: ordered, outer levels > 0, iff 0 < u < 1/2")
    print("  (G) level ratio (1/2+u)/(1/2-u) scale-invariant; absolute mass NOT fixed here")
    print("-" * 90)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 90)
    print("METHOD: under [H-Weak] L_Y supplies the line from which weights are built (exponents k, m open); the model")
    print("        operator carries the levels (read as E_Pi^2 only under [H-Res]); the mass comes after (Yukawa")
    print("        operator Y_Pi open: norm + sign(u) + mixing).")
    print("NOT TESTED: [H-Res], the level-to-generation map, [H-Spin], [H-Weak], the uniqueness enumeration of")
    print("            functors, membership of a generator in a structure group.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
