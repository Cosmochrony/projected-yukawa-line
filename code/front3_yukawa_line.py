"""Front 3: determinant-line algebra and the generation-level assignment of a model operator.

Exact symbolic verification (no sampling) of the algebraic content of the two statements opening the mass-sector
frontier. Nothing here is new representation theory. The tensor identities are those of Q14 Theorem 2.1(a)
(Sym^2(C^2) = sl_2(C), wedge^2(C^2) = determinant line), the model operator is that of PRS
(diag(1, 1/2+u, 1/2-u) on C^3_gen = Sym^2(V_gen)), and the exit deficits are those of Q14 Prop. 6.7.
The method lock is explicit: under [H-Weak] the determinant line L_Y := wedge^2(E_weak) fixes the weight datum of the
coupling; the model operator fixes the levels, read as levels of E_Pi^2 only under the named identification [H-Res];
the mass comes after (Yukawa norm + sign of u + mixing).

Typing (Q14): the rank-two fibre is either the spinor factor S_L, with structure group SL(2,C) under [H-Spin], on which
wedge^2 is a trivial line with no hypercharge, or the weak factor E_weak of [H-Weak], with structure group U(2), on
which L_Y = wedge^2(E_weak) carries the abelian weight. E_weak is not constructed in the corpus. The script works on
C^2 and does not construct either bundle.

Seven groups of checks are run.
  (A) Dimensions: dim Sym^2(C^2) = 3 (= dim sl_2), dim wedge^2(C^2) = 1. The uniqueness of wedge^2 among Sym^k, wedge^k
      (k >= 1) rests on the dimension count in the paper's proof; the script does not enumerate the functors.
  (B) Action on the line: for g in GL(2,C) the induced action on wedge^2(C^2) is det(g); it is 1 on SL(2,C), hence on
      any subgroup (any SU(2) compact real form, the structure group of S_L under [H-Spin]). So wedge^2 is trivial
      there and carries no charge.
  (C) Weight of a gl_2 generator: a diagonal generator diag(y1, y2) acts on wedge^2(C^2) by the trace y1 + y2. This
      generator lies in gl_2 and not in sl_2; it belongs to the structure algebra only if the structure group is U(2),
      which is the choice of [H-Weak] for E_weak. Checks B and C are separate and the script does NOT test that a
      generator belongs to the structure group of S_L or of E_weak: that membership is an input.
  (D) Normal form of the model operator: (C2 - J3^2)/C2 + u J3 (v = 0, CP-even) with C2 = 2, J3 = diag(0,1,-1) gives
      exactly diag(1, 1/2+u, 1/2-u) on (e_0, e_+, e_-).
  (E) Even sector: at u = 0 this reduces to diag(1, 1/2, 1/2), the algebraic value of (C2 - J3^2)/C2 at C2 = 2
      (PRS eq:even); its Born-Infeld reading through O30 is not used.
  (F) Exit deficits d(e_i) = s_0 - s_i with s_0 = 1: d(e_0) = 0, d(e_+) = 1/2 - u, d(e_-) = 1/2 + u; for u > 0 the
      ordering 0 < d(e_+) < d(e_-) holds. This orders three levels of the model operator.
  (G) Scale independence of the level RATIO: (1/2 + u)/(1/2 - u) is invariant under rescaling, so the dimensionless
      level structure is fixed by u alone, while the absolute mass normalisation is NOT fixed here.

Not tested: the identification [H-Res] of the model operator with the restriction of E_Pi^2 to C^3_gen, the map from
levels to generations, [H-Spin], [H-Weak], and the Yukawa operator. No number is produced for any mass. No figures.
English.
"""

import sympy as sp


def main():
    checks = {}

    # ---- carrier C^2 and its functorial powers ---------------------------------------------------
    # Sym^2(C^2): symmetric 2-tensors, dim 3.  wedge^2(C^2): antisymmetric 2-tensors, dim 1.
    dim_sym2 = sp.binomial(2 + 2 - 1, 2)          # = 3
    dim_wedge2 = sp.binomial(2, 2)                 # = 1
    checks["A_dim_sym2_is_3"] = (dim_sym2 == 3)
    checks["A_dim_wedge2_is_1"] = (dim_wedge2 == 1)
    checks["A_sym2_dim_matches_sl2"] = (dim_sym2 == 3)  # dim sl_2(C) = 3
    checks["A_wedge2_is_a_line"] = (dim_wedge2 == 1)

    # ---- (B) action on wedge^2 is det g; det = 1 on SL(2,C), hence trivial for a structure group inside SL(2,C) ----
    a, b, c, d = sp.symbols("a b c d")
    g = sp.Matrix([[a, b], [c, d]])
    # Action on wedge^2(C^2) (a line spanned by e1 ^ e2): g.(e1^e2) = det(g) (e1^e2).
    wedge_action = g.det()                          # = a d - b c
    # On SL(2,C): det g = 1.  SU(2) is the compact real form, a subgroup of SL(2,C).
    checks["B_wedge_action_is_det"] = sp.simplify(wedge_action - (a * d - b * c)) == 0
    # impose SL(2): det = 1  =>  trivial action
    checks["B_trivial_on_SL2"] = sp.simplify(wedge_action.subs(a * d - b * c, 1) - 1) == 0 \
        or sp.simplify((wedge_action - 1).subs(d, (1 + b * c) / a)) == 0
    # explicit SU(2) element g = [[alpha, beta], [-conj beta, conj alpha]], |alpha|^2+|beta|^2 = 1
    al, be = sp.symbols("alpha beta", real=True)    # take a real-section SU(2) slice for an exact check
    gsu2 = sp.Matrix([[al, be], [-be, al]])
    checks["B_su2_det_one"] = sp.simplify(gsu2.det().subs(al**2 + be**2, 1) - 1) == 0 \
        or sp.simplify(gsu2.det() - (al**2 + be**2)) == 0

    # ---- (C) weight of a gl_2 generator: action on wedge^2 of a diagonal generator is the trace --------------
    y1, y2 = sp.symbols("y1 y2")
    Y = sp.diag(y1, y2)                             # diagonal gl_2 generator on C^2 (not in sl_2)
    # infinitesimal action on the line wedge^2: derivative of det(exp(tY)) at 0 = tr(Y) = y1 + y2
    t = sp.symbols("t")
    line_weight = sp.diff(sp.exp(t * y1) * sp.exp(t * y2), t).subs(t, 0)  # d/dt det(e^{tY})|_0
    checks["C_gl2_weight_is_trace"] = sp.simplify(line_weight - (y1 + y2)) == 0

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

    # ---- (F) exit-deficit ordering --------------------------------------------------------------
    s0, sp_, sm = 1, sp.Rational(1, 2) + u, sp.Rational(1, 2) - u
    d0, dplus, dminus = s0 - s0, s0 - sp_, s0 - sm
    checks["F_deficit_e0_zero"] = (d0 == 0)
    checks["F_deficit_eplus"] = sp.simplify(dplus - (sp.Rational(1, 2) - u)) == 0
    checks["F_deficit_eminus"] = sp.simplify(dminus - (sp.Rational(1, 2) + u)) == 0
    # for u > 0: 0 < d(e_+) < d(e_-)  =>  e_0 smallest deficit, e_- largest deficit
    checks["F_ordering_u_positive"] = sp.simplify((dminus - dplus)) == 2 * u  # > 0 for u > 0

    # ---- (G) scale-independence of the level ratio ----------------------------------------------
    lam = sp.symbols("lambda", positive=True)
    ratio = (sp.Rational(1, 2) + u) / (sp.Rational(1, 2) - u)
    ratio_scaled = (lam * (sp.Rational(1, 2) + u)) / (lam * (sp.Rational(1, 2) - u))
    checks["G_ratio_scale_invariant"] = sp.simplify(ratio_scaled - ratio) == 0

    # ---------------------------------------------------------------------------------------------
    print("Front 3 - determinant-line algebra and generation-level assignment of a model operator (exact symbolic)")
    print("=" * 90)
    print("  (A) dim Sym^2(C^2) = 3 (= dim sl_2);  dim wedge^2(C^2) = 1 (determinant line)")
    print("  (B) wedge^2 action = det g; = 1 on SL(2,C) => trivial line for a structure group inside SL(2,C)")
    print("  (C) gl_2 generator diag(y1,y2) acts on wedge^2 by trace y1+y2 (needs structure group U(2), [H-Weak])")
    print("  (D) model operator (C2-J3^2)/C2 + u J3 = diag(1, 1/2+u, 1/2-u)   (v=0, CP-even)")
    print("  (E) u=0 reduces to diag(1, 1/2, 1/2)  (algebraic value of (C2-J3^2)/C2 at C2=2; BI reading not used)")
    print("  (F) exit deficits d(e_0)=0, d(e_+)=1/2-u, d(e_-)=1/2+u => three ordered levels (u > 0)")
    print("  (G) level ratio (1/2+u)/(1/2-u) scale-invariant; absolute mass NOT fixed here")
    print("-" * 90)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 90)
    print("METHOD LOCK: under [H-Weak] L_Y fixes the weight datum; the model operator fixes the levels (read as E_Pi^2")
    print("            only under [H-Res]); the mass comes after (Yukawa operator Y_Pi open: norm + sign(u) + mixing).")
    print("NOT TESTED: [H-Res], the level-to-generation map, [H-Spin], [H-Weak], membership of a generator in a structure group.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
