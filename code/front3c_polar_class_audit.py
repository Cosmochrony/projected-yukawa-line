"""Front 3c: canonical-vs-class audit of the chiral polar factor U_Pi.

Exact symbolic verification (no sampling). The projected Yukawa block on the generation carriers is written
    Y^gen_Pi = U_Pi H_Pi^{1/2},   H_Pi := Y^gen_Pi^dag Y^gen_Pi  on the left generation carrier,
with H_Pi|gen = lambda_Y^2 diag(1, 1/2+u, 1/2-u) taken as a model operator on C^3_gen (its reading as a restriction of
E_Pi^2 is not supplied; PYO Beau2026pyo / PRS Beau2026prs / PYL Beau2026pyl). The polar factor
U_Pi : G_L -> G_R between the left and right generation carriers is not determined by H_Pi. The audit asks: which data
of U_Pi survive the admissible basis changes,

    U_Pi  ~  V_R U_Pi V_L^{-1}        (V_L, V_R: admissible rephasings of the left and right generation carriers) ?

Admissible rephasings (as in PYO, Proposition on the polar class). The generation block lives on two carriers, each a
copy of C^3_gen with its own Cartan generator J_3 = diag(0, 1, -1) on the weight basis (e_0, e_+, e_-). The admissible
rephasings of a carrier are the unitaries that preserve its J_3 grading, that is, the unitaries commuting with J_3. The
weights 0, +1, -1 are distinct, so these are the diagonal unitaries, U(1)^3_L and U(1)^3_R. The groups are defined by
the carriers and their gradings, not by Y^gen_Pi, and not by the eigenbasis of Y^gen_Pi Y^gen_Pi^dag. H_Pi is
positive and commutes with J_3; for 0 < u < 1/2 its levels are distinct (under [H-Res] and [H-Sq]), so on the left
carrier the eigenbasis of H_Pi is the weight basis.

Results (all exact symbolic).
  (A) Unrestricted right basis: with V_R unrestricted in U(3), V_R = U_Pi removes U_Pi entirely
      (U_Pi^{-1} Y^gen_Pi = H_Pi^{1/2}, diagonal positive), so only Spec(H_Pi) would survive.
  (B) Admissible groups: the commutant of J_3 is the set of diagonal matrices ((W J_3 - J_3 W)_{ij} =
      W_{ij}(J_j - J_i) with J_j - J_i != 0 for i != j). H_Pi is invariant under the left rephasings (a diagonal
      matrix commutes with a diagonal one), and not under a generic rotation (negative control). The removing choice
      V_R = U_Pi commutes with J_3 only if U_Pi is diagonal, so it is not admissible for a generic U_Pi (the check
      fails for a non-diagonal U_Pi and passes for a diagonal one). Y^gen_Pi Y^gen_Pi^dag = U_Pi H_Pi U_Pi^dag has
      the same spectrum as H_Pi; this is a fact about square matrices and is not used to define the groups.
  (C) Class invariants under diagonal rephasing U_Pi -> V_R U_Pi V_L^{-1}: the moduli |(U_Pi)_{ij}| and the quartet
      phase J_CP := Im( U_11 U_22 conj(U_12) conj(U_21) ) are invariant.
  (D) Parameter count: U(3) has 9 real parameters; the rephasings remove 2*3 - 1 = 5; so 4 survive, three mixing
      moduli (angles) and one CP phase. U_Pi is a class, not a canonical matrix.
  (E) A real orthogonal U_Pi has J_CP = 0. H_Pi is blind to U_Pi: Y^gen_Pi^dag Y^gen_Pi is the same for every
      unitary U_Pi, so neither H_Pi nor, under [H-Res] and [H-Sq], E_Pi^2|gen constrains U_Pi.
  (F) Illustration with a generator chosen by hand (not the generator of Q14's step model, for which
      front3e_reality_structure.py shows A_Pi = 0): U_Pi(gamma) = exp(gamma A) with A the real antisymmetric matrix
      E_01 - E_10. It has |U_01| = |sin gamma|, non-zero for generic gamma, while every element of the class of I has
      |U_ij| = delta_ij; so a generator with a transverse part moves the class. Real gamma keeps J_CP = 0.

Conclusion (printed). The class [U_Pi] of the polar factor, not a matrix, is what the construction carries: three
moduli and one CP phase. Which generators of the step model move the class off [I] is the transverse-route question
(front3c_polar_class_nontriviality.py, front3d_transverse_route.py, front3e_reality_structure.py). No mass and no
mixing value is produced. No figures. English.
"""

import sympy as sp


def _zero(M):
    """Robust matrix-zero test: complex- and trig-aware simplification entrywise."""
    return all(sp.simplify(sp.trigsimp(sp.expand_complex(e))) == 0 for e in M)


def dag(M):
    return M.conjugate().T


def main():
    checks = {}

    u = sp.symbols("u", real=True)
    I3 = sp.eye(3)

    # the class statements are proved for GENERIC positive singular values (a, b, c); the level identification enters
    # only through the distinctness of the three levels on the left carrier.
    a, b, c = sp.symbols("a b c", positive=True)
    Hhalf = sp.diag(a, b, c)
    H = Hhalf * Hhalf
    levels = [sp.Integer(1), sp.Rational(1, 2) + u, sp.Rational(1, 2) - u]
    asm = sp.Q.positive(sp.Rational(1, 2) + u) & sp.Q.positive(sp.Rational(1, 2) - u) & sp.Q.positive(u)

    # a generic unitary U_Pi (a real generation rotation and a phase); Y^gen_Pi = U_Pi H_Pi^{1/2}
    th, ph = sp.symbols("theta phi", real=True)
    U_rot = sp.Matrix([[sp.cos(th), -sp.sin(th), 0],
                       [sp.sin(th), sp.cos(th), 0],
                       [0, 0, 1]])
    U_ph = sp.diag(1, 1, sp.exp(sp.I * ph))
    U_pi = U_ph * U_rot
    Y = U_pi * Hhalf

    # ---- (A) unrestricted right basis removes U_Pi ----------------------------------------------------------
    checks["A_unitary_Upi"] = _zero(dag(U_pi) * U_pi - I3)
    checks["A_free_right_removes_Upi"] = _zero(U_pi.inv() * Y - Hhalf)
    checks["A_invariant_is_spectrum"] = _zero(dag(Y) * Y - H)

    # ---- (B) admissible groups: the commutant of each carrier's J_3 ------------------------------------------
    J3 = sp.diag(0, 1, -1)
    w = sp.symbols("w0:9")
    W = sp.Matrix(3, 3, lambda i, j: w[3 * i + j])
    comm = W * J3 - J3 * W
    j3 = [J3[k, k] for k in range(3)]
    checks["B_commutant_entries_are_W_ij_times_weight_difference"] = all(
        sp.simplify(comm[i, j] - W[i, j] * (j3[j] - j3[i])) == 0 for i in range(3) for j in range(3))
    checks["B_weights_distinct_so_commutant_is_diagonal"] = all(
        j3[j] - j3[i] != 0 for i in range(3) for j in range(3) if i != j)
    # H_Pi (diagonal in the weight basis, distinct levels) is invariant under the diagonal rephasings
    pL = sp.symbols("pL0:3", real=True)
    pR = sp.symbols("pR0:3", real=True)
    VL = sp.diag(*[sp.exp(sp.I * pL[i]) for i in range(3)])
    VR = sp.diag(*[sp.exp(sp.I * pR[i]) for i in range(3)])
    checks["B_H_invariant_under_left_rephasing"] = _zero(VL * H * VL.inv() - H)
    # negative control: a generic rotation of the weight basis does not leave H invariant
    rot = U_rot.subs(th, sp.pi / 3)
    H_num = H.subs({a: 1, b: 2, c: 3})
    checks["B_H_not_invariant_under_generic_rotation"] = not _zero(rot * H_num * rot.T - H_num)
    checks["B_H_invariant_under_trivial_rotation"] = _zero(U_rot.subs(th, 0) * H_num * U_rot.subs(th, 0).T - H_num)
    # distinctness of the three levels for 0 < u < 1/2
    distinct = [sp.refine(sp.simplify(levels[i] - levels[j]), asm) for i, j in [(0, 1), (0, 2), (1, 2)]]
    interval = sp.Interval.open(0, sp.Rational(1, 2))
    checks["B_levels_distinct"] = all(sp.solveset(d, u, domain=interval) == sp.EmptySet for d in distinct)
    # negative control: at u = 1/2 two levels coincide-or-vanish, the test must see a zero of a difference at the edge
    checks["B_levels_distinct_fails_outside"] = any(
        sp.solveset(d, u, domain=sp.Interval(0, sp.Rational(1, 2))) != sp.EmptySet for d in distinct)
    # the removing choice V_R = U_Pi is admissible only when U_Pi commutes with J_3, i.e. is diagonal
    gen_U = U_pi.subs({th: sp.pi / 3, ph: sp.Rational(1, 2)})
    checks["B_removing_choice_not_admissible_for_nondiagonal_Upi"] = not _zero(gen_U * J3 - J3 * gen_U)
    diag_U = U_pi.subs(th, 0)
    checks["B_removing_choice_admissible_when_Upi_diagonal"] = _zero(diag_U * J3 - J3 * diag_U)
    # right level operator has the same spectrum as H (fact about square matrices; not used to define the groups)
    x = sp.symbols("x")
    H_right = Y * dag(Y)
    cpL = H.charpoly(x).as_expr()
    cpR = H_right.charpoly(x).as_expr()
    checks["B_right_level_operator_same_spectrum"] = sp.simplify(sp.expand(sp.expand_complex(cpL - cpR))) == 0

    # ---- (C) class invariants under diagonal rephasing ---------------------------------------------------------
    Uent = sp.symbols("U0:9")
    Uabs = sp.Matrix(3, 3, lambda i, j: Uent[3 * i + j])
    Urep = VR * Uabs * VL.inv()
    checks["C_moduli_invariant"] = all(
        sp.simplify(sp.Abs(Urep[i, j]) - sp.Abs(Uabs[i, j])) == 0 for i in range(3) for j in range(3))

    def jarl(M):
        return M[0, 0] * M[1, 1] * sp.conjugate(M[0, 1]) * sp.conjugate(M[1, 0])

    checks["C_JCP_invariant"] = sp.simplify(sp.expand_complex(sp.im(jarl(Urep)) - sp.im(jarl(Uabs)))) == 0

    # ---- (D) parameter count: 9 - (2*3 - 1) = 4 = 3 angles + 1 phase -----------------------------------------
    n = 3
    physical = n * n - (2 * n - 1)
    angles = n * (n - 1) // 2
    phases = (n - 1) * (n - 2) // 2
    checks["D_param_count_identity"] = (physical == 4 and angles == 3 and phases == 1 and angles + phases == physical)
    # testable counterpart: the real rank of the rephasing action (dl, dr) -> i (D_r U - U D_l) at a generic unitary
    # is 2n - 1 = 5 (the common phase acts trivially), so the orbit has dimension 5 inside the 9-dimensional U(3).
    dls = sp.symbols("dl0:3", real=True)
    drs = sp.symbols("dr0:3", real=True)
    def rot(i, j, ang):
        R = sp.eye(3)
        R[i, i] = R[j, j] = sp.cos(ang)
        R[i, j], R[j, i] = -sp.sin(ang), sp.sin(ang)
        return R

    # a unitary with all nine entries non-zero and non-aligned phases (not block-structured)
    Ugen = (sp.diag(1, sp.exp(sp.I * sp.pi / 4), sp.exp(sp.I * sp.pi / 3)) * rot(0, 1, sp.pi / 3) * rot(1, 2, sp.pi / 4)
            * sp.diag(1, 1, sp.I) * rot(0, 2, sp.pi / 6))
    checks["D_generic_unitary_all_entries_nonzero"] = all(sp.simplify(e) != 0 for e in Ugen)
    tang = sp.I * (sp.diag(*drs) * Ugen - Ugen * sp.diag(*dls))
    rows = []
    for i in range(3):
        for j in range(3):
            e = sp.expand_complex(tang[i, j])
            rows.append([sp.re(e).diff(v) for v in (*dls, *drs)])
            rows.append([sp.im(e).diff(v) for v in (*dls, *drs)])
    checks["D_rephasing_orbit_dimension_5"] = sp.Matrix(rows).rank() == 2 * n - 1

    # ---- (E) real U_Pi has J_CP = 0; H_Pi is blind to U_Pi -----------------------------------------------------
    checks["E_real_unitary_JCP_zero"] = sp.simplify(sp.expand_complex(sp.im(jarl(U_pi.subs(ph, 0))))) == 0
    other = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])      # another unitary polar factor
    checks["E_Hpi_blind_to_Upi"] = _zero(dag(other * Hhalf) * (other * Hhalf) - H) and \
        _zero(dag(U_pi * Hhalf) * (U_pi * Hhalf) - H)

    # ---- (F) illustration: a generator chosen by hand moves the class ------------------------------------------
    g = sp.symbols("gamma", real=True)
    A = sp.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])          # real antisymmetric, hand-chosen
    U_gamma = sp.simplify(sp.exp(A * g))
    checks["F_gamma0_identity"] = _zero(U_gamma.subs(g, 0) - I3)
    checks["F_moduli_move_off_class_of_identity"] = sp.simplify(sp.Abs(U_gamma[0, 1]) ** 2 - sp.sin(g) ** 2) == 0 \
        and sp.simplify(sp.Abs(U_gamma[0, 1]).subs(g, sp.Rational(1, 3))) != 0
    checks["F_real_gamma_JCP_zero"] = sp.simplify(sp.expand_complex(sp.im(jarl(U_gamma)))) == 0

    # ---- report ------------------------------------------------------------------------------------------------
    print("Front 3c - canonical-vs-class audit of the chiral polar factor U_Pi (exact symbolic, no sampling)")
    print("=" * 100)
    print("  Y^gen_Pi = U_Pi H_Pi^{1/2};  question: which data of U_Pi survive U_Pi ~ V_R U_Pi V_L^{-1}?")
    print("  (A) with V_R unrestricted in U(3), V_R = U_Pi removes U_Pi: only Spec(H_Pi) would survive")
    print("  (B) admissible V_L, V_R = unitaries commuting with each carrier's J_3 = diagonal unitaries (distinct")
    print("      weights); V_R = U_Pi is admissible only if U_Pi is diagonal; no use of the eigenbasis of Y Y^dag")
    print("  (C) under diagonal rephasing the moduli |(U_Pi)_{ij}| and the quartet phase J_CP are INVARIANT")
    print("  (D) count: 9 - (2*3-1) = 4 = 3 mixing moduli + 1 CP phase => a CLASS, not a canonical matrix")
    print("  (E) real U_Pi => J_CP = 0; H_Pi = Y^dag Y is blind to U_Pi")
    print("  (F) a hand-chosen generator with a transverse part moves the class (|U_01| = |sin gamma|)")
    print("-" * 100)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 100)
    print(f"  checks run: {len(checks)}")
    print("RESULT: the construction carries the class [U_Pi] under the rephasings that preserve each carrier's J_3")
    print("        grading, not a matrix: three moduli |(U_Pi)_{ij}| and one CP phase J_CP. H_Pi, hence E_Pi^2 under")
    print("        [H-Res] and [H-Sq], does not fix U_Pi. Which generators of the step model move the class off [I]")
    print("        is the transverse-route question (front3c_polar_class_nontriviality.py, front3d, front3e).")
    print("        No mass and no mixing value is produced.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
