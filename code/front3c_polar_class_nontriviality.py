"""Front 3c step 3: polar-class non-triviality audit of the generator of the chiral polar factor.

Bias-independent, exact symbolic verification (no sampling). The canonical-vs-class audit
(front3c_polar_class_audit.py) established that the chiral polar factor U_Pi of the projected Yukawa
Y^gen_Pi = U_Pi H_Pi^{1/2} is not a canonical observable but a rephasing CLASS

    [U_Pi]  in  U(3) / (U(1)^3_R x U(1)^3_L),     U_Pi ~ V_R U_Pi V_L^{-1}  (V_L,V_R diagonal).

The observable is therefore the class [U_Pi] and not a matrix. This audit asks:

    Does a generator A of U_Pi(gamma) = exp(gamma A) give a NON-TRIVIAL polar class  [U_Pi(gamma)] != [I],
    or only a representative diagonal-rephasing-equivalent to the identity?

Geometric criterion (proved below). The diagonal-rephasing orbit through I has tangent space exactly the
DIAGONAL anti-hermitian matrices (i d, d real). A unitary generator A := dU_Pi/dgamma|_0 is anti-hermitian;
its class moves off [I] at first order iff A has a component TRANSVERSE to that orbit tangent, i.e. iff its
OFF-DIAGONAL part is non-zero. The rephasing-invariant detectors are:
  * the first variation of the moduli  d|(U_Pi)_{ij}|  (i != j), which is |A_{ij}|;
  * the Jarlskog-type invariant  J_CP = Im( U_11 U_22 conj(U_12) conj(U_21) ), whose leading non-vanishing
    order in gamma is gamma^3 and is proportional to Im(A_12 A_23 A_31) -- the genuine 3-generation CP datum.

Two precautions.
  (i)  AOG/CHO treat the ORIENTATION caveat [H-orient] and the spin-Galois sqrt(5) factor; they do NOT settle
       the non-triviality of [U_Pi]. Under [H-Spin], [H-WS], [M], [B] and [C] of AOG, the chiral lift induces
       no action on Q(sqrt 5), so no spin-Galois factor orthogonal to zeta_q survives -- but that concerns orientation, not the polar class.
  (ii) u != 0 (the diagonal generation split) is NOT mixing. u lives in the DIAGONAL channel (a level split that
       commutes with the rephasings) and contributes nothing to the off-diagonal moduli or to J_CP. Mixing
       needs a generator with a non-zero TRANSVERSE part A_off. Whether the step model supplies one is the
       question of front3d_transverse_route.py and front3e_reality_structure.py: for the sl_2 lift under the
       internal antilinear parity J_Pi of Q14 Section 6 the projection A_Pi^odd vanishes, so its off-diagonal part is
       zero; A_Pi^odd is not the polar generator (which is defined in PYO and is, at the identity, the J_Pi-even
       antiherm(L(M)), non-zero), it does not determine the polar factor, and nothing about the polar class or
       mixing follows (A_Pi^odd is the projection defined in PYO, not the anomaly density of Q14).

Results (all exact symbolic). J_CP denotes the Jarlskog-type phase; J_Pi is reserved for the antilinear parity.
  (A) Orbit tangent at I is diagonal anti-hermitian: d/ds [V_R(s) V_L(s)^{-1}]|_0 = i(D_R - D_L), diagonal,
      anti-hermitian; off-diagonal part identically zero; and every diagonal anti-hermitian is reached.
  (B) Transversality detector: for U(gamma) = exp(gamma A), A anti-hermitian, the second-order coefficient of
      |(U)_{ij}|^2 (i != j) is |A_{ij}|^2; so d|U_{ij}| != 0 iff A_{ij} != 0, i.e. iff A is transverse.
  (C) Jarlskog order: J_CP vanishes at orders gamma^0, gamma^1, gamma^2; its leading gamma^3 coefficient is
      proportional to Im(A_12 A_23 A_31), independent of the diagonal phases d_k (rephasing-invariant).
  (D) Diagonal generator: if A_off = 0, the exponential series of A has zero off-diagonal entries and J_CP = 0
      through the computed order, so [U_Pi] = [I]. (The series is computed from A, not typed in.)
  (E) Transverse generator: A_off != 0 gives off-diagonal moduli != 0 ([U_Pi] != [I], CP-conserving mixing if A_off
      is real); a complex A_off with Im(A_12 A_23 A_31) != 0 gives J_CP != 0 at order gamma^3.
  (F) Diagonal level split: the model operator diag(1, 1/2+u, 1/2-u) commutes with every diagonal rephasing, and
      does not commute with a generic non-diagonal unitary (negative control): the diagonal split is not a class
      datum of U_Pi.

Conclusion (printed). The non-triviality of the polar class is controlled EXACTLY by the off-diagonal (transverse)
part of the generator of U_Pi. This script is generic in the anti-hermitian generator A and does not say which A the
step model supplies; that is the object of front3d_transverse_route.py and front3e_reality_structure.py. No mass and
no mixing value is produced. No figures. English.
"""

import sympy as sp


def anti_hermitian(d, off):
    """Build a 3x3 anti-hermitian matrix: diagonal i*d_k (d real); off[(i,j)] = A_{ij}, A_{ji} = -conj(A_{ij})."""
    A = sp.zeros(3, 3)
    for k in range(3):
        A[k, k] = sp.I * d[k]
    for (i, j), val in off.items():
        A[i, j] = val
        A[j, i] = -sp.conjugate(val)
    return A


def jarl(M):
    return sp.im(M[0, 0] * M[1, 1] * sp.conjugate(M[0, 1]) * sp.conjugate(M[1, 0]))


def expm_series(A, g, order):
    """Truncated matrix exponential exp(g A) up to g^order (exact, symbolic)."""
    U = sp.zeros(3, 3)
    term = sp.eye(3)
    for k in range(order + 1):
        U += term * (g**k) / sp.factorial(k)
        term = sp.simplify(term * A)
    return U


def series_coeff(expr, g, n):
    return sp.expand_complex(sp.expand(expr)).coeff(g, n)


def main():
    checks = {}
    g = sp.symbols("gamma", real=True)

    # real diagonal phases and generic complex off-diagonal couplings (anti-hermitian generator)
    d = sp.symbols("d0:3", real=True)
    a12, b12, a13, b13, a23, b23 = sp.symbols("a12 b12 a13 b13 a23 b23", real=True)
    off_full = {(0, 1): a12 + sp.I * b12, (0, 2): a13 + sp.I * b13, (1, 2): a23 + sp.I * b23}
    A = anti_hermitian(list(d), off_full)

    # ---- (A) orbit tangent at I = diagonal anti-hermitian ---------------------------------------
    s = sp.symbols("s", real=True)
    dl = sp.symbols("dl0:3", real=True)
    dr = sp.symbols("dr0:3", real=True)
    VL = sp.diag(*[sp.exp(sp.I * s * dl[k]) for k in range(3)])
    VR = sp.diag(*[sp.exp(sp.I * s * dr[k]) for k in range(3)])
    curve = VR * VL.inv()
    tangent = sp.diff(curve, s).subs(s, 0)
    checks["A_tangent_offdiag_zero"] = all(sp.simplify(tangent[i, j]) == 0
                                           for i in range(3) for j in range(3) if i != j)
    checks["A_tangent_diag_antiherm"] = all(
        sp.simplify(tangent[k, k] - sp.I * (dr[k] - dl[k])) == 0 for k in range(3))
    # every diagonal anti-hermitian is reached (set dl = 0): i*dr arbitrary
    checks["A_tangent_surjective_diag"] = all(
        sp.simplify(tangent.subs({dl[k]: 0 for k in range(3)})[k, k] - sp.I * dr[k]) == 0
        for k in range(3))

    # ---- (B) modulus first variation detects the transverse (off-diagonal) part -----------------
    U2 = expm_series(A, g, 2)
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        mod2 = sp.expand_complex(sp.expand(U2[i, j] * sp.conjugate(U2[i, j])))
        coeff_g2 = mod2.coeff(g, 2)
        target = sp.Abs(off_full[(i, j)])**2          # |A_{ij}|^2
        checks[f"B_modvar_{i}{j}"] = sp.simplify(coeff_g2 - target) == 0

    # ---- (C) Jarlskog order: 0 up to g^2, leading g^3 ~ Im(A_12 A_23 A_31) -----------------------
    U3 = expm_series(A, g, 3)
    J = sp.expand_complex(sp.expand(jarl(U3)))
    checks["C_J_order0"] = series_coeff(J, g, 0) == 0
    checks["C_J_order1"] = series_coeff(J, g, 1) == 0
    checks["C_J_order2"] = series_coeff(J, g, 2) == 0
    J3 = sp.simplify(series_coeff(J, g, 3))
    # A_31 = -conj(A_13); the genuine 3-gen CP datum is Im(A_12 A_23 A_31)
    A31 = -sp.conjugate(off_full[(0, 2)])
    cp_datum = sp.im(off_full[(0, 1)] * off_full[(1, 2)] * A31)
    ratio = sp.simplify(J3 / cp_datum)
    # the proportionality constant is computed above and tested against its value (-1), not only for independence
    checks["C_J3_equals_minus_cpdatum"] = sp.simplify(J3 + cp_datum) == 0
    checks["C_J3_indep_diag"] = all(sp.simplify(sp.diff(J3, d[k])) == 0 for k in range(3))

    # ---- (D) diagonal generator: A_off = 0 => series has no off-diagonal entries, J_CP = 0 ----------------------
    sub_zero_off = {a12: 0, b12: 0, a13: 0, b13: 0, a23: 0, b23: 0}
    A_diag_only = A.subs(sub_zero_off)
    U_d = expm_series(A_diag_only, g, 3)
    checks["D_diag_generator_offdiag_series_zero"] = all(sp.simplify(sp.expand_complex(U_d[i, j])) == 0
                                                         for i in range(3) for j in range(3) if i != j)
    checks["D_diag_generator_JCP_zero"] = sp.simplify(sp.expand_complex(sp.expand(jarl(U_d)))) == 0

    # ---- (E) transverse generator vs CP-conserving real channel ------------------------------------------------
    # real off-diagonal channel (A_off real): [U_Pi] != [I] but J_CP = 0 (CP conserving)
    A_real_off = anti_hermitian([0, 0, 0], {(0, 1): a12, (0, 2): a13, (1, 2): a23})
    U_re = expm_series(A_real_off, g, 2)
    mod_re = sp.expand_complex(sp.expand(U_re[0, 1] * sp.conjugate(U_re[0, 1])))
    checks["E_real_off_mixing"] = sp.simplify(mod_re.coeff(g, 2) - a12**2) == 0
    U_re3 = expm_series(A_real_off, g, 3)
    checks["E_real_off_JCP_zero"] = sp.simplify(series_coeff(jarl(U_re3), g, 3)) == 0
    # complex channel with Im(A_12 A_23 A_31) != 0 => J_CP != 0 at g^3
    sub_cp = {a12: 1, b12: 0, a23: 1, b23: 0, a13: 0, b13: 1, d[0]: 0, d[1]: 0, d[2]: 0}
    cp_val = sp.simplify(cp_datum.subs(sub_cp))
    checks["E_complex_off_JCP_nonzero"] = cp_val != 0 and sp.simplify(J3.subs(sub_cp)) != 0

    # ---- (F) the diagonal level split commutes with the rephasings --------------------------------------------
    # u enters the model operator diag(1,1/2+u,1/2-u) on C^3_gen (E_Pi^2|gen only under [H-Res]): a diagonal level
    # operator. It commutes with every diagonal rephasing; a generic non-diagonal unitary does not (negative control).
    u = sp.symbols("u", real=True)
    Hlev = sp.diag(1, sp.Rational(1, 2) + u, sp.Rational(1, 2) - u)
    Vd = sp.diag(*[sp.exp(sp.I * dl[k]) for k in range(3)])
    checks["F_level_split_commutes_with_rephasings"] = all(
        sp.simplify(e) == 0 for e in (Vd * Hlev - Hlev * Vd))
    th = sp.symbols("theta", real=True)
    Rot = sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])
    checks["F_level_split_not_invariant_under_generic_rotation"] = not all(
        sp.simplify(e) == 0 for e in (Rot * Hlev - Hlev * Rot).subs({th: sp.pi / 3, u: sp.Rational(1, 10)}))

    # ---------------------------------------------------------------------------------------------
    print("Front 3c step 3 - polar-class non-triviality of the generator of U_Pi (exact symbolic)")
    print("=" * 100)
    print("  Question: does a generator A of U_Pi(gamma) = exp(gamma A) give [U_Pi(gamma)] != [I] in")
    print("  U(3)/(U(1)^3_R x U(1)^3_L)?")
    print("  Criterion: [U_Pi] moves off [I] at first order  <=>  A has a non-zero OFF-DIAGONAL (transverse) part")
    print("             (the orbit tangent at I is diagonal anti-hermitian).")
    print("  (A) orbit tangent at I = diagonal anti-hermitian (off-diagonal part identically zero)")
    print("  (B) modulus first variation: coeff of g^2 in |U_{ij}|^2 is |A_{ij}|^2  => detects A_off")
    print("  (C) J_CP: 0 up to g^2; leading g^3 term ~ Im(A_12 A_23 A_31), independent of diagonal phases")
    print("  (D) diagonal generator (A_off = 0): series stays diagonal, J_CP = 0 => [U_Pi] = [I]")
    print("  (E) transverse generator: [U_Pi] != [I]; Im(A_12 A_23 A_31) != 0 => J_CP != 0")
    print("  (F) the diagonal level split commutes with the rephasings (negative control: a rotation does not)")
    print("-" * 100)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 100)
    print(f"  checks run: {len(checks)}")
    print("RESULT: polar-class non-triviality is controlled EXACTLY by the off-diagonal (transverse) part of the")
    print("        generator of U_Pi. The script is generic in the anti-hermitian generator A: which A the step model")
    print("        supplies is the question of front3d_transverse_route.py and front3e_reality_structure.py (for the")
    print("        sl_2 lift under the internal antilinear parity J_Pi of Q14 Sec. 6, the projection A_Pi^odd")
    print("        vanishes; it is not the polar generator and determines no polar class).")
    print("        AOG/CHO treat [H-orient], not this.")
    print("        No mass and no mixing value is produced.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
