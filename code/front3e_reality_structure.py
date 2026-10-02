"""Front 3e: reality structure of the metaplectic step under the antilinear parity J_Pi of Q14.

Exact symbolic verification (no sampling). Object under audit: the sl_2(C) step generator M = p E + q F + r H with
complex coefficients (p, q, r), its derived Sym^2 lift L(M) on C^3_gen = Sym^2(V_gen), and the generator

    A_Pi := antiherm( J_Pi-odd( L(M) ) )                                   (definition, taken as given)

of the chiral unitary polar factor. This script tests what that definition yields; it does not test why the
generator is defined through the J_Pi-ODD part (see the open modelling question below).

Definitions (Q14 Section 6; Q14 Theorem 3.7 for the spinorial lift).
  * V = C^2 = <v_+, v_->, e_0 = sqrt(2) v_+ v_-, e_+ = v_+ v_+, e_- = v_- v_-, J_3 = diag(0, 1, -1).
  * J_Pi on V is the antilinear map z -> eps conj(z) with eps = [[0, 1], [-1, 0]] (J_Pi^2 = -1 on V). Lifted to
    C^3_gen = Sym^2(V) it is the ANTILINEAR map  J_Pi z = S conj(z)  with S = Sym^2(eps), which maps e_0 to -e_0 and
    exchanges e_+ <-> e_- (Q14 Sec. 6). S is computed below from eps, not typed in.
  * Conjugation of an operator by the antilinear J_Pi:  J_Pi X J_Pi^{-1} = S conj(X) S^{-1}. The J_Pi-odd part is
    (X - S conj(X) S^{-1}) / 2 and the J_Pi-even part is (X + S conj(X) S^{-1}) / 2. The map X -> S conj(X) S^{-1}
    is real-linear and not complex-linear in X.

Results (all exact symbolic, generic complex (p, q, r)).
  (1) Key identity: S conj(L(M)) S^{-1} = - L(M)^dagger for every M in sl_2(C).
  (2) Hence the J_Pi-odd part of L(M) is its hermitian part and the J_Pi-even part is its anti-hermitian part.
  (3) Hence A_Pi = antiherm(J_Pi-odd(L(M))) = 0 identically in (p, q, r) in C^3: the reality of p, q, r plays no
      role. The J_3 coefficient of the odd part is 2 Re r, its R_mix coefficient is 0, and the diagonal
      i Im(r) diag(0, 2, -2) is J_Pi-EVEN, so it is not in A_Pi either.
  (4) Negative control, the LINEAR involution X -> S X S^{-1} (the matrix S without the conjugation). It is not
      Q14's J_Pi. It gives a transverse part proportional to Im(p + q), not "Im p or Im q" (Im p = - Im q gives
      zero), and it coincides with the antilinear parity on real data. The two differ exactly on the imaginary part
      of X, which the antilinear parity exchanges between odd and even.
  (5) Scope of the vanishing: it is specific to A_Pi as defined. The J_Pi-odd anti-hermitian operators of u(3) form a
      6-dimensional real space that is Hilbert-Schmidt orthogonal to L(sl_2(C)); the internal block e_0 <-> e_+/-
      and the external block R_mix (e_+ <-> e_-) both live in it, so a J_Pi-odd anti-hermitian source of either block
      (A_Pi as defined) needs a generator outside the sl_2 image, in the spin-2 sector of End(Sym^2 V_gen). No
      source supplies such a generator. The external block R_mix is not reached by the image of sl_2 at all (Q14
      Remark 6.4); the internal block is reached through the J_Pi-even part (see the open modelling question).
  (6) Frame covariance: in a complex unitary frame W with the transported antilinear parity S' = W S W^T, the
      generator is W A_Pi W^dagger = 0 (true and vacuous).
  (7) A central imaginary scalar i c I_3 is J_Pi-odd and anti-hermitian; it is diagonal, hence class-trivial.

Open modelling question (not settled here). The J_Pi-EVEN anti-hermitian part of L(M), antiherm(L(M)), is not zero:
its internal entry is sqrt(2) (q - conj p) / 2, non-zero for real p != q. The script computes it and reports it;
so the image of sl_2 does reach the internal block e_0 <-> e_+/- through that part (the external block R_mix is not
reached). Whether the polar generator is rightly the J_Pi-odd part is a modelling choice of the companion note that
the sources used here do not justify; Q14 excludes the internal block only by hypothesis (Prop. 6.3 (i)-(ii)). No
statement that the internal block is excluded as a physical matter follows.

Not tested: the identification [H-Res], the map from levels to generations, [H-Spin], [H-Weak], the existence of
any generator in the spin-2 sector. No mass and no mixing value is produced. No figures. English.
"""

import sympy as sp

s2 = sp.sqrt(2)


def zero(X):
    """Exact zero test of a matrix with symbolic complex entries."""
    return all(sp.simplify(sp.expand_complex(e)) == 0 for e in X)


def herm(X):
    return (X + X.H) / 2


def antiherm(X):
    return (X - X.H) / 2


def sym2_of_matrix(g):
    """Group action Sym^2(g) of a 2x2 matrix g on the orthonormal basis (e_0, e_+, e_-) of Sym^2(C^2).

    A symmetric tensor is represented by a symmetric 2x2 matrix T (e_+ -> diag(1,0), e_- -> diag(0,1),
    e_0 -> offdiag(1,1)/sqrt2); g acts by T -> g T g^T.
    """
    basis = [sp.Matrix([[0, 1], [1, 0]]) / s2, sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    cols = []
    for T in basis:
        Tn = g * T * g.T
        cols.append([Tn[0, 1] * s2, Tn[0, 0], Tn[1, 1]])
    return sp.Matrix(cols).T


def sym2_lift(M):
    """Derived (Lie algebra) action d/dt Sym^2(1 + t M) at t = 0 of a 2x2 matrix M."""
    t = sp.Symbol("t")
    return sp.simplify(sym2_of_matrix(sp.eye(2) + t * M).diff(t).subs(t, 0))


def main():
    checks = {}
    E = sp.Matrix([[0, 1], [0, 0]])
    F = sp.Matrix([[0, 0], [1, 0]])
    H = sp.Matrix([[1, 0], [0, -1]])
    I3 = sp.eye(3)
    eps = sp.Matrix([[0, 1], [-1, 0]])
    S = sym2_of_matrix(eps)                      # Sym^2(eps), the matrix of the antilinear lift z -> S conj(z)
    Sinv = S.inv()

    def jconj(X):
        """J_Pi X J_Pi^{-1} for the antilinear J_Pi z = S conj(z)."""
        return S * X.conjugate() * Sinv

    def lin(X):
        """The LINEAR involution X -> S X S^{-1} (negative control, not Q14's J_Pi)."""
        return S * X * Sinv

    def odd(X):
        return (X - jconj(X)) / 2

    pR, pI, qR, qI, rR, rI = sp.symbols("pR pI qR qI rR rI", real=True)
    p, q, r = pR + sp.I * pI, qR + sp.I * qI, rR + sp.I * rI
    M = p * E + q * F + r * H
    L = sym2_lift(M)

    # ---- (0) definitions, computed -------------------------------------------------------------------------
    checks["0a_S_from_eps_is_e0_flip_and_exchange"] = S == sp.Matrix([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    checks["0b_S_real_unitary_involution"] = zero(S * S - I3) and zero(S - S.conjugate()) and zero(S * S.H - I3)
    checks["0c_J_squared_minus_one_on_V"] = zero(eps * eps.conjugate() + sp.eye(2))   # (z -> eps conj z)^2 = -1
    checks["0d_lift_explicit_form"] = L == sp.Matrix([[0, s2 * q, s2 * p], [s2 * p, 2 * r, 0], [s2 * q, 0, -2 * r]])
    checks["0e_lift_is_star_rep"] = zero(L.H - sym2_lift(M.H))
    checks["0f_H_lifts_to_2J3"] = sym2_lift(H) == sp.diag(0, 2, -2)

    # ---- (1) key identity ------------------------------------------------------------------------------------
    checks["1_S_conjL_Sinv_equals_minus_Ldagger"] = zero(jconj(L) + L.H)

    # ---- (2) odd part = hermitian part, even part = anti-hermitian part -------------------------------------
    checks["2a_odd_part_is_hermitian_part"] = zero(odd(L) - herm(L))
    checks["2b_even_part_is_antihermitian_part"] = zero((L + jconj(L)) / 2 - antiherm(L))
    expected_odd = sp.Matrix([
        [0, s2 * (q + sp.conjugate(p)) / 2, s2 * (p + sp.conjugate(q)) / 2],
        [s2 * (p + sp.conjugate(q)) / 2, 2 * rR, 0],
        [s2 * (q + sp.conjugate(p)) / 2, 0, -2 * rR]])
    checks["2c_odd_part_explicit"] = zero(odd(L) - expected_odd)

    # ---- (3) A_Pi vanishes identically ----------------------------------------------------------------------
    A_pi = antiherm(odd(L))
    checks["3a_A_Pi_zero_all_complex_pqr"] = zero(A_pi)
    checks["3b_A_Pi_zero_real_step"] = zero(A_pi.subs({pI: 0, qI: 0, rI: 0}))
    J3 = sp.diag(0, 1, -1)
    Rmix = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])        # Q14 Sec. 6: antisymmetric on (e_+, e_-)

    def hs(X, Y):
        return (X.H * Y).trace() / (Y.H * Y).trace()

    checks["3c_J3_coefficient_of_odd_part_is_2Re_r"] = sp.simplify(sp.expand_complex(hs(odd(L), J3) - 2 * rR)) == 0
    checks["3d_Rmix_coefficient_of_odd_part_zero"] = sp.simplify(sp.expand_complex(hs(odd(L), Rmix))) == 0
    cartan = sp.I * rI * sp.diag(0, 2, -2)
    checks["3e_Im_r_diagonal_is_J_even"] = zero(jconj(cartan) - cartan)
    even_ah = antiherm(L)
    checks["3f_even_antiherm_internal_entry"] = sp.simplify(
        sp.expand_complex(even_ah[0, 1] - s2 * (q - sp.conjugate(p)) / 2)) == 0
    checks["3g_even_antiherm_nonzero_for_real_p_ne_q"] = not zero(even_ah.subs({pI: 0, qI: 0, rI: 0, pR: 1, qR: 0}))

    # ---- (4) negative control: the linear involution (NOT Q14's J_Pi) ---------------------------------------
    O_lin = (L - lin(L)) / 2
    A_lin = antiherm(O_lin)
    checks["4a_linear_internal_entries_are_Im_p_plus_q"] = (
        sp.simplify(sp.expand_complex(A_lin[0, 1] - sp.I * s2 * (pI + qI) / 2)) == 0
        and sp.simplify(sp.expand_complex(A_lin[0, 2] - sp.I * s2 * (pI + qI) / 2)) == 0)
    checks["4b_linear_external_zero"] = zero(sp.Matrix([A_lin[1, 2], A_lin[2, 1]]))
    checks["4c_linear_diagonal_is_Im_r_rephasing"] = zero(
        sp.diag(A_lin[0, 0], A_lin[1, 1], A_lin[2, 2]) - sp.I * rI * sp.diag(0, 2, -2))
    zero_data = {pR: 0, qR: 0, rR: 0, rI: 0}
    checks["4d_linear_Im_p_nonzero_gives_transverse"] = not zero(
        sp.Matrix([A_lin[0, 1], A_lin[0, 2]]).subs({**zero_data, pI: 1, qI: 0}))
    checks["4e_linear_Im_p_equals_minus_Im_q_gives_zero"] = zero(A_lin.subs({**zero_data, pI: 1, qI: -1}))
    checks["4f_both_parities_agree_on_real_data"] = zero((odd(L) - O_lin).subs({pI: 0, qI: 0, rI: 0}))
    checks["4g_linear_differs_from_antilinear_off_real_data"] = not zero((lin(L) - jconj(L)).subs({pI: 1}))
    XR = sp.Matrix(3, 3, sp.symbols("a0:9", real=True))
    XI = sp.Matrix(3, 3, sp.symbols("b0:9", real=True))
    X = XR + sp.I * XI
    odd_anti = (X - jconj(X)) / 2
    odd_lin = (X - lin(X)) / 2
    checks["4h_parities_differ_exactly_on_imaginary_part"] = (
        zero(odd_anti - ((XR - lin(XR)) / 2 + sp.I * (XI + lin(XI)) / 2))
        and zero(odd_lin - ((XR - lin(XR)) / 2 + sp.I * (XI - lin(XI)) / 2)))

    # ---- (5) the vanishing is specific to the sl_2 image ----------------------------------------------------
    aR, aI, bR, bI, za, ya = sp.symbols("aR aI bR bI za ya", real=True)
    a, b = aR + sp.I * aI, bR + sp.I * bI
    Xoa = sp.Matrix([[sp.I * za, a, sp.conjugate(a)],
                     [-sp.conjugate(a), sp.I * ya, b],
                     [-a, -sp.conjugate(b), sp.I * ya]])
    checks["5a_family_is_antiherm_and_J_odd"] = zero(Xoa + Xoa.H) and zero(Xoa + jconj(Xoa))
    image = [sym2_lift(E), sym2_lift(F), sym2_lift(H)]
    checks["5b_family_HS_orthogonal_to_sl2_image"] = all(
        sp.simplify(sp.expand_complex((g.H * Xoa).trace())) == 0 for g in image)
    xs = sp.symbols("x0:18", real=True)
    Xg = sp.Matrix(3, 3, lambda i, j: xs[2 * (3 * i + j)] + sp.I * xs[2 * (3 * i + j) + 1])
    eqs = []
    for e in list(Xg + Xg.H) + list(Xg + jconj(Xg)):
        e = sp.expand_complex(e)
        eqs += [sp.re(e), sp.im(e)]
    dim_real = len(xs) - sp.Matrix([[sp.diff(e, x) for x in xs] for e in eqs]).rank()
    checks["5c_odd_antiherm_space_has_real_dimension_6"] = (dim_real == 6)
    # the internal block (parameter a) and R_mix (parameter b) are both in that space
    internal = Xoa.subs({aR: 1, aI: 0, bR: 0, bI: 0, za: 0, ya: 0})
    external = Xoa.subs({aR: 0, aI: 0, bR: 1, bI: 0, za: 0, ya: 0})
    checks["5d_internal_block_generator_exists_outside_image"] = (
        internal[0, 1] != 0 and zero(internal + internal.H) and zero(internal + jconj(internal))
        and sp.Matrix.hstack(*[sp.Matrix(list(g)) for g in image], sp.Matrix(list(internal))).rank() == 4)
    checks["5e_external_Rmix_generator_exists_outside_image"] = (
        zero(external - Rmix) and sp.Matrix.hstack(*[sp.Matrix(list(g)) for g in image],
                                                     sp.Matrix(list(external))).rank() == 4)

    # ---- (6) frame covariance ---------------------------------------------------------------------------------
    inv2 = 1 / s2
    W = sp.Matrix([[1, 0, 0], [0, inv2, sp.I * inv2], [0, sp.I * inv2, inv2]])
    checks["6a_W_unitary"] = zero(W * W.H - I3)
    Sp = W * S * W.T                              # transported antilinear parity z -> Sp conj(z)
    Lp = W * L * W.H
    odd_p = (Lp - Sp * Lp.conjugate() * Sp.inv()) / 2
    checks["6b_frame_has_complex_entries"] = any(
        sp.simplify(sp.im(e)) != 0 for e in (W * sym2_lift(E + F) * W.H))
    checks["6c_frame_generator_is_W_A_Wdagger_zero"] = zero(odd_p - W * odd(L) * W.H) and zero(antiherm(odd_p))

    # ---- (7) central imaginary scalar -------------------------------------------------------------------------
    cR, cI = sp.symbols("cR cI", real=True)
    central = (cR + sp.I * cI) * I3
    A_central = antiherm(odd(central))
    checks["7a_central_generator_is_i_Im_c_identity"] = zero(A_central - sp.I * cI * I3)
    checks["7b_central_generator_diagonal"] = all(sp.simplify(A_central[i, j]) == 0
                                                  for i in range(3) for j in range(3) if i != j)

    # ---- report -----------------------------------------------------------------------------------------------
    print("Front 3e: reality structure of the sl_2 step under the antilinear parity J_Pi of Q14 (exact symbolic)")
    print("=" * 104)
    print("  J_Pi z = S conj(z), S = Sym^2(eps); X -> S conj(X) S^{-1};")
    print("  A_Pi := antiherm(J_Pi-odd(L(M))), M in sl_2(C)")
    print("-" * 104)
    print("  (1) S conj(L(M)) S^{-1} = - L(M)^dagger for every complex (p, q, r)")
    print("  (2) J_Pi-odd part of L(M) = hermitian part; J_Pi-even part = anti-hermitian part")
    print("  (3) A_Pi = 0 identically in (p, q, r) in C^3; odd J_3 coefficient 2 Re r; odd R_mix coefficient 0")
    print("  (4) negative control, LINEAR involution S X S^{-1}: transverse part ~ Im(p+q), differs from J_Pi")
    print("  (5) odd anti-hermitian operators: 6 real dimensions, HS-orthogonal to the sl_2 image (spin-2 sector)")
    print("  (6) complex unitary frame with transported parity: W A_Pi W^dagger = 0 (vacuous)")
    print("  (7) central imaginary scalar: odd, anti-hermitian, diagonal (class-trivial)")
    print("-" * 104)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 104)
    print(f"  checks run: {len(checks)}")
    print("RESULT (proved, exact symbolic): A_Pi, the anti-hermitian part of the J_Pi-odd part of the sl_2 lift of the")
    print("  step generator, vanishes identically for all complex coefficients, because with Q14's antilinear J_Pi the")
    print("  J_Pi-odd part of the lift is its hermitian part. This is a statement about A_Pi as defined. The diagonal")
    print("  i Im(r) diag(0, 2, -2) is J_Pi-even. The J_Pi-even anti-hermitian part of the lift has the non-zero")
    print("  internal entry (sqrt2/2)(q - conj p), non-zero for real p != q, so the image of sl_2 does reach the")
    print("  internal block e_0 <-> e_+/- through that part; the external block R_mix is not reached by the image of")
    print("  sl_2 (Q14 Remark 6.4, Prop. 6.3). The linear involution is not J_Pi; it appears only as a negative")
    print("  control (condition Im(p+q)).")
    print("OPEN MODELLING QUESTION: no physical exclusion of the internal block is claimed. Whether A_Pi, rather than")
    print("  the J_Pi-even part, is the right object is a modelling choice of the companion note that no source")
    print("  justifies. No mass and no mixing value is produced.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
