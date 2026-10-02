r"""Front 3c -> Front 2 bridge: transverse route of the metaplectic step generator.

Exact symbolic verification (no sampling). The audit of the chiral polar class (front3c_polar_class_nontriviality.py)
reduces the non-triviality of the class [U_Pi] to the OFF-DIAGONAL (transverse) part of the anti-hermitian generator A
of the polar factor. This script asks which transverse blocks the Sym^2 lift of an sl_2(C) step can source.

Generator (definition, taken as given). A_Pi := antiherm( J_Pi-odd( L(M) ) ), M = p E + q F + r H in sl_2(C) with
complex coefficients, L(M) the derived Sym^2 lift on C^3_gen, and J_Pi the ANTILINEAR parity of Q14 Section 6,
J_Pi z = S conj(z) with S = Sym^2(eps), eps = [[0, 1], [-1, 0]] (J_Pi: e_0 -> -e_0, e_+ <-> e_- antilinearly). The
conjugation of an operator is X -> S conj(X) S^{-1}. The exact audit of this definition is
front3e_reality_structure.py; this script reuses it on the blocks.

Two transverse blocks on C^3_gen = Sym^2(C^2), basis (e_0, e_+, e_-):
    INTERNAL block (e_0 <-> e_+, e_0 <-> e_-): central generation mixing with the outer pair;
    EXTERNAL block (e_+ <-> e_-) = R_mix:      mixing of the two outer generations (antisymmetric, Q14 Sec. 6).

Results (all exact symbolic).
  (A) For a real sl_2 element A_Pi = 0 (both blocks).
  (B) For a complex sl_2 element A_Pi = 0 as well: the internal entry and the external entry vanish for every
      complex (p, q, r). The J_Pi-odd part of L(M) is the hermitian part of L(M).
  (C) The N_A channel: the J_3 coefficient of the J_Pi-odd part is 2 Re r; on the cascade step
      g = exp(tE) exp(sF) it is the oriented area alpha = t s (AAR). It lies in the diagonal channel.
  (D) The lift has identically zero external entries for any complex M, and R_mix is linearly independent of
      span{L(E), L(F), L(H), I_3} (rank 4 -> 5). The internal-block generator is also outside the sl_2 image
      (rank 3 -> 4) and Hilbert-Schmidt orthogonal to it. Both blocks need a generator in the spin-2 sector of
      End(Sym^2 V_gen), outside the sl_2(C) image; none is supplied.
  (E) A central scalar c I_3: its J_Pi-odd anti-hermitian part is i Im(c) I_3, diagonal, hence class-trivial.
  (F) A diagonal generator diag(d_0, d_+, d_-) (complex) has a zero off-diagonal J_Pi-odd anti-hermitian part.

Verdict (printed, scoped to A_Pi as defined). Within the sl_2 image neither transverse block is sourced, for real or
complex coefficients. Q14 Remark 6.4 concerns R_mix and rests on the absence of a matrix element between the weights
+1 and -1; the vanishing of the internal block rests on the antilinear parity (front3e). A source of either block
needs a generator outside the sl_2 image. Whether the polar generator should be the J_Pi-odd part is a property of the
definition of A_Pi that is not justified by the sources used here (see front3e). No mass and no mixing value is
produced. No figures. English.
"""

import sympy as sp

s2 = sp.sqrt(2)


def fundamental_generators():
    E = sp.Matrix([[0, 1], [0, 0]])
    F = sp.Matrix([[0, 0], [1, 0]])
    H = sp.Matrix([[1, 0], [0, -1]])
    return E, F, H


def sym2_of_matrix(g):
    """Group action Sym^2(g) on the orthonormal basis (e_0, e_+, e_-); T -> g T g^T on symmetric 2x2 matrices."""
    basis = [sp.Matrix([[0, 1], [1, 0]]) / s2, sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    cols = []
    for T in basis:
        Tn = g * T * g.T
        cols.append([Tn[0, 1] * s2, Tn[0, 0], Tn[1, 1]])
    return sp.Matrix(cols).T


def sym2_lift(M):
    """Derived Sym^2 representation of a 2x2 matrix M: d/dt Sym^2(1 + t M) at t = 0."""
    t = sp.Symbol("t")
    return sp.simplify(sym2_of_matrix(sp.eye(2) + t * M).diff(t).subs(t, 0))


S = sym2_of_matrix(sp.Matrix([[0, 1], [-1, 0]]))      # Sym^2(eps)


def jpi_odd_part(A):
    """J_Pi-odd part (A - S conj(A) S^{-1}) / 2 for the antilinear J_Pi z = S conj(z)."""
    return (A - S * A.conjugate() * S.inv()) / 2


def antiherm(A):
    return (A - A.H) / 2


def Upi_generator(M):
    """A_Pi: J_Pi-odd AND anti-hermitian part of the Sym^2 lift of M."""
    return antiherm(jpi_odd_part(sym2_lift(M)))


def hs_inner(A, B):
    return sp.trace(A.H * B)


def vec(M):
    return sp.Matrix([M[i, j] for i in range(3) for j in range(3)])


def main():
    checks = {}
    E, F, H = fundamental_generators()
    J3 = sp.diag(0, 1, -1)
    Rmix = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
    I3 = sp.eye(3)
    p, q, r = sp.symbols("p q r")                       # generic COMPLEX sl_2 coefficients
    pR, pI, qR, qI, rR, rI = sp.symbols("pR pI qR qI rR rI", real=True)
    rep = {p: pR + sp.I * pI, q: qR + sp.I * qI, r: rR + sp.I * rI}

    def is_zero(X):
        return all(sp.simplify(sp.expand_complex(e)) == 0 for e in X)

    # ---- (A) real step: A_Pi = 0 -------------------------------------------------------------------------------
    A_real = Upi_generator(pR * E + qR * F + rR * H)
    checks["A_real_generator_zero"] = is_zero(A_real)

    # ---- (B) complex step: A_Pi = 0 as well (both blocks) ------------------------------------------------------
    A_cpx = Upi_generator((p * E + q * F + r * H).subs(rep))
    checks["B_internal_entries_zero_any_complex_M"] = is_zero(sp.Matrix([A_cpx[0, 1], A_cpx[0, 2], A_cpx[1, 0],
                                                                         A_cpx[2, 0]]))
    checks["B_external_entries_zero_any_complex_M"] = is_zero(sp.Matrix([A_cpx[1, 2], A_cpx[2, 1]]))
    checks["B_whole_generator_zero_any_complex_M"] = is_zero(A_cpx)
    L = sym2_lift((p * E + q * F + r * H).subs(rep))
    checks["B_odd_part_is_hermitian_part"] = is_zero(jpi_odd_part(L) - (L + L.H) / 2)

    # ---- (C) N_A channel: J_3 coefficient of the J_Pi-odd part -------------------------------------------------
    t, s = sp.symbols("t s", real=True)
    odd_step = jpi_odd_part(sym2_lift(t * E + s * F + (t * s / 2) * H))
    alpha = sp.simplify(hs_inner(J3, odd_step) / hs_inner(J3, J3))
    checks["C_NA_is_real_area_ts"] = sp.simplify(alpha - t * s) == 0
    odd_cpx = jpi_odd_part(L)
    checks["C_J3_coefficient_is_2Re_r"] = sp.simplify(
        sp.expand_complex(hs_inner(J3, odd_cpx) / hs_inner(J3, J3) - 2 * rR)) == 0

    # ---- (D) both blocks need a generator outside the sl_2 image -----------------------------------------------
    LE, LF, LH = sym2_lift(E), sym2_lift(F), sym2_lift(H)
    full = sym2_lift(p * E + q * F + r * H)
    checks["D_lift_external_entries_zero"] = sp.simplify(full[1, 2]) == 0 and sp.simplify(full[2, 1]) == 0
    basis4 = sp.Matrix.hstack(vec(LE), vec(LF), vec(LH), vec(I3))
    basis5 = sp.Matrix.hstack(basis4, vec(Rmix))
    checks["D_Rmix_outside_image"] = basis4.rank() == 4 and basis5.rank() == 5
    checks["D_J3_perp_Rmix"] = sp.simplify(hs_inner(J3, Rmix)) == 0
    internal = sp.Matrix([[0, 1, 1], [-1, 0, 0], [-1, 0, 0]])          # J_Pi-odd, anti-hermitian, internal block
    checks["D_internal_generator_is_odd_antiherm"] = (
        is_zero(internal + internal.H) and is_zero(internal + S * internal.conjugate() * S.inv()))
    basis3 = sp.Matrix.hstack(vec(LE), vec(LF), vec(LH))
    checks["D_internal_outside_image"] = basis3.rank() == 3 and sp.Matrix.hstack(basis3, vec(internal)).rank() == 4
    checks["D_internal_HS_orthogonal_to_image"] = all(
        sp.simplify(hs_inner(g, internal)) == 0 for g in (LE, LF, LH))

    # ---- (E) central scalar --------------------------------------------------------------------------------------
    cR, cI = sp.symbols("cR cI", real=True)
    A_central = antiherm(jpi_odd_part((cR + sp.I * cI) * I3))
    checks["E_central_generator_is_i_Im_c_identity"] = is_zero(A_central - sp.I * cI * I3)

    # ---- (F) diagonal data: no transverse part -------------------------------------------------------------------
    d0, dp, dm = sp.symbols("d0 dp dm")
    A_diag = antiherm(jpi_odd_part(sp.diag(d0, dp, dm)))
    checks["F_diag_no_transverse"] = all(sp.simplify(sp.expand_complex(A_diag[i, j])) == 0
                                         for i in range(3) for j in range(3) if i != j)

    # ---- report ----------------------------------------------------------------------------------------------------
    print("Front 3c -> Front 2 bridge: transverse route of the sl_2 step generator (exact symbolic)")
    print("=" * 100)
    print("  A_Pi = antiherm(J_Pi-odd(L(M))), J_Pi antilinear (z -> S conj z), M in sl_2(C);")
    print("  transverse blocks: INTERNAL (e_0 <-> e_+/-), EXTERNAL R_mix (e_+ <-> e_-)")
    print("  (A) real step: A_Pi = 0 (both blocks)")
    print("  (B) complex step: A_Pi = 0 as well; the J_Pi-odd part of the lift is its hermitian part")
    print("  (C) N_A: the J_3 coefficient of the odd part is 2 Re r (= t s on the cascade step), a diagonal channel")
    print("  (D) R_mix and the internal block are outside the sl_2 image (ranks 4 -> 5, 3 -> 4), HS-orthogonal to it")
    print("  (E) central scalar: i Im(c) I_3, diagonal, class-trivial")
    print("  (F) diagonal data: zero off-diagonal generator")
    print("-" * 100)
    allok = True
    for k, val in checks.items():
        ok = bool(val)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print("=" * 100)
    print(f"  checks run: {len(checks)}")
    print("VERDICT (scoped to A_Pi as defined): within the sl_2 image neither transverse block is sourced, for real")
    print("  or complex coefficients. Q14 Remark 6.4 concerns R_mix only; the vanishing of the internal block rests on")
    print("  the antilinear parity (front3e). A source of either block needs a generator outside the sl_2 image, in")
    print("  the spin-2 sector of End(Sym^2 V_gen); none is supplied. Whether the polar generator should be the")
    print("  J_Pi-odd part is a property of the definition, not justified by the sources used here (front3e).")
    print("  No mass and no mixing value is produced.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
