# Sign check: does additive (Poisson) positivity transfer through Connes' E to Weil positivity? (Queen, 2026-10-28)

Daniel Toupin / Golden Physics Project, with AI assistance from Grok.
RH is OPEN. This note proves nothing about RH; it fixes a sign/space in the GAP statement.

Sources read this session: Connes 1999, "Trace formula in noncommutative geometry and the zeros of the
Riemann zeta function", Selecta Math. 5 (1999) 29-106, arXiv math/9811068: Sec. III (18)-(19), (31)-(34);
Sec. V Thm 3; Sec. VII Thm 4; Sec. VIII Lemma 1, Cor. 2, (16), Thm 5 with proof (17)-(26), Lemma 3.
Meyer, Duke Math. J. 127 (2005) 519-595, arXiv math/0311468: abstract only (virtual representation whose
character is the Weil distribution); no positivity statement is taken from Meyer.
CCM: Connes-Consani-Moscovici arXiv:2310.18423 (Ann. Funct. Anal. 2024), intro: for each window the
semilocal trace formula writes Weil's Q_n as the trace of an operator expression.

## Connes' positivity statement (proof of Thm 5, a => b), written out

- E : L^2(X)_0 -> L^2(C_k) is an isometry, E U(a) = |a|^{1/2} V(a) E, (Ef)(g) = |g|^{1/2} sum_{q in k^*} f(qg).
- S_Lambda = projection onto {xi in L^2(C_k) : xi(g) = 0 unless |g| in [1/Lambda, Lambda]}.
- B_{Lambda,0} = span of f in S(A)_0 with f and f^ vanishing for |x| > Lambda; Q'_{Lambda,0} = E Q_{Lambda,0} E^{-1}.
- Since E(f)(g) = E(f^)(1/g) (this is Theta = F_A o inv transported by E), E(B_{Lambda,0}) is inside S_Lambda,
  so Q'_{Lambda,0} <= S_Lambda (operator inequality between projections).
- Hence Delta_Lambda(h) := Tr((S_Lambda - Q'_{Lambda,0}) V(h)) is of positive type: Delta_Lambda(h * h^*) >= 0
  for every Lambda. [Connes (23)-(25)]
- Tr(S_Lambda V(h)) = 2 h(1) log' Lambda exactly. So Delta_Lambda = (window trace) - (trace on the E-range
  of the doubly-localized space).
- The global trace formula (16) says exactly that lim Delta_Lambda = Weil distribution. Thm 5: (16) <=> RH
  (function fields: proof written; for Q Connes states the analogues with the prolate B_Lambda, Sec. VIII end).
- Lemma 3 (function fields, written): lim Delta_Lambda(f) = sum over ALL zeros of int f^ d mu_rho, mu_rho = harmonic
  measure (Poisson kernel) of rho w.r.t. the critical line; it equals the point mass at rho iff rho is on the line.

## Verdict

1. Additive positivity lives on the RANGE of E: the Poisson comb is positive-definite, so the range projection
   Q'_{Lambda,0} >= 0. That range enters Delta_Lambda with a MINUS sign. Transferring additive positivity
   through E therefore gives only the trivial side (Tr Q' V(h)V(h)^* >= 0), which bounds the Weil form from
   ABOVE by the window term, not from below. "Additive positivity survives the transfer": WRONG SIGN.
2. The positive object is the COKERNEL/complement within the window, S_Lambda - Q'_{Lambda,0}. Its positivity is
   free (projection inequality, using only Theta-symmetry of the doubly-localized space). PROVED (Connes).
   For Q it is stated by Connes with prolate B_Lambda (only approximately space-and-band limited): stated, not
   re-checked here.
3. So the GAP is NOT a positivity statement in Connes' picture; it is an IDENTIFICATION: Weil = lim Delta_Lambda,
   i.e. the residual delta_Lambda(f) := Delta_Lambda(f*f^*) - Q(f) -> 0. By Lemma 3 the residual in the limit is
   exactly sum_{rho off-line} m_rho [ int g^ d mu_rho - g^(gamma_rho) ]: off-line zeros enter the positive
   side smeared by a Poisson kernel ("resonances") and the Weil side as point evaluations; the gap between the two is
   the whole of RH.
4. Consistency with ARP G.4d: the GAP "Fl_L >= -MF_L uniformly in L" and "delta_Lambda -> 0" are both RH-equivalent,
   so they are consistent. Two wording corrections to G.4d:
   (a) the positive part is Tr((S_Lambda - Q'_{Lambda,0}) V(h)V(h)^*), a trace over the complement of the range,
       not "||P_Lambda^perp E(phi)||^2" (that is a norm of a RANGE vector, the wrong space);
   (b) the 2 log Lambda term is not part of delta_Lambda: it is cancelled exactly by Tr(S_Lambda V(h)).
   No identification of MF_L/Fl_L with Connes' (Delta_Lambda, delta_Lambda) is claimed.

Labels: items 1-2 PROVED (Connes 1999, function-field proof checked line by line; Q case as stated by Connes);
item 3 PROVED for function fields (Connes Thm 5 + Lemma 3), for Q as stated by Connes (Sec. VIII, last part), not
re-derived; item 4 is bookkeeping (PROVED, elementary) plus a proposed wording change.
