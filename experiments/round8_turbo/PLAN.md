# Round 8 plan (Turbo 6890080172): (T2) explicit constants, rho_j decay
Heart (Grok / GPPDiscovery5). Written before computing. L=0.5, l=1/4, N=12. No zero data.
- Exact spectral identity: psi_hat_j = S_j - w_j sin(tl) = (t/omega)^2 S_j  (since omega l = j pi).
- Second boundary functional: -delta''_l has F = t^2 sin(tl), sigma = -sinh(l/2)/4.  q_x = B(x, -delta''_l).
- Remainder: rho_hat_j = psi_hat_j - (w_j/omega^2) t^2 sin(tl) = (t/omega)^4 S_j, sigma_rho = sig(j)/(16 omega^4).
  So omega^4 rho_j(x) = (16/pi) int sin^4(t/2) V F_x S_j dt - (bI2^2/8) sigma_x sig(j)   (Psi^2 t^4 = 16 sin^4(t/2)).
- Old float growth hypothesis: truncation at T=8000 of the non-decaying part of psi_hat (~1/omega at t>>omega).
- Compute q, rho_j at j=73,123,173,273,472 with T up to 1e6 (chunked GL) + analytic mean tail; check T-convergence;
  check v_j = (w_j/omega^2) q + rho_j against direct v_j.
- Summability for (T3): |v_j| = O(omega^-3) suffices since D_j ~ log(omega)/omega^4.
Backup: if tail correction does not converge, use kernel-domain IBP form rho_j = -sqrt(2/l) omega^-3 int h''' cos.
