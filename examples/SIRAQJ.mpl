read "../IdentifiabilityODE.mpl";

sys := [
  diff(S(t), t) = b * N(t) - S(t) * (In(t) * lam + lam * Q(t) * eps_a * eps_q + lam * eps_a * A(t) + lam * eps_j * Jj(t) + d1),
  diff(In(t), t) = k1 * A(t) - (g1 + mu2 + d2) * In(t), 
  diff(R(t), t) = g1 * In(t) + g2 * Jj(t) - d3 * R(t),
  diff(A(t), t) = S(t) * (In(t) * lam + lam * Q(t) * eps_a * eps_q + lam * eps_a * A(t) + lam * eps_j * Jj(t)) - (k1 + mu1 + d4) * A(t),
  diff(Q(t), t) = mu1 * A(t) - (k2 + d5) * Q(t),
  diff(Jj(t), t) = k2 * Q(t) + mu2 * In(t) - (g2 + d6) * Jj(t), 
  diff(N(t), t) = 0,
  y1(t) = Q(t),
  y2(t) = Jj(t)
];

IdentifiabilityODE(sys, GetParameters(sys), substitute_tr_basis=true, optimize_tr_basis=true, infolevel=2):

