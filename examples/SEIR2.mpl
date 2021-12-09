# Taken from
# N. Tuncer, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.2) with cumulative incidence observations
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(S(t), t) = -b * S(t) * In(t) / N,
  diff(E(t), t) = b * S(t) * In(t) / N - nu * E(t),
  diff(In(t), t) = nu * E(t) - a * In(t),
  diff(Cu(t), t) = nu * E(t),
  y1(t) = Cu(t),
  y2(t) = N
];

IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true):
