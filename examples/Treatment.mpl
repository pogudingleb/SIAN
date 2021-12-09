# Taken from
# N. Tucker, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.3) with observed treatment
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(S(t), t) = -b * S(t) * In(t) / N - d * b * S(t) * Tr(t) / N,
  diff(In(t), t) = b * S(t) * In(t) / N + d * b * S(t) * Tr(t) / N - (a + g) * In(t),
  diff(Tr(t), t) = g * In(t) - nu * Tr(t),
  y1(t) = Tr(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true):