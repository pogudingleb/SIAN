# Taken from
# N. Tunker, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.3) with observed treatment
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(S(t), t) = -b * S(t) * In(t) / N(t) - d * b * S(t) * Tr(t) / N(t),
  diff(In(t), t) = b * S(t) * In(t) / N(t) + d * b * S(t) * Tr(t) / N(t) - (a + g) * In(t),
  diff(Tr(t), t) = g * In(t) - nu * Tr(t),
  diff(N(t), t) = 0,
  y1(t) = Tr(t),
  y2(t) = N(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
