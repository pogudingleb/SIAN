# Taken from
# N. Tunker, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.3)
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(S(t), t) = -b * S(t) * In(t) / (S(t) + Tr(t) + In(t) + R(t)) - d * b * S(t) * Tr(t) / (S(t) + Tr(t) + In(t) + R(t)),
  diff(In(t), t) = b * S(t) * In(t) / (S(t) + Tr(t) + In(t) + R(t)) + d * b * S(t) * Tr(t) / (S(t) + Tr(t) + In(t) + R(t)) - (a + g) * In(t),
  diff(Tr(t), t) = g * In(t) - nu * Tr(t),
  diff(R(t), t) = a * In(t) + nu * Tr(t),
  y(t) = Tr(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
