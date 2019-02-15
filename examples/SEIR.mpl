# Taken from
# N. Tunker, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.2)
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(S(t), t) = -b * S(t) * In(t) / (S(t) + E(t) + In(t) + R(t)),
  diff(E(t), t) = b * S(t) * In(t) / (S(t) + E(t) + In(t) + R(t)) - nu * E(t),
  diff(In(t), t) = nu * E(t) - a * In(t),
  diff(R(t), t) = a * In(t),
  y(t) = In(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
