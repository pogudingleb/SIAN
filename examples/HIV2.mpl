# The system is taken from
# Wodarz, D., Nowak, M.
# Specific therapy regimes could lead to long-term immunological control of HIV
# https://doi.org/10.1073/pnas.96.25.14464
# Page 1
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(x(t), t) = lm - d * x(t) - beta * x(t) * v(t),
  diff(y(t), t) = beta * x(t) * v(t) - a * y(t),
  diff(v(t), t) = k * y(t) - u * v(t),
  diff(w(t), t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
  diff(z(t), t) = c * q * y(t) * w(t) - h * z(t),
  y1(t) = w(t),
  y2(t) = z(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), weighted_ordering=true, infolevel=2):
