# Taken from
# Capistran M., Moreles M., Lara B.
# "Parameter Estimation of Some Epidemic Models. The Case of Recurrent Epidemics Caused by Respiratory Syncytial Virus"
# doi.org/10.1007/s11538-009-9429-3
# Equations (7)-(11)
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(s(t), t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
  diff(i(t), t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
  diff(r(t), t) = nu * i(t) - (mu + g) * r(t),
  diff(x1(t), t) = -M * x2(t),
  diff(x2(t), t) = M * x1(t),
  y1(t) = i(t),
  y2(t) = r(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
