read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = x2(t) + a,
  diff(x2(t), t) = 0,
  y(t) = x1(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), known_initial_values = [x2(0)], infolevel = 2);
