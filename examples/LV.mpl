read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = a * x1(t) - b * x1(t) * x2(t),
  diff(x2(t), t) = -c * x2(t) + d * x1(t) * x2(t),
  y(t) = x1(t) + u(t)
];

IdentifiabilityODE(sigma, [a, b, c, d, x1(0), x2(0)]);
