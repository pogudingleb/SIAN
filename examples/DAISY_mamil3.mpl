read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = -(a21 + a31 + a01) * x1(t) + a12 * x2(t) + a13 * x3(t) + u(t),
  diff(x2(t), t) = a21 * x1(t) - a12 * x2(t),
  diff(x3(t), t) = a31 * x1(t) - a13 * x3(t),
  y(t) = x1(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
