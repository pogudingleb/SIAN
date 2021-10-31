read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = -1 * p1 * x1(t) + p2 * x2(t) + u(t),
  diff(x2(t), t) = p3 * x1(t) - p4 * x2(t) + p5 * x3(t),
  diff(x3(t), t) = p6 * x1(t) - p7 * x3(t),
  y(t) = x1(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), sub_transc=true):
