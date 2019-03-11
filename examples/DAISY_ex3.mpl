read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = -1 * p1 * x1(t) + x2(t) + u0(t),
  diff(x2(t), t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  diff(x3(t), t) = p6 * x1(t) - p7 * x3(t),
  diff(u0(t), t) = 1,
  y(t) = x1(t),
  y2(t) = u0(t)
 
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
