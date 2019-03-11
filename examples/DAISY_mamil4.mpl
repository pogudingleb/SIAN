read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = -k01 * x1(t) + k12 * x2(t) + k13 * x3(t) + k14 * x4(t) - k21 * x1(t) - k31 * x1(t) - k41 * x1(t) + u(t),
  diff(x2(t), t) = -k12 * x2(t) + k21 * x1(t),
  diff(x3(t), t) = -k13 * x3(t) + k31 * x1(t),
  diff(x4(t), t) = -k14 * x4(t) + k41 * x1(t),
  y(t) = x1(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
