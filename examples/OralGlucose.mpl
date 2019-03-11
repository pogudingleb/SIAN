# Example (with initial conditions assumed being unknown) from Section III of "DAISY: an Efficient Tool to Test Global Identifiability. Some Case Studies"
# by G. Bellu, M.P. Saccomani

read "../IdentifiabilityODE.mpl";

sigma := [
  diff(G(t), t) = -(p1 + X(t)) * G(t) + p1 * Gb(t) + v * R(t),
  diff(X(t), t) = -p2 * X(t) + p3 * (u(t) - Ib(t)),
  diff(R(t), t) = k,
  diff(Ib(t), t) = 0,
  diff(Gb(t), t) = 0,
  y1(t) = G(t),
  y2(t) = Ib(t),
  y3(t) = Gb(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma));
