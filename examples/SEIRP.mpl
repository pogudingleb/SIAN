read "../IdentifiabilityODE.mpl";

sys := [
  diff(s(t), t) = - a_e * s(t) * e(t) - a_i * s(t) * i(t),
  diff(e(t), t) = a_e * s(t) * e(t) + a_i * s(t) * i(t) - k * e(t) - rho * e(t),
  diff(i(t), t) = k * e(t) - b * i(t) - mu * i(t),
  diff(r(t), t) = b * i(t) + rho * e(t),
  diff(p(t), t) = mu * i(t),
  y1(t) = i(t) + s(t)
];

IdentifiabilityODE(sys, GetParameters(sys), substitute_tr_basis=true):

