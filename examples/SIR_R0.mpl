read "../IdentifiabilityODE.mpl";

sigma := [
  diff(S(t), t) = -b * In(t) * S(t),
  diff(In(t), t) = b * In(t) * S(t) - g * In(t),
  diff(R(t), t) = g * In(t),
  diff(aux(t), t) = 0,
  y1(t) = In(t),
  y2(t) = b / g + aux(t)
];

IdentifiabilityODE(sigma, [aux(0)], infolevel = 3):
