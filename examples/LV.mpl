read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = a * x1(t) - b * x1(t) * x2(t),
  diff(x2(t), t) = -c * x2(t) + d * x1(t) * x2(t),
  y(t) = x1(t)
];

<<<<<<< HEAD
IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true, p=0.99);
||||||| e3ff1e1
IdentifiabilityODE(sigma, [a, b, c, d, x1(0), x2(0)]);
=======
IdentifiabilityODE(sigma, GetParameters(sigma),p=0.01);
>>>>>>> master
