read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [a12, a13, a21, a31, a01],
  x_vars = [x1_, x2_, x3_],
  y_vars = [y],
  u_vars = [u],
  x_eqs = [
    x1_1 = -(a21 + a31 + a01) * x1_0 + a12 * x2_0 + a13 * x3_0 + u0,
    x2_1 = a21 * x1_0 - a12 * x2_0,
    x3_1 = a31 * x1_0 - a13 * x3_0
  ],
  y_eqs = [
    y0 = x1_0
  ]
]);

theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), x1_0, x2_0, x3_0], 0.99, 2);
print(theta_g);
