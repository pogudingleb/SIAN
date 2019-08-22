read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [p1, p2, p3, p4, p5, p6, p7],
  x_vars = [x1_, x2_, x3_],
  y_vars = [y],
  u_vars = [u],
  x_eqs = [
    x1_1 = -1 * p1 * x1_0 + p2 * x2_0 + u0,
    x2_1 = p3 * x1_0 - p4 * x2_0 + p5 * x3_0,
    x3_1 = p6 * x1_0 - p7 * x3_0
  ],
  y_eqs = [
    y0 = x1_0
  ]
]);

theta_g := GlobalIdentifiability(sigma, [p1, p4, p7, x1_0, x2_0, x3_0], 0.99, 1);
print(theta_g);
