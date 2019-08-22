read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [k01, k12, k21, k13, k31, k14, k41],
  x_vars = [x1_, x2_, x3_, x4_],
  y_vars = [y],
  u_vars = [u],
  x_eqs = [
    x1_1 = -k01 * x1_0 + k12 * x2_0 + k13 * x3_0 + k14 * x4_0 - k21 * x1_0 - k31 * x1_0 - k41 * x1_0 + u0,
    x2_1 = -k12 * x2_0 + k21 * x1_0,
    x3_1 = -k13 * x3_0 + k31 * x1_0,
    x4_1 = -k14 * x4_0 + k41 * x1_0
  ],
  y_eqs = [
    y0 = x1_0
  ]
]);

theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), x1_0, x2_0, x3_0, x4_0], 0.99, 1);
print(theta_g);
