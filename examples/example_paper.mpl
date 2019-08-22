read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [mu1, mu2],
  x_vars = [x],
  y_vars = [y],
  u_vars = [],
  x_eqs = [
    x1 = mu2 * x0 + mu1
  ],
  y_eqs = [
    y0 = x0^2
  ]
]);

theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), x0], 0.8);
print(theta_g);
