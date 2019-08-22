# Example (with initial conditions assumed being unknown) from Section IV of "DAISY: an Efficient Tool to Test Global Identifiability. Some Case Studies"
# by G. Bellu, M.P. Saccomani

read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [beta, q1, q2, mu1, mu2, k1, k2, c, s, d],
  x_vars = [x1_, x2_, x3_, x4_],
  y_vars = [y1_,y2_],
  u_vars = [],
  x_eqs = [
    x1_1 = -beta * x1_0 * x4_0 - d * x1_0 + s,
    x2_1 = beta * q1 * x1_0 * x4_0 - k1 * x2_0 - mu1 * x2_0,
    x3_1 = beta * q2 * x1_0 * x4_0 + k1 * x2_0 - mu2 * x3_0,
    x4_1 = -c * x4_0 + k2 * x3_0
  ],
  y_eqs = [
    y1_0 = x1_0,
    y2_0 = x4_0
  ]
]);

# only beta, s, d, x1_0, x4_0 are locally identifiable
theta_g := GlobalIdentifiability(sigma, [beta, s, d, x1_0, x4_0], 0.99);
print(theta_g);
