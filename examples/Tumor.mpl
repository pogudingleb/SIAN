# Example (with initial conditions assumed being unknown) from Section 3 of "Examples of testing global identifiability of biological and biomedical models with the DAISY software"
# by M.P. Saccomani, S. Audoly, G. Bellu, L. D'Angio

read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [k3, k4, k5, k6, k7],
  x_vars = [x1_, x2_, x3_, x4_, x5_, a, b, d],
  y_vars = [y1_, y2_, y3_, y4_],
  u_vars = [],
  x_eqs = [
    x1_1 = -(k3 + k7) * x1_0 + k4 * x2_0,
    x2_1 = k3 * x1_0 - (k4 + a0 * k5 + b0 * d0 * k5) * x2_0 + k6 * x3_0 + k6 * x4_0 + k5 * x2_0 * x3_0 + k5 * x2_0 * x4_0,
    x3_1 = a0 * k5 * x2_0 - k6 * x3_0 - k5 * x2_0 * x3_0,
    x4_1 = b0 * d0 * k5 * x2_0 - k6 * x4_0 - k5 * x2_0 * x4_0,
    x5_1 = k7 * x1_0,
    a1 = 0,
    b1 = 0,
    d1 = 0
  ],
  y_eqs = [
    y1_0 = x5_0,
    y2_0 = a0,
    y3_0 = b0,
    y4_0 = d0
  ]
]);

# x3_0 and x4_0 are not locally identifiable
theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), x1_0, x2_0, x5_0], 0.99);
print(theta_g);
