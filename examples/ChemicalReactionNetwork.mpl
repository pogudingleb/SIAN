# Example 9 in the paper
# Taken from 
# Conradi, C., Shiu, A.,
# Dynamics of post-translational modification systems: recent progress and future directions
# Eq. 3.4
read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [k1, k2, k3, k4, k5, k6],
  x_vars = [x1_, x2_, x3_, x4_, x5_, x6_],
  y_vars = [y1_, y2_],
  u_vars = [],
  x_eqs = [
    x1_1 = -k1 * x1_0 * x2_0 + k2 * x4_0 + k4 * x6_0,
    x2_1 = k1 * x1_0 * x2_0 + k2 * x4_0 + k3 * x4_0,
    x3_1 = k3 * x4_0 + k5 * x6_0 - k6 * x3_0 * x5_0,
    x4_1 = k1 * x1_0 * x2_0 - k2 * x4_0 - k3 * x4_0,
    x5_1 = k4 * x6_0 + k5 * x6_0 - k6 * x3_0 * x5_0,
    x6_1 = -k4 * x6_0 - k5 * x6_0 + k6 * x3_0 * x5_0
  ],
  y_eqs = [
    y1_0 = x3_0,
    y2_0 = x2_0
  ]
]);

theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), x1_0, x2_0, x3_0, x4_0, x5_0, x6_0], 0.99);
print(theta_g);
