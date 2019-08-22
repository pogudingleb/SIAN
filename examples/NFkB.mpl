# Example 11 from the paper, taken from
# Balsa-Canto, E., Alonso, A. A., Banga, J. R., 
# An iterative identification procedure for dynamic modeling of biochemical networks
read "../GlobalIdentifiability.mpl";

# Data from Table 1 in 
# An iterative identification procedure for dynamic modeling of biochemical networks
known_data := [
  a1 = 1 / 2,
  a2 = 1 / 5,
  a3 = 1,
  c_1a = 5 / 10^(7),
  c_2a = 0,
  c_5a = 1 / 10^(4),
  c_6a = 2 / 10^(5),
  c1 = 5 / 10^(7),
  c2 = 0,
  c3 = 4 / 10^(4),
  c4 = 1 / 2,
  kv = 5,
  e_1a = 5 / 10^(4),
  c_1c = 5 / 10^(7),
  c_2c = 0,
  c_3c = 4 / 10^(4)
];

sigma := table([
  mu = [t1, t2, c_3a, c_4a, c5, k1, k2, k3, k_prod, k_deg, i1, e_2a, i_1a],
  x_vars = [x1_, x2_, x3_, x4_, x5_, x6_, x7_, x8_, x9_, x10_, x11_, x12_, x13_, x14_, x15_],
  y_vars = [y1_, y2_, y3_, y4_, y6_, y5_],
  u_vars = [u],
  x_eqs = subs(known_data, [
    x1_1 = k_prod - k_deg * x1_0 - k1 * x1_0 * u0,
    x2_1 = -k3 * x2_0 - k_deg * x2_0 - a2 * x2_0 * x10_0 + t1 * x4_0 - a3 * x2_0 * x13_0 + t2 * x5_0 + (k1 * x1_0 - k2 * x2_0 * x8_0) * u0,
    x3_1 = k3 * x2_0 - k_deg * x3_0 + k2 * x2_0 * x8_0 * u0,
    x4_1 = a2 * x2_0 * x10_0 - t1 * x4_0,
    x5_1 = a3 * x2_0 * x13_0 - t2 * x5_0,
    x6_1 = c_6a * x13_0 - a1 * x6_0 * x10_0 + t2 * x5_0 - i1 * x6_0,
    x7_1 = i1 * kv * x6_0 - a1 * x11_0 * x7_0,
    x8_1 = c4 * x9_0 - c5 * x8_0,
    x9_1 = c2 + c1 * x7_0 - c3 * x9_0,
    x10_1 = -a2 * x2_0 * x10_0 - a1 * x10_0 * x6_0 + c_4a * x12_0 - c_5a * x10_0 - i_1a * x10_0 + e_1a * x11_0,
    x11_1 = -a1 * x11_0 * x7_0 + i_1a * kv * x10_0 - e_1a * kv * x11_0,
    x12_1 = c_2a + c_1a * x7_0 - c_3a * x12_0,
    x13_1 = a1 * x10_0 * x6_0 - c_6a * x13_0 - a3 * x2_0 * x13_0 + e_2a * x14_0,
    x14_1 = a1 * x11_0 * x7_0 - e_2a * kv * x14_0,
    x15_1 = c_2c + c_1c * x7_0 - c_3c * x15_0
  ]),
  y_eqs = [
    y1_0 = x2_0,
    y2_0 = x10_0 + x13_0,
    y3_0 = x9_0,
    y4_0 = x1_0 + x2_0 + x3_0,
    y6_0 = x12_0,
    y5_0 = x7_0
  ]
]);

theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0, x10_0, x11_0, x12_0, x13_0, x14_0], 0.99);

print(theta_g);

