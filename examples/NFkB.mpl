# Example 6.3 from the paper "Global Identifiability of Differential Models", taken from
# Balsa-Canto, E., Alonso, A. A., Banga, J. R., 
# An iterative identification procedure for dynamic modeling of biochemical networks
read "../IdentifiabilityODE.mpl";

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

sigma := subs(known_data, [
  diff(x1(t), t) = k_prod - k_deg * x1(t) - k1 * x1(t) * u(t),
  diff(x2(t), t) = -k3 * x2(t) - k_deg * x2(t) - a2 * x2(t) * x10(t) + t1 * x4(t) - a3 * x2(t) * x13(t) + t2 * x5(t) + (k1 * x1(t) - k2 * x2(t) * x8(t)) * u(t),
  diff(x3(t), t) = k3 * x2(t) - k_deg * x3(t) + k2 * x2(t) * x8(t) * u(t),
  diff(x4(t), t) = a2 * x2(t) * x10(t) - t1 * x4(t),
  diff(x5(t), t) = a3 * x2(t) * x13(t) - t2 * x5(t),
  diff(x6(t), t) = c_6a * x13(t) - a1 * x6(t) * x10(t) + t2 * x5(t) - i1 * x6(t),
  diff(x7(t), t) = i1 * kv * x6(t) - a1 * x11(t) * x7(t),
  diff(x8(t), t) = c4 * x9(t) - c5 * x8(t),
  diff(x9(t), t) = c2 + c1 * x7(t) - c3 * x9(t),
  diff(x10(t), t) = -a2 * x2(t) * x10(t) - a1 * x10(t) * x6(t) + c_4a * x12(t) - c_5a * x10(t) - i_1a * x10(t) + e_1a * x11(t),
  diff(x11(t), t) = -a1 * x11(t) * x7(t) + i_1a * kv * x10(t) - e_1a * kv * x11(t),
  diff(x12(t), t) = c_2a + c_1a * x7(t) - c_3a * x12(t),
  diff(x13(t), t) = a1 * x10(t) * x6(t) - c_6a * x13(t) - a3 * x2(t) * x13(t) + e_2a * x14(t),
  diff(x14(t), t) = a1 * x11(t) * x7(t) - e_2a * kv * x14(t),
  diff(x15(t), t) = c_2c + c_1c * x7(t) - c_3c * x15(t),
  y1(t) = x2(t),
  y2(t) = x10(t) + x13(t),
  y3(t) = x9(t),
  y4(t) = x1(t) + x2(t) + x3(t),
  y6(t) = x12(t),
  y5(t) = x7(t)
]);

IdentifiabilityODE(sigma, GetParameters(sigma)):
