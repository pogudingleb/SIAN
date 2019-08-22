# This example was constructed in 
# Saccomani, M., Audoly, S., Bellu, G., D'Angio, L.
# Examples of testing global identifiability of biological and biomedical models with the daisy software
# Section 6
read "../GlobalIdentifiability.mpl";

# The full example is the following
#
#sigma := table([
#  mu = [v_max, km, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20],
#  x_vars = [x1_, x2_, x3_, x4_, x5_, x6_, x7_, x8_, x9_, x10_, x11_, x12_, x13_, x14_, x15_, x16_, x17_, x18_, x19_, x20_],
#  y_vars = [y1_, y2_, y3_, y4_, y5_, y6_, y7_, y8_, y9_, y10_, y11_, y12_, y13_, y14_, y15_, y16_, y17_, y18_, y19_, y20_],
#  u_vars = [u],
#  x_eqs = [
#    x1_1 = -v_max * x1_0 / (km + x1_0) - p1 * x1_0 + u0,
#    x2_1 = p1 * x1_0 - p2 * x2_0,
#    x3_1 = p2 * x2_0 - p3 * x3_0,
#    x4_1 = p3 * x3_0 - p4 * x4_0,
#    x5_1 = p4 * x4_0 - p5 * x5_0,
#    x6_1 = p5 * x5_0 - p6 * x6_0,
#    x7_1 = p6 * x6_0 - p7 * x7_0,
#    x8_1 = p7 * x7_0 - p8 * x8_0,
#    x9_1 = p8 * x8_0 - p9 * x9_0,
#    x10_1 = p9 * x9_0 - p10 * x10_0,
#    x11_1 = p10 * x10_0 - p11 * x11_0,
#    x12_1 = p11 * x11_0 - p12 * x12_0,
#    x13_1 = p12 * x12_0 - p13 * x13_0,
#    x14_1 = p13 * x13_0 - p14 * x14_0,
#    x15_1 = p14 * x14_0 - p15 * x15_0,
#    x16_1 = p15 * x15_0 - p16 * x16_0,
#    x17_1 = p16 * x16_0 - p17 * x17_0,
#    x18_1 = p17 * x17_0 - p18 * x18_0,
#    x19_1 = p18 * x18_0 - p19 * x19_0,
#    x20_1 = p19 * x19_0 - p20 * x20_0
#  ],
#  y_eqs = [
#    y1_0 = x1_0,
#    y2_0 = x2_0,
#    y3_0 = x3_0,
#    y4_0 = x4_0,
#    y5_0 = x5_0,
#    y6_0 = x6_0,
#    y7_0 = x7_0,
#    y8_0 = x8_0,
#    y9_0 = x9_0,
#    y10_0 = x10_0,
#    y11_0 = x11_0,
#    y12_0 = x12_0,
#    y13_0 = x13_0,
#    y14_0 = x14_0,
#    y15_0 = x15_0,
#    y16_0 = x16_0,
#    y17_0 = x17_0,
#    y18_0 = x18_0,
#    y19_0 = x19_0,
#    y20_0 = x20_0
#  ]
#]);
#
#theta_g := GlobalIdentifiability(sigma, sigma[mu], 0.99);
#print(theta_g);

# However, one can see that the example is reducible and one can first check identifiability for a subsystem

sigma_reduced := table([
  mu = [v_max, km, p1],
  x_vars = [x1_],
  y_vars = [y1_],
  u_vars = [u],
  x_eqs = [
    x1_1 = -v_max * x1_0 / (km + x1_0) - p1 * x1_0 + u0
  ],
  y_eqs = [
    y1_0 = x1_0
  ]
]);

theta_g_reduced := GlobalIdentifiability(sigma_reduced, sigma_reduced[mu], 0.99);
print(theta_g_reduced);

# And the identifiability of the other parameters straightforwardly follows from
#
# p_{i + 1} = ( y_i * p_i - y_{i + 1}' ) / y_{i + 1}
# for 1 <= i <= 19
