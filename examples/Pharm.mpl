# Example 12 in the paper, taken from
# Demignot, S., D., D., 
# Effect of prosthetic sugar groups on the pharmacokinetics of glucose-oxidase
read("../GlobalIdentifiability.mpl");

sigma := table([
  mu = [a1, b1, b2, kc, ka, n],
  x_vars = [x1_, x2_, x3_, x4_],
  y_vars = [y1_],
  u_vars = [],
  x_eqs = [
    x1_1 = a1 * (x2_0 - x1_0) - (ka * n * x1_0) / (kc * ka + kc * x3_0 + ka * x1_0),
    x2_1 = a1 * (x1_0 - x2_0),
    x3_1 = b1 * (x4_0 - x3_0) - (kc * n * x3_0) / (kc * ka + kc * x3_0 + ka * x1_0),
    x4_1 = b2 * (x3_0 - x4_0)
  ],
  y_eqs = [
    y1_0 = x1_0
  ]
]);

p := GlobalIdentifiability(sigma, [op(sigma[mu]), x1_0, x2_0, x3_0, x4_0], 0.99);

print(p);
