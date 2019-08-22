# Example 10 in the paper, taken from
# Lee, E. C., Kelly, M. R., Ochocki, B. M., Akinwumi, S. M., Hamre, K. E., Tien, J. H., Eisenberg, M. C.,
# Model distinguishability and inference robustness in mechanisms of cholera transmission and loss of immunity
# Eq. (3)
read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [mu, bi, bw, al, k, g, dz],
  x_vars = [s, i, w, r],
  y_vars = [y1_, y2_],
  u_vars = [],
  x_eqs = [
    s1 = mu - bi * s0 * i0 - bw * s0 * w0 - mu * s0 + al * r0,
    i1 = bw * s0 * w0 + bi * s0 * i0 - g * i0 - mu * i0,
    w1 = dz * (i0 - w0),
    r1 = g * i0 - mu * r0 - al * r0
  ],
  y_eqs = [
    y1_0 = k * i0,
    y2_0 = i0 + r0 + s0
  ]
]);

theta_g := GlobalIdentifiability(sigma, [op(sigma[mu]), s0, i0, w0, r0], 0.99);
print(theta_g);
