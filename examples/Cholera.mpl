# Example 6.2 in the paper "Global Identifiability of Differential Models", taken from
# Lee, E. C., Kelly, M. R., Ochocki, B. M., Akinwumi, S. M., Hamre, K. E., Tien, J. H., Eisenberg, M. C.,
# Model distinguishability and inference robustness in mechanisms of cholera transmission and loss of immunity
# Eq. (3)
read "../IdentifiabilityODE.mpl";

sigma := [
  diff(s(t), t) = mu - bi * s(t) * i(t) - bw * s(t) * w(t) - mu * s(t) + al * r(t),
  diff(i(t), t) = bw * s(t) * w(t) + bi * s(t) * i(t) - g * i(t) - mu * i(t),
  diff(w(t), t) = dz * (i(t) - w(t)),
  diff(r(t), t) = g * i(t) - mu * r(t) - al * r(t),
  y1(t) = k * i(t),
  y2(t) = i(t) + r(t) + s(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):
