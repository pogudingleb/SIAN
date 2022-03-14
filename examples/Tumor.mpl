# Example (with initial conditions assumed being unknown) from Section 3 of "Examples of testing global identifiability of biological and biomedical models with the DAISY software"
# by M.P. Saccomani, S. Audoly, G. Bellu, L. D'Angio

# read "../IdentifiabilityODE.mpl";
read "../generate_tr_bases.mpl";


sigma := [
  diff(x1(t), t) = -(k3 + k7) * x1(t) + k4 * x2(t),
  diff(x2(t), t) = k3 * x1(t) - (k4 + a * k5 + b * d * k5) * x2(t) + k6 * x3(t) + k6 * x4(t) + k5 * x2(t) * x3(t) + k5 * x2(t) * x4(t),
  diff(x3(t), t) = a * k5 * x2(t) - k6 * x3(t) - k5 * x2(t) * x3(t),
  diff(x4(t), t) = b * d * k5 * x2(t) - k6 * x4(t) - k5 * x2(t) * x4(t),
  diff(x5(t), t) = k7 * x1(t),
  y1(t) = x5(t)
  # y2(t) = a,
  # y3(t) = b,
  # y4(t) = d
];

# IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true):
IdentifiabilityODE(sigma, GetParameters(sigma), "new_logs/Tumor", sub_transc=true):
quit;
