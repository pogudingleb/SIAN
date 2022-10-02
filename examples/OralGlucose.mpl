# Example (with initial conditions assumed being unknown) from Section III of "DAISY: an Efficient Tool to Test Global Identifiability. Some Case Studies"
# by G. Bellu, M.P. Saccomani

# read "../IdentifiabilityODE.mpl";
read "../generate_tr_bases.mpl":

sigma := [
  diff(G(t), t) = -(p1 + X(t)) * G(t) + p1 * Gb + v * R(t),
  diff(X(t), t) = -p2 * X(t) + p3 * (u(t) - Ib),
  diff(R(t), t) = k,
  y1(t) = G(t),
  y2(t) = Ib,
  y3(t) = Gb
];

# IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true):
IdentifiabilityODE(sigma, GetParameters(sigma), "new_logs/OralGlucose", sub_transc=true):
