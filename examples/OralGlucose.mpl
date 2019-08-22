# Example (with initial conditions assumed being unknown) from Section III of "DAISY: an Efficient Tool to Test Global Identifiability. Some Case Studies"
# by G. Bellu, M.P. Saccomani

read "../GlobalIdentifiability.mpl";

sigma := table([
  mu = [p1, p2, p3, k, v],
  x_vars = [G, X, R, Ib, Gb],
  y_vars = [y1_,y2_, y3_],
  u_vars = [I_],
  x_eqs = [
    G1 = -(p1 + X0) * G0 + p1 * Gb0 + v * R0,
    X1 = -p2 * X0 + p3 * (I_0 - Ib0),
    R1 = k,
    Ib1 = 0,
    Gb1 = 0
    # Gb and Id are assumed to be known
  ],
  y_eqs = [
    y1_0 = G0,
    y2_0 = Ib0,
    y3_0 = Gb0
  ]
]);

# k and v are not locally identifiable
theta_g := GlobalIdentifiability(sigma, [p1, p2, p3, G0, X0, R0], 0.99);
print(theta_g);
