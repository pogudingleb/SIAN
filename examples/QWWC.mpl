# read "../IdentifiabilityODE.mpl";
#read "../generate_tr_bases_with_dc.mpl":
read "../generate_tr_bases.mpl":

sigma := [
diff(x(t), t) = -x(t)*a + z(t)*y(t) + a*y(t),
diff(w(t), t) = e*z(t) - w(t)*f + x(t)*y(t),
diff(z(t), t) = -c*z(t) - w(t)*d + x(t)*y(t),
diff(y(t), t) = b*x(t) + b*y(t) - x(t)*z(t),
g(t) = x(t)
];
# IdentifiabilityODE(sys, GetParameters(sys), substitute_tr_basis=true, optimize_tr_basis=true, infolevel=2):
IdentifiabilityODE(sigma, GetParameters(sigma), "new_logs/QWWC", sub_transc=true):
quit;