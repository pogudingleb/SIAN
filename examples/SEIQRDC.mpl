read "../IdentifiabilityODE.mpl";

sys := [
diff(s(t), t) = mu * n - a * s(t) - b * n * s(t) * i(t) * n - mu * s(t),
diff(e(t), t) = b * s(t) * i(t) * n - mu * e(t) - g * e(t),
diff(i(t), t) = g * e(t) - dlt * i(t) - mu * i(t) * s(t), # originally has mu * i(t) * mu * s(t)
diff(q(t), t) = dlt * i(t) - l(t) * q(t) - k(t) * q(t) - mu * q(t),
diff(r(t), t) = l(t) * q(t) - mu * s(t),
diff(d0(t), t) = k(t) * q(t),
diff(c(t), t) = a * s(t) - mu * c(t) - tau0 * c(t),
y1(t) = c(t)
]: 

IdentifiabilityODE(sys, GetParameters(sys), substitute_tr_basis=true, optimize_tr_basis=true, infolevel=2):

