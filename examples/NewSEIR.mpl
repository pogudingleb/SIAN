# Taken from
# N. Tuncer, T. Le
# "Structural and practical identifiability analysis of outbreak models"
# https://doi.org/10.1016/j.mbs.2018.02.004
# Equation (2.2) with prevalence observations
read "../IdentifiabilityODE.mpl":

sigma := [
diff(s(t), t) = Lam - r0 * b * s(t) * i(t) / n - mu * s(t),
diff(e(t), t) = b * s(t) * i(t) / n - eps * e(t) - mu * e(t),
diff(i(t), t) = eps * e(t) - g * i(t) - mu * i(t),
diff(r(t), t) = g * i(t) - mu * r(t),
y1(t) = i(t) + r(t)
]: 

IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true):
