read "../generate_tr_bases.mpl":

sigma := [
diff(s(t), t) = Lam - r0 * b * s(t) * i(t) / n - mu * s(t),
diff(e(t), t) = b * s(t) * i(t) / n - eps * e(t) - mu * e(t),
diff(i(t), t) = eps * e(t) - g * i(t) - mu * i(t),
diff(r(t), t) = g * i(t) - mu * r(t),
y1(t) = i(t) + r(t)
]: 

IdentifiabilityODE(sigma, GetParameters(sigma), "new_logs/new_SEIR", sub_transc=true):
quit;