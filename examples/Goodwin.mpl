# read "../IdentifiabilityODE.mpl";
read "../generate_tr_bases_with_dc.mpl":

sigma := [
diff(x2(t), t) = alpha*x1(t) - beta*x2(t),
diff(x4(t), t) = (gama*sgm*x2(t)*x4(t) - delta*sgm*x3(t)*x4(t)) / (x3(t)),
diff(x1(t), t) = (-b*c*x1(t) - b*x1(t)*x4(t) + 1) / (c + x4(t)),
diff(x3(t), t) = gama*x2(t) - delta*x3(t),
y(t) = x1(t)
];

# output :=IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true, optimize_tr_basis=true):
IdentifiabilityODE(sigma, GetParameters(sigma), "new_logs/Goodwin", sub_transc=true):