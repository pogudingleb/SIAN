kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl":
# substitutions := [k4=1, b1=1, b2=1, b3=1, b4=1, b5=1, b6=1, b7=1, b8=1, b9=1, d1=1,
# d2=1, d3=1, d4=1, d5=1, d6=1, d7=1, d9=1, k1=1, k10=1, k11=1, k12=1, k13=1, k14=1, k15=1, k16=1, k17=1, k18=1,
# d10=1, d11=1, d12=1, d13=1, d14=1, d15=1, d16=1, d17=1, k14=1, k15=1, k16=1, k17=1, k18=1, k19=1, k2=1
# ];

sigma := [
    diff(s(t), t) = (1 - s(t)) * (k1 + k2 * c22(t) + k3) / ( (1 + b1 * c4(t)) * (1 + b2 * c13(t)) * (1 + b3 * c17(t)) * (1 + b4 * c22(t)) * (1 + b5 * c31(t))) - s(t) * (d1 * (1 + d3 * p(t)) + d2 * c31(t)),
    diff(p(t), t) = k4 / (1+b6*s(t)) - p(t) * ( ((1+d4*p(t))*(1+d5*c17(t))*(1+d6*c22(t))*(1+d7*cg(t))) / ((1+b7*c4(t))*(1+b8*c13(t))) ),
    diff(c4(t), t) = k11*ct2(t) + k12 - d10*c4(t),
    diff(c13(t), t) = k13*ct2(t) + k14 - d11*c13(t),
    diff(c17(t), t) = k15*ct17(t) + k16 - d12*c17(t),
    diff(c22(t), t) = k17*ct22(t) + k18 - d13*c22(t),
    diff(c31(t), t) = k19*ct2(t) + k20 - d14*c31(t),
    diff(cg(t), t) = k21*ct1(t) + k22 - d15*cg(t),
    diff(cTS(t), t) = k23*p(t) + k24 - d16*cTS(t),
    diff(cOX(t), t) = k25*cTS(t) + k26 - d17*cOX(t),
    diff(ct1(t), t) = k5*p(t) * (1+k9*cg(t)) / (4 +k9*cg(t)+k10*c4(t)) - d9*ct1(t)/(1+b9*cOX(t)),
    diff(ct2(t), t) = k6*p(t) * (1+k10*c4(t)) / (4+k9*cg(t)+k10*c4(t)),
    diff(ct17(t), t) = k7*p(t) / (4+k9*cg(t)+k10*c4(t)) - d9*ct17(t) / (1+b9*cOX(t)),
    diff(ct22(t), t) = k8*p(t) / (4+k9*cg(t)+k10*c4(t)) - d9*ct22(t) / (1+b9*cOX(t)),
    y1(t) = c4(t),
    y2(t) = s(t),
    y3(t) = c13(t),
    y4(t) = c17(t),
    y5(t) = c22(t),
    y6(t) = c31(t),
    y7(t) = cg(t),
    y8(t) = p(t),
    y9(t) = cg(t),
    y10(t) = cTS(t),
    y11(t) = cOX(t),
    y12(t) = ct1(t),
    y13(t) = ct2(t),
    y14(t) = ct17(t),
    y15(t) = ct22(t),
    ye(t) = 72 * (2*p(t) + 2*(1 - s(t))) / 4
];
# sigma := subs(substitutions, sigma);
output :=IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true, optimize_tr_basis=true):