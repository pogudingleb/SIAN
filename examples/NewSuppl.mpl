kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl":

substitutions := [
    CLtot = CLlin + (kiny * V1 * Rtot(t))/(Kss + 1/2 * ((Atot(t) / V1 - Rtot(t) - Kss) + g(t))), 
    c = 1/2 * ((Atot(t) / V1 - Rtot(t) - Kss) + g(t))
];

sigma := [
    diff(AD(t), t) = -ka * AD(t),
    diff(Atot(t), t) = ka * AD(t) - CLtot * c - Q * (c - AP(t)/V2),
    diff(AP(t), t) = Q * (c - AP(t)/V2),
    diff(Rtot(t), t) = ksyn - kdeg * Rtot(t) - (kint - kdeg) * Rtot * c / (Kss + c),
    diff(g(t), t) = g(t) * ((Atot(t) / V1 - Rtot(t) - Kss) * (a * AD(t) - CLtot * c - Q * (c - AP(t)/V2) - (ksyn - kdeg * Rtot(t) - (kint - kdeg) * Rtot * c / (Kss + c))) + 4 * Kss / V1 * (ka * AD(t) - CLtot * c - Q * (c - AP(t)/V2))),
    y1(t) = Rtot(t)
];
output :=IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true, optimize_tr_basis=true):