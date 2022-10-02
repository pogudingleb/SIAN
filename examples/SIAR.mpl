# read "../IdentifiabilityODE.mpl";
read "../generate_tr_bases.mpl";


sys := [
diff(E(t), t) = -alphaEIs*E(t) + betaIa*Ia(t)*S(t)*Ninv + Is(t)*S(t)*Ninv*betaIs + betaT*S(t)*Ninv*T(t) + S(t)*Ninv*betaH*H(t) - E(t)*alphaEIa,
diff(H(t), t) = -alphaHT*H(t) + Ia(t)*xi + Is(t)*alphaIsH - alphaHRd*H(t),
diff(Is(t), t) = alphaEIs*E(t) + Ia(t)*alphaIaIs - Is(t)*alphaIsT - Is(t)*alphaIsD - Is(t)*alphaIsH - Is(t)*alphaIsRu,
diff(D(t), t) = Is(t)*alphaIsD + T(t)*alphaTD,
diff(S(t), t) = -betaIa*Ia(t)*S(t)*Ninv - Is(t)*S(t)*Ninv*betaIs - betaT*S(t)*Ninv*T(t) - S(t)*Ninv*betaH*H(t),
diff(Ia(t), t) = -Ia(t)*alphaIaRu - Ia(t)*xi - Ia(t)*alphaIaIs + E(t)*alphaEIa,
diff(Rd(t), t) = alphaTRd*T(t) + alphaHRd*H(t),
diff(T(t), t) = alphaHT*H(t) + Is(t)*alphaIsT - alphaTRd*T(t) - T(t)*alphaTD,
y4(t) = D(t),
y2(t) = T(t),
y3(t) = Rd(t),
y5(t) = Ninv,
y1(t) = H(t)
];
# IdentifiabilityODE(sys, GetParameters(sys), substitute_tr_basis=true):
IdentifiabilityODE(sys, GetParameters(sys), "new_logs/SIAR", sub_transc=true):
quit;