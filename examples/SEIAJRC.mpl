read "../IdentifiabilityODE.mpl";

sys := [
diff(S(t), t) = -b*S(t)*Ninv*A(t)*q - b*S(t)*Ninv*II(t) - b*S(t)*Ninv*J(t),
diff(A(t), t) = -E(t)*k*r + E(t)*k - A(t)*g1,
diff(E(t), t) = b*S(t)*Ninv*A(t)*q + b*S(t)*Ninv*II(t) + b*S(t)*Ninv*J(t) - E(t)*k,
diff(C(t), t) = alpha*II(t),
diff(J(t), t) = alpha*II(t) - g2*J(t),
diff(II(t), t) = -alpha*II(t) + E(t)*k*r - g1*II(t),
y2(t) = Ninv,
y(t) = C(t)
];

IdentifiabilityODE(sys, GetParameters(sys), sub_transc=true):
quit;