kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl":
sigma := [
diff(x(t), t) = -x(t)*a + z(t)*y(t) + a*y(t),
diff(w(t), t) = e*z(t) - w(t)*f + x(t)*y(t),
diff(z(t), t) = -c*z(t) - w(t)*d + x(t)*y(t),
diff(y(t), t) = b*x(t) + b*y(t) - x(t)*z(t),
g(t) = x(t)
];

out := IdentifiabilityODE(sigma, GetParameters(sigma), infolevel=-1);
quit;