kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl":
sigma := [
    diff(x1(t), t) = (-b*c*x1(t) - b*x1(t)*x4(t) + 1) / (c + x4(t)),
    diff(x2(t), t) = alpha*x1(t) - beta*x2(t),
    diff(x3(t), t) = gama*x2(t) - delta*x3(t),
    diff(x4(t), t) = (gama*sgm*x2(t)*x4(t) - delta*sgm*x3(t)*x4(t)) / (x3(t)),
    y(t) = x1(t)
];

out := IdentifiabilityODE(sigma, GetParameters(sigma), infolevel=-1);
quit;