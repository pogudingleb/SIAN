read "../IdentifiabilityODE.mpl";

sigma := [
 diff(x1(t),t) = -p1*x1(t)*x3(t)-u(t),
 diff(x2(t),t) = p1*x1(t)*x3(t)-p2*x2(t),
 diff(x3(t),t) = p2*x2(t)+p3*x5(t)-p4*x3(t),
 diff(x4(t),t) = u(t)+p5*p4*x3(t),
 diff(x5(t),t) = (1-p5)*p4*x3(t)-p3*x5(t),
 y(t) = p2*x2(t)+p6*p3*x5(t)
]: 

IdentifiabilityODE(sigma, GetParameters(sigma)):