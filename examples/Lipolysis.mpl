# Taken from
# R. Munoz-Tamayo, L. Puillet, J.B. Daniel, D. Sauvant, O. Martin, M. Taghipoor, P. Blavy
# Review: To be or not to be an identifiable model. Is this a relevant question in animal science modelling?
# doi.org/10.1017/S1751731117002774
# System (1) in Supplementary Material 2, initial conditions are assumed to be unknown
# brought to the rational function form by introducing new state variable x5(t) = k1 e^(-k3 t)
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(x1(t), t) = -x1(t) * x5(t) / (k2 + x1(t)),
  diff(x2(t), t) = 2 * x1(t) * x5(t) / ((k2 + x1(t)) * 3) - k4 * x2(t),
  diff(x3(t), t) = k4*(x2(t))/2 - k4*x3(t),
  diff(x4(t), t) = x1(t) * x5(t) / (3 * (k2 + x1(t))) + k4 * (x2(t))/2 + k4 * x3(t),
  diff(x5(t), t) = -k3 * x5(t),
  y1(t) = x1(t),
  y2(t) = x2(t) + x3(t),
  y3(t) = x4(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma)):

