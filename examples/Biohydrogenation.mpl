# Taken from
# R. Munoz-Tamayo, L. Puillet, J.B. Daniel, D. Sauvant, O. Martin, M. Taghipoor, P. Blavy
# Review: To be or not to be an identifiable model. Is this a relevant question in animal science modelling?
# doi.org/10.1017/S1751731117002774
# System (3) in Supplementary Material 2, initail conditions are assumed to be unknown
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(x4(t), t) = - k5 * x4(t) / (k6 + x4(t)),
  diff(x5(t), t) = k5 * x4(t) / (k6 + x4(t)) - k7 * x5(t)/(k8 + x5(t) + x6(t)),
  diff(x6(t), t) = k7 * x5(t) / (k8 + x5(t) + x6(t)) - k9 * x6(t) * (k10 - x6(t)) / k10,
  diff(x7(t), t) = k9 * x6(t) * (k10 - x6(t)) / k10,
  y1(t) = x4(t),
  y2(t) = x5(t)
]:

IdentifiabilityODE(sigma, GetParameters(sigma), p = 0.999):
