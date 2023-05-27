# The system is taken from equation (6)
# Wodarz, D., Nowak, M.
# Mathematical models of HIV pathogenesis and treatment 
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/bies.10196
# The first term in the w' equation is corrected based on the text after (6) and
# the original paper by the same authors
# Specific therapy regimes could lead to long-term immunological control of HIV
# https://doi.org/10.1073/pnas.96.25.14464
# Page 1
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(x(t), t) = lm - d * x(t) - beta * x(t) * v(t),
  diff(y(t), t) = beta * x(t) * v(t) - a * y(t),
  diff(v(t), t) = k * y(t) - u * v(t),
  diff(w(t), t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
  diff(z(t), t) = c * q * y(t) * w(t) - h * z(t),
  y1(t) = w(t),
  y2(t) = z(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), weighted_ordering=true):
