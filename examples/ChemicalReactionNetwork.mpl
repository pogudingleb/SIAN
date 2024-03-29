# Example 6.1 in the paper "Global Identifiability of Differential Models"
# Taken from 
# Conradi, C., Shiu, A.,
# Dynamics of post-translational modification systems: recent progress and future directions
# Eq. (1) in the published version (https://doi.org/10.1016/j.bpj.2017.11.3787), eq. (1.1) on arxiv (https://arxiv.org/pdf/1705.10913.pdf)
read "../IdentifiabilityODE.mpl";

sigma := [
  diff(x1(t), t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k4 * x6(t),
  diff(x2(t), t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k3 * x4(t),
  diff(x3(t), t) = k3 * x4(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
  diff(x4(t), t) = k1 * x1(t) * x2(t) - k2 * x4(t) - k3 * x4(t),
  diff(x5(t), t) = k4 * x6(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
  diff(x6(t), t) = -k4 * x6(t) - k5 * x6(t) + k6 * x3(t) * x5(t),
  y1(t) = x3(t),
  y2(t) = x2(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), weighted_ordering=true):
