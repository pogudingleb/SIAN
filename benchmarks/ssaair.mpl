# Ex 37 from  https://arxiv.org/abs/2006.14295
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl":
sigma := [
  diff(Sd(t), t) = -eps_s * b_a * (An(t) + eps_a * Ad(t)) * Sd(t) - h1 * Sd(t) + h2 * Sn(t) - eps_s * b_i * Sd(t) * In(t),
  diff(Sn(t), t) = -b_i * Sn(t) * In(t) - b_a * (An(t) + eps_a*Ad(t)) * Sn(t) + h1 * Sd(t) - h2 * Sn(t),
  diff(Ad(t), t) = eps_s * b_i * Sd(t) * In(t) + eps_s * b_a * (An(t) + eps_a * Ad(t)) * Sn(t) + h2 * An(t) - g_ai * Ad(t) - h1 * Ad(t),
  diff(An(t), t) = b_i * Sn(t) * In(t) + b_a * (An(t) + eps_a * Ad(t)) * Sn(t) + h1 * Ad(t) - g_ai * An(t) - h2 * An(t),
  diff(In(t), t) = f * g_ai * (Ad(t) + An(t)) - dlt * In(t) - g_ir * In(t),
  diff(R(t), t) = (1-f) * g_ai * (Ad(t) + An(t)) + g_ir * In(t),
  y1(t) = Sd(t),
  y2(t) = In(t)
]:

out := IdentifiabilityODE(sigma, GetParameters(sigma), infolevel=-1);
quit;