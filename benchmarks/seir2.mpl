# Ex 16 from  https://arxiv.org/abs/2006.14295
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl":
sigma := [
diff(s(t), t) = -b * s(t) * i(t),
diff(e(t), t) = b * s(t) * i(t) - eps * e(t),
diff(i(t), t) = eps * e(t) - (rho + mu) * i(t),
diff(r(t), t) = rho * i(t) - d0 * r(t),
y1(t) = i(t) + r(t)
]:


out := IdentifiabilityODE(sigma, GetParameters(sigma), infolevel=-1);
quit;