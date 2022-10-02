# Example (with initial conditions assumed being unknown) from Section 3 of "Examples of testing global identifiability of biological and biomedical models with the DAISY software"
# by M.P. Saccomani, S. Audoly, G. Bellu, L. D'Angio

# read "../IdentifiabilityODE.mpl";
read "../generate_tr_bases.mpl";


sigma := [
  diff(x1(t), t) = -(k3 + k7) * x1(t) + k4 * x2(t),
  diff(x2(t), t) = k3 * x1(t) - (k4 + a * k5 + b * d * k5) * x2(t) + k6 * x3(t) + k6 * x4(t) + k5 * x2(t) * x3(t) + k5 * x2(t) * x4(t),
  diff(x3(t), t) = a * k5 * x2(t) - k6 * x3(t) - k5 * x2(t) * x3(t),
  diff(x4(t), t) = b * d * k5 * x2(t) - k6 * x4(t) - k5 * x2(t) * x4(t),
  diff(x5(t), t) = k7 * x1(t),
  y1(t) = x5(t)
  # y2(t) = a,
  # y3(t) = b,
  # y4(t) = d
];

<<<<<<< HEAD
# IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true):
IdentifiabilityODE(sigma, GetParameters(sigma), "new_logs/Tumor", sub_transc=true):
quit;
||||||| e3ff1e1
IdentifiabilityODE(sigma, GetParameters(sigma));
=======
sigma:=[
  diff(q1(t), t) = k4*q3(t)-(k3+k7)*q1(t)+k2*q12(t)-k1*q1(t)*(PP*V1-q12(t)),
  diff(q12(t), t) = k1*q1(t)*(PP*V1-q12(t))-k2*q12(t),
  diff(q3(t), t) = k3*q1(t)-k4*q3(t)-k1*q3(t)*(PE*V3-q34(t))+k2*q34(t)-k5*q3(t)*(R*V3-q35(t))+k6*q35(t)-5*k5*V36/V3*q3(t)*(S*V36-q36(t))+k6*q36(t),
  diff(q34(t), t) = k1*q3(t)*(PE*V3-q34(t))-k2*q34(t),
  diff(q35(t), t) = k5*q3(t)*(R*V3-q35(t))-k6*q35(t),
  diff(q36(t), t) = 5*k5*V36/V3*q3(t)*(S*V36-q36(t))-k6*q36(t),
  diff(x5(t), t) = k7*q1(t),
  y1(t) = x5(t)
]:

IdentifiabilityODE(sigma, GetParameters(sigma), infolevel=3);

# dq1/dt = k4 * q3 - (k3 + k7) * q1 + k2 * q12 - k1 * q1 * (PP * V1 - q12);
# dq12/dt = k1 * q1 * (PP * V1 - q12) - k2 * q12;
# dq3/dt = k3 * q1 - k4 * q3 - k1 * q3 * (PE * V3 - q34) + k2 * q34 - k5 * q3 * (R*V3 - q35) + k6 * q35 - k5 * (5 * V36 / V3) * q3 * (S * V36 - q36) + k6 * q36;
# dq34/dt = k1 * q3 * (PE * V3 - q34) - k2 * q34;
# dq35/dt = k5 * q3 * (R*V3 - q35) - k6 * q35;
# dq36/dt = k5 * (5 * V36 / V3) * q3 * (S * V36 - q36) - k6 * q36;
# dx5/dt = k7 * q1;
# y1 = x5;
>>>>>>> master
