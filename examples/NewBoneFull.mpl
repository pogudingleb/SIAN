kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
read "../IdentifiabilityODE.mpl";
substitutions := [
    nu14 = (H14 * (A2(t)*2/(A1(t) + delta21)) + k14) * A1(t),
    nu124 = k124 * (1 - phi124) + k412 * phi124 * H18124 * (A24(t)*A180/(A240 * A18(t)) ^ gamma24124),
    nu412 = k412 * (A4(t) / A40) * ((1-phi412) + phi412 * ((A17a(t) + A17b(t)) / (A17a0 + A17b0))),
    nu4u = (2 - H64) * ( 3 / 10 * GFR * A4(t) - H4u * H74u),
    nu514 = 464/1000 * k412 * (A4(t) / A40) * ((1-phi412) + phi412 * ((A17a(t) + A17b(t)) / (A17a0 + A17b0))),
    nu145 = 464/1000 * (k124 * (1 - phi124) + k412 * phi124 * H18124 * (A24(t)*A180/(A240 * A18(t)) ^ gamma24124)),
    nu5u = 88 / 100 * GFR * A5(t) - 88 / 100 * GFR * phi5u,
    nu35 = k35 * A3(t),
    nu58 = k58 * A5(t),
    nu85 = k85 * A8(t),
    k17aD = (k17D * (A17a0 + A17b0) + (phi2 * k7D * k17D * H2017minus / pi0c - H2817D) * phik17D * A17a0 * phi17a-(phi2 * k7D * k17D * H2017minus / pi0c - H2817D) * phik17D * (A17a0 + A17b0)) / (A17a0 + A17b0),
    k17Dprime = phi2 * k7D * k17D * H2017minus / pi0c - H2817D,
    k22s = k22D * A220*((A17a(t) + A17b(t)) / (A17a0 + A17b0))^gamma1722 * alpha722 * (A7(t) / Vvasc) / (delta722 * ((A17a(t) + A17b(t)) / (A17a0 + A17b0))^gamma1722 + (A7(t) / Vvasc))
    # T64plus = 1 + (exp(bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410)) - exp(-bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410))) / (exp(bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410)) + exp(-bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410))),
    # T64minus = 1 - (exp(bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410)) - exp(-bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410))) / (exp(bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410)) + exp(-bT64 * (A6(t) / Vvasc - deltaT64 * (A40 / A4(y))^gamma410)))
];

gammas := [
    gamma1719S = 2, # added manually
    gamma1920 = 2,
    gamma181920 = 2,
    gamma1722 = 2,
    gamma2021 = 2,
    gamma24124 = 2
];

sigma:= [
    # Gut
    diff(A1(t), t) = D1 * H21plus - nu14,
    diff(A2(t), t) = H62plus * (1 - A2(t)) - H62minus *A2(t),
    diff(A3(t), t) = D3 * F3 - k35 * A3(t),
    # Vasculature
    diff(A4(t), t) = nu124 - nu412 - nu4u + nu14,
    diff(A5(t), t) = nu514 - nu145 - nu5u - nu35 - nu58 + nu85, 
    diff(A6(t), t) = A9(t) - k6D * A6(t),
    diff(A7(t), t) = H4107minus * (A10(t) * 2) * A11(t) - k7D * A7(t),
    # Intracellular phosphate
    diff(A8(t), t) = nu58 - nu85,
    # Kidney
    diff(A9(t), t) = k9s * H79 * H59minus - k9D * A9(t),
    # Parathyroid gland
    diff(A10(t), t) = (1 - A10(t)) * alpha10 * (85/100 * T64minus + 15/100) - A10(t) * alpha10 * (85/100 * T64plus + 15/100),
    diff(A11(t), t) = k11 * H611minus - k11* A11(t),
    # Bone
    diff(A12(t), t) = nu412 - nu124 + k1312 *A13(t) - k1213 *A12(t),
    diff(A13(t), t) = -k1312*A13(t) + k1213 * A12(t),
    diff(A14(t), t) = nu514 - nu145 + k1514 * A15(t) - k1415 * A14(t),
    diff(A15(t), t) = k1415 * A14(t) - k1514 * A15(t),
    diff(A16(t), t) = k17D * (A17a0 + A17b0)/pi0c * H2016plus - k17D * (A17a0 + A17b0) * pi0c / (A160 * H2017plus) * A16(t),
    diff(A17a(t), t) = (K17D * (A17a0 + A17b0) * pi) / (A160 * H2017plus) * A16(t) * phi17a * (k17aD / k17Dprime) - k17aD * A17a(t),
    diff(A17b(t), t) = (K17D * (A17a0 + A17b0) * pi) / (A160 * H2017plus) * A16(t) * (1 - phi17a) * phik17D - k17Dprime * phik17D * A17b(t),
    diff(A18(t), t) = k18D * pi0c * A180 * H2418Splus - k18D * H2018Dplus * H2218Dminus * A18(t),
    diff(A19(t), t) = k1920 * A190 * ((A17a(t) + A17b(t)) / (A17a0 + A17b0))^gamma1719S - k1920 * (A19(t) / A190)^gamma1920 * (A18(t) / A180)^gamma181920 * A19(t),
    diff(A20(t), t) = k1920 * (A19(t) / A190)^gamma1920 * (A18(t) / A180)^gamma181920 * A19(t) - 1000 * k1920 * A20(t),
    diff(A21(t), t) = (( k21D * A210 + k2124 * A210 * A220 - k2124 * A240 ) / A200^gamma2021) * A20(t)^gamma2021 + k21D * A21(t) - k2124 * A21(t) * A22(t) + k2421 * A24(t),
    diff(A22(t), t) = k22s - k22D * A22(t) - k2124 * ( A23(t) * A22(t) + A21(t) * A22(t) ) - k2421 * (A24(t) + A25(t)),
    diff(A23(t), t) = k23D * A230 * (A16(t)/A160) * ( (A7(t) / Vvasc) + delta723 * A16(t) / A160)/(2*(A7(t) / Vvasc)) - k2124 * A23(t) * A22(t) + k2421 * A25(t) -  k23D*A23(t),
    diff(A24(t), t) = k2124 * A21(t) * A22(t) - k2421 * A24(t),
    diff(A25(t), t) = k2124 * A23(t) * A22(t) - k2421 * A25(t),
    # Osteoblast intracellular components
    diff(A26(t), t) = k26s - H726Dplus * A26(t),
    diff(A27(t), t) = k27s * H727splus - k26D * A27(t),
    diff(A28(t), t) = k28D * A26(t) * A27(t) + k28D * A28(t),
    y1(t) = A1(t),
    y2(t) = A2(t),
    y3(t) = A3(t),
    y4(t) = A4(t),
    y5(t) = A5(t),
    y6(t) = A6(t),
    y7(t) = A7(t),
    y8(t) = A8(t),
    y9(t) = A9(t),
    y10(t) = A10(t),
    y11(t) = A11(t),
    y12(t) = A12(t),
    y13(t) = A13(t),
    y14(t) = A14(t),
    y15(t) = A15(t),
    y16(t) = A16(t),
    y17a(t) = A17a(t),
    y17b(t) = A17b(t),
    y18(t) = A18(t),
    y19(t) = A19(t),
    y20(t) = A20(t),
    y21(t) = A21(t),
    y22(t) = A22(t),
    y23(t) = A23(t),
    y24(t) = A24(t),
    y25(t) = A25(t),
    y26(t) = A26(t),
    y27(t) = A27(t),
    y28(t) = A28(t)
];

sigma:= subs(gammas, subs(substitutions, sigma));
output :=IdentifiabilityODE(sigma, GetParameters(sigma), substitute_tr_basis=true, optimize_tr_basis=true):