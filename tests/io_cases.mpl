read "../IdentifiabilityODE.mpl":

with(CodeTools):

cases := [
    [
        [
            diff(x1(t), t) = a * x1(t) - b * x1(t) * x2(t),
            diff(x2(t), t) = -c * x2(t) + d * x1(t) * x2(t),
            y(t) = x1(t) + u(t)
        ],
        table([globally = {a, c, d, x1(0)}, locally_not_globally = {}, non_identifiable = {b, x2(0)}])
    ],
    [
        [
            diff(xA(t), t) = -k1 * xA(t),
            diff(xB(t), t) = k1 * xA(t) - k2 * xB(t),
            diff(xC(t), t) = k2 * xB(t),
            diff(eA(t), t) = 0,
            diff(eC(t), t) = 0,
            y1(t) = xC(t),
            y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
            y3(t) = eA(t),
            y4(t) = eC(t)
        ],
        table([globally = {xC(0), eA(0), eC(0)}, locally_not_globally = {eB, k1, k2, xA(0), xB(0)}, non_identifiable = {}])
    ],
    [
        [
            diff(S(t), t) = -b * S(t) * In(t) / N(t),
            diff(E(t), t) = b * S(t) * In(t) / N(t) - nu * E(t),
            diff(In(t), t) = nu * E(t) - a * In(t),
            diff(N(t), t) = 0,
            y1(t) = In(t),
            y2(t) = N(t)
        ],
        table([globally = {b, In(0), N(0)}, locally_not_globally = {a, nu, S(0), E(0)}, non_identifiable = {}])
    ],
    [
        [
            diff(x1(t), t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k4 * x6(t),
            diff(x2(t), t) = k1 * x1(t) * x2(t) + k2 * x4(t) + k3 * x4(t),
            diff(x3(t), t) = k3 * x4(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
            diff(x4(t), t) = k1 * x1(t) * x2(t) - k2 * x4(t) - k3 * x4(t),
            diff(x5(t), t) = k4 * x6(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
            diff(x6(t), t) = -k4 * x6(t) - k5 * x6(t) + k6 * x3(t) * x5(t),
            y1(t) = x3(t),
            y2(t) = x2(t)
        ],
        table([globally = {k1, k2, k3, k4, k5, k6, x1(0), x2(0), x3(0), x4(0), x5(0), x6(0)}, locally_not_globally = {}, non_identifiable = {}])
    ],
    [
        [
            diff(x(t), t) = 0,
            y1(t) = x(t),
            y2(t) = mu1 * x(t) + mu2
        ],
        table([globally = {x(0)}, locally_not_globally = {}, non_identifiable = {mu1, mu2}])
    ]
]:

num_passed := 0:
num_failed := 0:
num_passed_w := 0:
num_failed_w := 0:

for i from 1 to nops(cases) do
    sigma := cases[i][1]:
    correct_result := cases[i][2]:
    for method from 1 to 2 do
        result := IdentifiabilityODE(sigma, GetParameters(sigma), method = method, infolevel = 0, count_solutions = false):
        if verify(correct_result, result, table) then
            print("PASSED");
            num_passed := num_passed + 1:
        else
            print("FAILED");
            print(sigma, correct_result, result, method);
            num_failed := num_failed + 1:
        end if;
    end do;
    for method from 1 to 2 do
        result := IdentifiabilityODE(sigma, GetParameters(sigma), method = method, weighted_ordering=true, infolevel = 0, count_solutions = false):
        if verify(correct_result, result, table) then
            print("PASSED");
            num_passed_w := num_passed_w + 1:
        else
            print("FAILED");
            print(sigma, correct_result, result, method);
            num_failed_w := num_failed_w + 1:
        end if;
    end do;
end do:

printf("Passed: %a, failed %a \n", num_passed, num_failed);
printf("Passed with weights: %a, failed with weights %a \n", num_passed_w, num_failed_w);

