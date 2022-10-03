read "common.mpl":

with(CodeTools):


cases := [
    table(
        [
            input = [
                    diff(x1(t), t) = a*x1(t) + b^2 * x2(t), 
                    diff(x2(t), t) = c*x1(t) + d^2 * x2(t),
                    y1(t) = x1(t)
                ], 
            output = table(
                [
                    x1_0 = x1_0, x1_1 = x1_1, x1_2 = x1_2, x1_3 = x1_3, x1_4 = x1_4,
                    x2_0 = x2_0^2, x2_1 = x2_1^2, x2_2 = x2_2^2, x2_3 = x2_3^2,
                    z_aux = z_aux^2
                ]
            )
        ]
    ),
    table(
        [
            input = [
                    diff(x1(t), t) = a*x1(t), 
                    diff(x2(t), t) = c*x1(t) + d^2 * x2(t),
                    y1(t) = x1(t)
                ], 
            output = table(
                [
                    x1_2 = x1_2, x1_1 = x1_1, x1_0 = x1_0, z_aux = z_aux, a = a^2,
                    z_aux = z_aux
                ]
            )
        ]
    ),
    table(
        [
            input = [
                    diff(x1(t), t) = a*x1(t), 
                    diff(x2(t), t) = c*x1(t) + d^2 * x2(t),
                    y1(t) = x1(t) + x2(t)
                ], 
            output = table(
                [
                    x1_3 = x1_3, x1_2 = x1_2, x1_4 = x1_4, x1_1 = x1_1, x2_0 = x2_0, x2_4 = x2_4, x1_0 = x1_0, z_aux = z_aux, x2_1 = x2_1,
                    d = d^2 , x2_2 = x2_2, a = a^2, x2_3 = x2_3
                ]
            )
        ]
    ),
    table(
        [
            input = [
                    diff(x1(t), t) = x1(t), 
                    diff(x2(t), t) = x1(t) * x2(t),
                    y1(t) = x2(t)
                ], 
            output = table(
                [
                    x1_1 = x1_1^2, x2_0 = x2_0, x1_0 = x1_0^2, z_aux = z_aux, x2_1 = x2_1, x2_2 = x2_2
                ]
            )
        ]
    ),
    table(
        [
            input = [
                    diff(x1(t), t) = a, 
                    diff(x2(t), t) = a * b,
                    y1(t) = x1(t)
                ], 
            output = table(
                [
                     z_aux = z_aux, a = a^2
                ]
            )
        ]
    )
]:


num_passed := 0:
num_failed := 0:

for i from 1 to nops(cases) do
    in_ := cases[i][input]:
    out_ := GetWeights(in_, GetParameters(in_)):
    if verify([entries(out_, `pairs`)], [entries(cases[i][output], `pairs`)]) then
        print("PASSED");
        num_passed := num_passed + 1:
    else
        print("FAILED"):
        printf("Expected: %a,\nGot: %a\n\n", [entries(out_, `pairs`)], [entries(cases[i][output], `pairs`)]):
        num_failed := num_failed + 1:
    end if;
end do:

printf("Passed %d out of %d tests\n\n", num_passed, num_passed + num_failed):