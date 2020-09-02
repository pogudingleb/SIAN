read "../IdentifiabilityODE.mpl":

with(CodeTools):

cases := [
    [
        [[x(t) = diff(y(t), t) + a * b, diff(z(t), t) = z(t) * b - c(t)], 2],
        {x_r1(t) = diff(y_r1(t), t) + a * b, diff(z_r1(t), t) = z_r1(t) * b - c_r1(t), x_r2(t) = diff(y_r2(t), t) + a * b, diff(z_r2(t), t) = z_r2(t) * b - c_r2(t)}
    ],
    [
        [[diff(x1(t), t) = a * x1(t) - b * x1(t) * x2(t), x2(t) = u(t)], 3],
        {
            diff(x1_r1(t), t) = a * x1_r1(t) - b * x1_r1(t) * x2_r1(t), x2_r1(t) = u_r1(t),
            diff(x1_r2(t), t) = a * x1_r2(t) - b * x1_r2(t) * x2_r2(t), x2_r2(t) = u_r2(t),
            diff(x1_r3(t), t) = a * x1_r3(t) - b * x1_r3(t) * x2_r3(t), x2_r3(t) = u_r3(t)
        }
    ]
]:

for i from 1 to nops(cases) do
    input := cases[i][1]:
    replica := GenerateReplica(input[1], input[2]):
    replica := {op(replica)}:
    if verify(replica, cases[i][2], set) then
        print("PASSED");
    else
        print("FAILED")
    end if;
end do:
