# The system is taken from
# Wodarz, D., Nowak, M.
# Specific therapy regimes could lead to long-term immunological control of HIV
# https://doi.org/10.1073/pnas.96.25.14464
# Page 1
read "../IdentifiabilityODE.mpl":

sigma := [
  diff(x(t), t) = lm - d * x(t) - beta_ * x(t) * v(t),
  diff(y(t), t) = beta_ * x(t) * v(t) - a * y(t),
  diff(v(t), t) = k * y(t) - u * v(t),
  diff(w(t), t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
  diff(z(t), t) = c * q * y(t) * w(t) - h * z(t),
  y1(t) = w(t),
  y2(t) = z(t)
]:

# Et, x_theta_vars, u_hat, y_hat, all_subs := 
IdentifiabilityODE(sigma, GetParameters(sigma), "HIV2", sub_transc=true, infolevel=1):

# JacX := VectorCalculus[Jacobian](subs({op(u_hat), op(y_hat)}, Et), x_theta_vars = subs(all_subs, x_theta_vars)): 
# rrefJacX := LinearAlgebra[ReducedRowEchelonForm](JacX):
# pivots := {}:
# print(x_theta_vars);
# x_theta_vars := [a, b, d, h, u, w_0, z_0, z_6, z_5, z_4, z_3, z_2, z_1, y_6, y_5, y_4, y_3, y_2, y_1, y_0, x_6, x_5, x_4, x_3, x_2, x_1, x_0, w_7, w_6, w_5, w_4, w_3, w_2, w_1, v_5, v_4, v_3, v_2, v_1, v_0, k, c, lm, q, beta_]:
# for row_idx from 1 to nops(Et) do #nops(theta) do
#   row := rrefJacX[row_idx]:
#   pivot_idx := 1:
#   while row[pivot_idx]=0 and add(row)<>0 do
#     pivot_idx := pivot_idx + 1:
#   end do:
#   if pivot_idx <= numelems(row) then
#     pivots := {op(pivots), x_theta_vars[pivot_idx]}:
#   end if:
# end do:
# alg_indep := {op(x_theta_vars)} minus pivots:
# print(alg_indep):

# rhs_cols := []:
# idxs := []:
# for parameter in alg_indep do
#   idx := ListTools[Search](parameter, x_theta_vars):
#   idxs := [op(idxs), idx]:
#   rhs_cols := [op(rhs_cols), JacX[..,idx]]:
# end do:
# rhs_cols := Matrix(rhs_cols):
# print(idxs);
# A := LinearAlgebra[DeleteColumn](JacX, idxs):

# # test column space membership:
# for col_idx from 1 to numelems(alg_indep) do
#   solution := LinearAlgebra[LinearSolve](A, rhs_cols[.., col_idx]):
#   if add(solution)<>0 then
#     printf("\n%s %a %s\n", "Parameter", alg_indep[col_idx], "is in the span of columns" ):
#   end if;
# end do:
