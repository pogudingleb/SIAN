#===============================================================================
GlobalIdentifiability := proc(sigma, theta_l, p := 0.99, method := 1, num_nodes := 5) 
#===============================================================================
 local i, j, k, n, m, s, all_params, all_vars, eqs, Q, X, Y, poly, d0, D1, 
        sample, all_subs,alpha, beta, Et, x_theta_vars, prolongation_possible, 
        eqs_i, JacX, vars, vars_to_add, ord_var, var_index, deg_variety, D2, 
        y_hat, u_hat, theta_hat, Et_hat, Q_hat, theta_g, gb, v, X_eq, Y_eq, poly_d, 
        separant, leader,vars_local:

  #----------------------------------------------
  # 1. Construct the maximal system.
  #----------------------------------------------

  # (a) ---------------
  n := nops(sigma[x_vars]):
  m := nops(sigma[y_vars]):
  s := nops(sigma[mu]) + n:
  all_params := [op(sigma[mu]), op(map(x -> MakeDerivative(x, 0), sigma[x_vars] ))]:
  all_vars := [ op(sigma[x_vars]), op(sigma[y_vars]), op(sigma[u_vars]) ]:
  eqs := [op(sigma[x_eqs]), op(sigma[y_eqs])]:
  Q := foldl( (f, g) -> lcm(f, g), op( map(f -> denom(rhs(f)), eqs) )):
  

  # (b,c) ---------------
  X := []:
  X_eq := []:
  for i from 1 to n do
    X := [op(X), []]:
    poly := numer(lhs(sigma[x_eqs][i]) - rhs(sigma[x_eqs][i])):
    for j from 0 to s do
      poly_d := Differentiate(poly, all_vars, s, j):
      leader := MakeDerivative(sigma[x_vars][i], j + 1):
      separant := diff(poly_d, leader):
      X[i] := [op(X[i]), poly_d]:
      X_eq := [op(X_eq), leader = -(poly_d - separant * leader) / separant]:
    end do:
  end do:
  
  # (d,e) ---------------
  Y := []:
  Y_eq := []:
  for i from 1 to m do
    Y := [op(Y), []]:
    poly := numer(lhs(sigma[y_eqs][i]) - rhs(sigma[y_eqs][i])):
    for j from 0 to s do
      poly_d := Differentiate(poly, all_vars, s, j):
      leader := MakeDerivative(sigma[y_vars][i], j):
      separant := diff(poly_d, leader):
      Y[i] := [op(Y[i]), poly_d]:
      Y_eq := [op(Y_eq), leader = -(poly_d - separant * leader) / separant]:
    end do:
  end do:


  #----------------------------------------------
  # 2. Truncate.
  #----------------------------------------------

  # (a) ---------------
  d0 := max(op( map(f -> degree( simplify(Q * rhs(f)) ), eqs) ), degree(Q)):

  # (b) ---------------
  D1 := floor( 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p) ):
  print("Bound D_1  ", D1);

  # (c, d) ---------------
  sample := SamplePoint(D1, sigma, X_eq, Y_eq):
  all_subs := sample[4]:
  while subs(all_subs, Q) = 0 do
    sample := SamplePoint(D1, sigma, X_eq, Y_eq):
    all_subs := sample[4]:
  end do:
 
  # (e) ------------------
  alpha := [seq(1, i = 1..n)]:
  beta := [seq(0, i = 1..m)]:
  Et := [];
  x_theta_vars := all_params:
  prolongation_possible := [seq(1, i = 1..m)]:

  # (f) ------------------
  while add(prolongation_possible) > 0 do
    for i from 1 to m do
      if prolongation_possible[i] = 1 then
        eqs_i := [op(Et), Y[i][beta[i] + 1]]:
        JacX := subs(all_subs, VectorCalculus[Jacobian](eqs_i, x_theta_vars = subs(all_subs, x_theta_vars)));
        if LinearAlgebra[Rank](JacX) = nops(eqs_i) then
          Et := [op(Et), Y[i][beta[i] + 1]]:
          beta[i] := beta[i] + 1:
          for j from 1 to s + 1 do
            vars := {};
            for poly in [op(Et), seq(Y[k][beta[k] + 1], k=1..m)] do
              vars := vars union { op(GetVars(poly, sigma[x_vars], s + 1)) }:
            end do:
            vars_to_add := { op(remove(v -> evalb(v in x_theta_vars), vars)) };
            for v in vars_to_add do
              x_theta_vars := [op(x_theta_vars), v];
              ord_var := GetOrderVar(v, all_vars, s + 1);
              var_index := ListTools[Search](ord_var[2], sigma[x_vars]):
              poly := X[ var_index ][ ord_var[1] ]:
              Et := [op(Et), poly]:
              alpha[ var_index ] := max(alpha[ var_index ], ord_var[1] + 1):
            end do:
          end do:
        else
          prolongation_possible[i] := 0;
        end if:
      end if: 
    end do:
    print("Beta now", beta);
  end do:

  # (g) --------------
  for i from 1 to m do
    for j from beta[i] + 1 to nops(Y[i]) do
      to_add := true:
      for v in GetVars(Y[i][j], sigma[x_vars], s + 1) do
        if not (v in vars) then
          to_add := false:
        end if:
      end do:
      if to_add = true then
        beta[i] := beta[i] + 1:
        Et := [op(Et), Y[i][j]]:
      end if:
    end do:
  end do:
 
  print("Beta ", beta);
  print("Alpha ", alpha);
  deg_variety := foldl(`*`, op( map(e -> degree(e), Et) )):
  
  #----------------------------------------------
  # 3. Randomize.
  #----------------------------------------------

  # (a) ------------
  D2 := floor( 6 * nops(theta_l) * deg_variety * (1 + 2 * d0 * max(op(beta))) / (1 - p) ):
  print("Bound D_2 ", D2):

  # (b, c) ---------
  sample := SamplePoint(D2, sigma, X_eq, Y_eq):
  while subs(sample[4], Q) = 0 do
    sample := SamplePoint(D2, sigma, X_eq, Y_eq):
  end do:    
  y_hat := sample[1]:
  u_hat := sample[2]:
  theta_hat := sample[3]:

  # (d) ------------
  Et_hat := map(e -> subs([op(y_hat), op(u_hat)], e), Et):
  vars := { op(sigma[mu]) };
  for poly in Et_hat do
    vars := vars union { op(GetVars(poly, sigma[x_vars], s + 1)) }:
  end do:
  print("We finally have ", nops(Et_hat), "equations in ", nops(vars), "variables");
  Q_hat := subs(u_hat, Q):

  #----------------------------------------------
  # 4. Determine.
  #----------------------------------------------

  theta_g := []:
  if method = 1 then
    at_node := proc(var, args_node)
      local gb_loc;
      gb_loc := Groebner[Basis](op(args_node)):
      #print("Groebner basis for ", var, " is ", gb_loc);
      gb_loc;
    end proc:

    if nops(theta_l) > 1 then
      Grid[Setup]("local", numnodes = num_nodes):
      Grid[Set](at_node):
      gb := Grid[Seq](
        at_node(theta_l[i], [
          [op(Et_hat), z * Q_hat - 1, (theta_l[i] - subs(theta_hat, theta_l[i])) * w - 1],
          tdeg(op(vars), z, w)
        ]),
        i = 1..nops(theta_l)
      ):
    else
      # This is needed because of a bug in Grid[Seq]
      gb := [ at_node(theta_l[1], [
        [op(Et_hat), z * Q_hat - 1, (theta_l[1] - subs(theta_hat, theta_l[1])) * w - 1],
        tdeg(op(vars), z, w)
      ]) ]:
    end if:

    for i from 1 to nops(theta_l) do
      if gb[i] = [1] then
        theta_g := [op(theta_g), theta_l[i]]:
      end if:
    end do:     
  elif method = 2 then
    gb := Groebner[Basis]([op(Et_hat), z * Q_hat - 1], tdeg(op(vars), z));
    for i from 1 to nops(theta_l) do
      if Groebner[NormalForm](theta_l[i], gb, tdeg(op(vars), z)) = subs(theta_hat, theta_l[i]) then
        theta_g := [ op(theta_g), theta_l[i] ]:
      end if:
    end do:
  else
    print("No such method"):
  end if:
  theta_g;
end proc:


#===============================================================================
MakeDerivative := proc(var_name, der_order):
#===============================================================================
  cat(var_name, der_order):
end proc:


#===============================================================================
DifferentiateOnce := proc(diff_poly, var_list, max_ord) 
#===============================================================================
  local result, i, j:
  result := 0:
  for i from 1 to nops(var_list) do
    for j from 0 to max_ord do
      result := result + diff(diff_poly, MakeDerivative(var_list[i], j)) * MakeDerivative(var_list[i], j + 1):
    end do:
  end do:
  simplify(result):
end proc:


#===============================================================================
Differentiate := proc(diff_poly, var_list, max_ords, ord := 1) 
#===============================================================================
  local result, i;
  result := diff_poly:
  for i from 1 to ord do
    result := DifferentiateOnce(result, var_list, max_ords):
  end do:
  result:
end proc:


#===============================================================================
GetVars := proc(diff_poly, var_list, max_ord) 
#===============================================================================
  local all_vars, result;
  all_vars := map( v -> op(map( i -> MakeDerivative(v, i), [`$`(0..max_ord)] )), var_list):
  result := select(v -> evalb(diff(diff_poly, v) <> 0), all_vars):
  result:
end proc:


#===============================================================================
GetOrderVar := proc(diff_var, var_list, max_ord)
#===============================================================================
  local result, v, h;
  result := -1:
  for v in var_list do
    for h from 0 to max_ord do
      if diff(diff_var, MakeDerivative(v, h)) <> 0 then
        result := [h, v];
      end if;
    end do:
  end do:
  result;
end proc:


#===============================================================================
SamplePoint := proc(bound, sigma, X_eq, Y_eq)
#===============================================================================
  local n, m, s, all_params, all_vars, theta_hat, u_variables, 
        u_hat, x_hat, y_hat, eq, eq_num, eq_denom, 
        v, poly, i, j, all_subs, roll;
  n := nops(sigma[x_vars]):
  m := nops(sigma[y_vars]):
  s := nops(sigma[mu]) + n:
  all_params := [op(sigma[mu]), op(map(x -> MakeDerivative(x, 0), sigma[x_vars] ))]:
  all_vars := [ op(sigma[x_vars]), op(sigma[y_vars]), op(sigma[u_vars]) ]:

  roll := rand(0 .. bound):
  theta_hat := map(p -> p = roll(), all_params): 
  u_variables := [];
  for i from 1 to nops(sigma[u_vars]) do
    u_variables := [ op(u_variables), seq(MakeDerivative(sigma[u_vars][i], j), j = 0..s) ]:
  end do:
  u_hat := map(p -> p = roll(), u_variables) :   

  x_hat := X_eq;
  y_hat := Y_eq;
  all_subs := [op(theta_hat), op(u_hat)]:
  for i from 1 to s + 1 do
    x_hat := map(e -> lhs(e) = subs(all_subs, rhs(e)), x_hat):
    y_hat := map(e -> lhs(e) = subs(all_subs, rhs(e)), y_hat):
    all_subs := [ op(all_subs), op(select(e -> type(rhs(e), numeric), [op(x_hat), op(y_hat)])) ]:
  end do:

  [y_hat, u_hat, theta_hat, all_subs];
end proc:
