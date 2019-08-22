#===============================================================================
IdentifiabilityODE := proc(system_ODEs, params_to_assess, {p := 0.99, infolevel := 1, method := 1, num_nodes := 6}) 
#===============================================================================
 local i, j, k, n, m, s, all_params, all_vars, eqs, Q, X, Y, poly, d0, D1, 
        sample, all_subs,alpha, beta, Et, x_theta_vars, prolongation_possible, 
        eqs_i, JacX, vars, vars_to_add, ord_var, var_index, deg_variety, D2, 
        y_hat, u_hat, theta_hat, Et_hat, Q_hat, theta_l, theta_g, gb, v, X_eq, Y_eq, 
        poly_d, separant, leader,vars_local, x_functions, y_functions, u_functions,
        all_symbols_rhs, mu, x_vars, y_vars, u_vars, theta, subst_first_order,
        subst_zero_order, x_eqs, y_eqs, param, other_params, to_add, at_node,
        prime, max_rank, R, tr, e, p_local:

  #----------------------------------------------
  # 0. Extract inputs, outputs, states, and parameters from the system
  #----------------------------------------------

  if infolevel > 0 then
    PrintHeader("0. Extracting states, inputs, outputs, and parameters from the system"):
  end if:

  x_functions := map(f -> int(f, t), select( f -> type(int(f, t), function(name)), map(lhs, system_ODEs) )):
  y_functions := select( f -> not type(int(f, t), function(name)), map(lhs, system_ODEs) ):
  all_symbols_rhs := foldl(`union`, op( map(e -> indets(rhs(e)), system_ODEs) )) minus {t}:
  u_functions := select( f -> type(f, function), convert(all_symbols_rhs minus {op(x_functions), op(y_functions)}, list)):
  mu := convert(all_symbols_rhs minus {op(x_functions), op(y_functions), op(u_functions)}, list):

  x_vars := map(FunctionToVariable, x_functions):
  y_vars := map(FunctionToVariable, y_functions):
  u_vars := map(FunctionToVariable, u_functions):
  theta := map(ParamToInner, params_to_assess):
  subst_first_order := {seq(diff(x_functions[i], t) = MakeDerivative(x_vars[i], 1), i = 1 .. nops(x_vars))}:
  subst_zero_order := {
    seq(x_functions[i] = MakeDerivative(x_vars[i], 0), i = 1 .. nops(x_vars)),
    seq(y_functions[i] = MakeDerivative(y_vars[i], 0), i = 1 .. nops(y_vars)),
    seq(u_functions[i] = MakeDerivative(u_vars[i], 0), i = 1 .. nops(u_vars))
  }:
  x_eqs := subs(subst_zero_order, subs(subst_first_order, select(e -> type(int(lhs(e), t), function(name)), system_ODEs))):
  y_eqs := subs(subst_zero_order, select(e -> not type(int(lhs(e), t), function(name)), system_ODEs)):

  # taking into account that fact that Groebner[Basis] is Monte Carlo with probability of error 
  # at most 10^(-18) (for Maple 2017)
  p_local := p + nops(params_to_assess) * 10^(-18):
  if p_local >= 1 then
    printf("The probability of success cannot exceed 1 - #params_to_assess 10^{-18}. We reset it to 0.99");
    p_local := 0.99:
  end if:

  if infolevel > 0 then
    printf("\n=== Input info ===\n"):
    printf("%s %a\n", `State variables:         `, x_functions):
    printf("%s %a\n", `Output variables:        `, y_functions):
    printf("%s %a\n", `Input variables:         `, u_functions):
    printf("%s %a\n", `Parameters in equations: `, mu):
    printf("===================\n\n"):
  end if:

  #----------------------------------------------
  # 1. Construct the maximal system.
  #----------------------------------------------

  if infolevel > 0 then
    PrintHeader("1. Constructing the maximal polynomial system"):
  end if:

  # (a) ---------------
  n := nops(x_vars):
  m := nops(y_vars):
  s := nops(mu) + n:
  all_params := [op(mu), op(map(x -> MakeDerivative(x, 0), x_vars ))]:
  all_vars := [ op(x_vars), op(y_vars), op(u_vars) ]:
  eqs := [op(x_eqs), op(y_eqs)]:
  Q := foldl( (f, g) -> lcm(f, g), op( map(f -> denom(rhs(f)), eqs) )):
  

  # (b,c) ---------------
  X := []:
  X_eq := []:
  for i from 1 to n do
    X := [op(X), []]:
    poly := numer(lhs(x_eqs[i]) - rhs(x_eqs[i])):
    for j from 0 to s do
      poly_d := Differentiate(poly, all_vars, s, j):
      leader := MakeDerivative(x_vars[i], j + 1):
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
    poly := numer(lhs(y_eqs[i]) - rhs(y_eqs[i])):
    for j from 0 to s do
      poly_d := Differentiate(poly, all_vars, s, j):
      leader := MakeDerivative(y_vars[i], j):
      separant := diff(poly_d, leader):
      Y[i] := [op(Y[i]), poly_d]:
      Y_eq := [op(Y_eq), leader = -(poly_d - separant * leader) / separant]:
    end do:
  end do:


  #----------------------------------------------
  # 2. Truncate.
  #----------------------------------------------

  if infolevel > 0 then
    PrintHeader("2. Truncating the polynomial system based on the Jacobian condition"):
  end if:

  # (a) ---------------
  d0 := max(op( map(f -> degree( simplify(Q * rhs(f)) ), eqs) ), degree(Q)):

  # (b) ---------------
  # extra factor nops(theta) + 1 compared to the formula in the paper is to
  # provide probability gaurantee to the local identifiability test
  D1 := floor( (nops(theta) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p_local) ):
  # prime := nextprime(D1):
  if infolevel > 1 then
    printf("%s %a\n", `Bound D_1 for testing the rank of the Jacobian probabilistically: `, D1);
  end if:

  # (c, d) ---------------
  sample := SamplePoint(D1, x_vars, y_vars, u_vars, mu, X_eq, Y_eq):
  all_subs := sample[4]:
  while subs(all_subs, Q) = 0 do
    sample := SamplePoint(D1, x_vars, y_vars, u_vars, mu, X_eq, Y_eq):
    all_subs := sample[4]:
  end do:
  u_hat := sample[2]:
  y_hat := sample[1]:
 
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
        JacX := VectorCalculus[Jacobian](subs({op(u_hat), op(y_hat)}, eqs_i), x_theta_vars = subs(all_subs, x_theta_vars)):
        if LinearAlgebra[Rank](JacX) = nops(eqs_i) then
          Et := [op(Et), Y[i][beta[i] + 1]]:
          beta[i] := beta[i] + 1:
          for j from 1 to s + 1 do
            vars := {};
            for poly in [op(Et), seq(Y[k][beta[k] + 1], k=1..m)] do
              vars := vars union { op(GetVars(poly, x_vars, s + 1)) }:
            end do:
            vars_to_add := { op(remove(v -> evalb(v in x_theta_vars), vars)) };
            for v in vars_to_add do
              x_theta_vars := [op(x_theta_vars), v];
              ord_var := GetOrderVar(v, all_vars, s + 1);
              var_index := ListTools[Search](ord_var[2], x_vars):
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
  end do:
  # is used for assessing local identifiabilty
  max_rank := nops(Et):

  # (g) --------------
  for i from 1 to m do
    for j from beta[i] + 1 to nops(Y[i]) do
      to_add := true:
      for v in GetVars(Y[i][j], x_vars, s + 1) do
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
 
  if infolevel > 1 then
    printf("%s %a\n", `Orders of prolongations of the outputs (beta) = `, beta):
    printf("%s %a\n", `Orders of prolongations of the state variables (alpha) = `, alpha):
  end if:
 
  ##############################

  if infolevel > 0 then
    PrintHeader("3. Assessing local identifiability"):
  end if:
 
  theta_l := []:
  for param in theta do
    other_params := subs(param = NULL, x_theta_vars):
    JacX := VectorCalculus[Jacobian]( 
        subs( { op(u_hat), param = subs(all_subs, param), op(y_hat) }, Et), 
        other_params = subs(all_subs, other_params)
    ):
    if LinearAlgebra[Rank](JacX) <> max_rank then
      theta_l := [op(theta_l), param]:
    end if:
  end do:
 
  if infolevel > 1 then
    printf("%s %a\n", `Locally identifiable paramters: `, map(ParamToOuter, theta_l));
  end if:
  
  #----------------------------------------------
  # 3. Randomize.
  #----------------------------------------------

  if infolevel > 0 then
    PrintHeader("4. Randomizing the truncated system"):
  end if:

  # (a) ------------
  deg_variety := foldl(`*`, op( map(e -> degree(e), Et) )):
  D2 := floor( 6 * nops(theta_l) * deg_variety * (1 + 2 * d0 * max(op(beta))) / (1 - p_local) ):
  if infolevel > 1 then
    printf("%s %a\n", `Bound D_2 for assessing global identifiability: `, D2):
  end if:

  # (b, c) ---------
  sample := SamplePoint(D2, x_vars, y_vars, u_vars, mu, X_eq, Y_eq):
  while subs(sample[4], Q) = 0 do
    sample := SamplePoint(D2, x_vars, y_vars, u_vars, mu, X_eq, Y_eq):
  end do:    
  y_hat := sample[1]:
  u_hat := sample[2]:
  theta_hat := sample[3]:
  if infolevel > 1 then
    printf("%s %a\n", `Random sample for the outputs and inputs is generated from `, theta_hat):
  end if:

  # (d) ------------
  Et_hat := map(e -> subs([op(y_hat), op(u_hat)], e), Et):
  vars := { op(mu) }:
  for poly in Et_hat do
    vars := vars union { op(GetVars(poly, x_vars, s + 1)) }:
  end do:
  if infolevel > 1 then
    printf("%s %a %s %a %s\n", `The polynomial system \widehat{E^t} contains `, nops(Et_hat), `equations in `, nops(vars), ` variables`);
  end if:
  Q_hat := subs(u_hat, Q):

  #----------------------------------------------
  # 4. Determine.
  #----------------------------------------------

  if infolevel > 0 then
    PrintHeader("5. Assessing global identifiability"):
  end if:

  theta_g := []:
  if method = 1 then
    at_node := proc(var, args_node)
      local gb_loc, fname;
      if infolevel > 2 then
        fname := cat("SIAN_GB_computation_separated_for_", var):
        writedata(fname, 
          [
            "with(Groebner):",
            cat("polys := ", convert(args_node[1], string), ":"),
            cat("ordering := tdeg(op(", convert(convert(args_node[2], list), string), ")):"),
            "Basis(polys, ordering);"
          ], 
          string):
      end if:
      gb_loc := Groebner[Basis](op(args_node)):
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
    elif nops(theta_l) = 1 then
      # This is needed because of a bug in Grid[Seq]
      gb := [ at_node(theta_l[1], [
        [op(Et_hat), z * Q_hat - 1, (theta_l[1] - subs(theta_hat, theta_l[1])) * w - 1],
        tdeg(op(vars), z, w)
      ]) ]:
    end if:

    for i from 1 to nops(theta_l) do
      if gb[i] = [1] then
        theta_g := [op(theta_g), theta_l[i]]:
      else
        if infolevel > 1 then
          printf("%s %a %s %a\n", `Groebner basis corresponding to the parameter `, theta_l[i], ` is `, gb[i]):
        end if:
      end if:
    end do:    
  elif method = 2 then
    gb := Groebner[Basis]([op(Et_hat), z * Q_hat - 1], tdeg(op(vars), z));
    for i from 1 to nops(theta_l) do
      if Groebner[NormalForm](theta_l[i], gb, tdeg(op(vars), z)) = subs(theta_hat, theta_l[i]) then
        theta_g := [ op(theta_g), theta_l[i] ]:
      end if:
    end do:
  elif method = 3 then
    R := RegularChains[PolynomialRing](convert(vars, list)):
    for i from 1 to nops(theta_l) do
      tr := [RegularChains[Triangularize](Et_hat, [Q_hat, theta_l[i] - subs(theta_hat,theta_l[i])], R)]:
      for e in tr do
        print(RegularChains[Equations](e, R)):
      end do:
    end do:
  else
    print(`No such method`):
  end if:

  if infolevel > 0 then
    printf("\n=== Summary ===\n"):
    printf("%s %a\n", `Globally identifiable parameters:                 `, map(ParamToOuter, theta_g)):
    printf("%s %a\n", `Locally but not globally identifiable parameters: `, map(ParamToOuter, select(p -> not p in theta_g, theta_l))):
    printf("%s %a\n", `Not identifiable parameters:                      `, map(ParamToOuter, select(p -> not p in theta_l, theta))):
    printf("===============\n\n"):
  end if:

  table([
    globally = map(ParamToOuter, theta_g),
    locally_not_globally = map(ParamToOuter, select(p -> not p in theta_g, theta_l)),
    non_identifiable = map(ParamToOuter, select(p -> not p in theta_l, theta))
  ]):
end proc:

#===============================================================================
PrintHeader := proc(text):
#===============================================================================
  printf("\n=======================================================\n"):
  printf(text):
  printf("\n=======================================================\n"):
end proc:

#===============================================================================
GetParameters := proc(system_ODEs) local initial_values, all_symbols_rhs, mu:
#===============================================================================
  initial_values := map(f -> subs({t = 0}, int(f, t)), select( f -> type(int(f, t), function(name)), map(lhs, system_ODEs) )):
  all_symbols_rhs := foldl(`union`, op( map(e -> indets(rhs(e)), system_ODEs) )) minus {t}:
  mu := select(s -> not type(s, function), all_symbols_rhs):
  [op(mu), op(initial_values)]:
end proc:

#===============================================================================
FunctionToVariable := proc(f):
#===============================================================================
  convert(cat(convert(f, string)[1..-4], "_"), symbol):
end proc:

#===============================================================================
ParamToInner := proc(p) local s:
#===============================================================================
  s := convert(p, string):
  if length(s) > 3 and s[-3..-1] = "(0)" then
    MakeDerivative(FunctionToVariable(p), 0):
  else
    p:
  end if:
end proc:

#===============================================================================
ParamToOuter := proc(p) local s:
#===============================================================================
  s := convert(p, string):
  if length(s) > 2 and s[-2..-1] = "_0" then
    parse(cat(s[1..-3], "(0)")):
  else
    p:
  end if:
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
SamplePoint := proc(bound, x_vars, y_vars, u_vars, mu, X_eq, Y_eq)
#===============================================================================
  local n, m, s, all_params, all_vars, theta_hat, u_variables, 
        u_hat, x_hat, y_hat, eq, eq_num, eq_denom, 
        v, poly, i, j, all_subs, roll;
  n := nops(x_vars):
  m := nops(y_vars):
  s := nops(mu) + n:
  all_params := [op(mu), op(map(x -> MakeDerivative(x, 0), x_vars ))]:
  all_vars := [ op(x_vars), op(y_vars), op(u_vars) ]:

  roll := rand(0 .. bound):
  theta_hat := map(p -> p = roll(), all_params): 
  u_variables := [];
  for i from 1 to nops(u_vars) do
    u_variables := [ op(u_variables), seq(MakeDerivative(u_vars[i], j), j = 0..s) ]:
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
