#===============================================================================
IdentifiabilityODE := proc(system_ODEs, params_to_assess, {p := 0.99, count_solutions:=true, use_weights:=true, infolevel := 1, method := 2, num_nodes := 6, char:=0}) 
#===============================================================================
 local i, j, k, n, m, s, all_params, all_vars, eqs, Q, X, Y, poly, d0, D1, 
        sample, all_subs,alpha, beta, Et, x_theta_vars, prolongation_possible, 
        eqs_i, JacX, vars, vars_to_add, ord_var, var_index, deg_variety, D2, 
        y_hat, u_hat, theta_hat, Et_hat, Q_hat, theta_l, theta_g, gb, v, X_eq, Y_eq, 
        poly_d, separant, leader,vars_local, x_functions, y_functions, u_functions,
        all_symbols_rhs, mu, x_vars, y_vars, u_vars, theta, subst_first_order,
        subst_zero_order, x_eqs, y_eqs, param, other_params, to_add, at_node,
        prime, max_rank, R, tr, e, p_local, xy_ders, polys_to_process, new_to_process, solutions_table,
        Et_x_vars, var, G, P, output:

  #----------------------------------------------
  # 0. Extract inputs, outputs, states, and parameters from the system
  #----------------------------------------------

  if SearchText(".", convert(system_ODEs, string)) <> 0 then
    PrintHeader("WARNING: It looks like your system involves floating-point numbers. This may result into a non-meaninful result, please convert them to rationals (e.g., 0.2 -> 1/5)"):
  end if:

  if not verify(indets(system_ODEs, name), indets(system_ODEs), `subset`) then
    PrintHeader(cat("ERROR: you are using reserved maple symbols:", convert(indets(system_ODEs, name) minus indets(system_ODEs), string))):
    return;
  end if:

  randomize():

  if infolevel > 0 then
    PrintHeader("0. Extracting states, inputs, outputs, and parameters from the system"):
  end if:

  x_functions := map(f -> int(f, t), select( f -> type(int(f, t), function(name)), map(lhs, system_ODEs) )):
  y_functions := select( f -> not type(int(f, t), function(name)), map(lhs, system_ODEs) ):
  all_symbols_rhs := foldl(`union`, op( map(e -> indets(rhs(e)), system_ODEs) )) minus {t}:
  xy_ders := {op(x_functions), op(y_functions), op(select(f -> (f in all_symbols_rhs), map(lhs, system_ODEs)))}:
  u_functions := select( f -> type(f, function), convert(all_symbols_rhs minus xy_ders, list)):
  mu := convert(all_symbols_rhs minus {op(xy_ders), op(u_functions)}, list):

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

  if nops(y_functions) = 0 then
    PrintHeader("ERROR: no outputs in the model");
    return;
  end:

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
    poly_d := numer(lhs(x_eqs[i]) - rhs(x_eqs[i])):
    for j from 0 to s + 1 do
      leader := MakeDerivative(x_vars[i], j + 1):
      separant := diff(poly_d, leader):
      X[i] := [op(X[i]), poly_d]:
      X_eq := [op(X_eq), leader = -(poly_d - separant * leader) / separant]:
      poly_d := Differentiate(poly_d, all_vars):
    end do:
  end do:
  
  # (d,e) ---------------
  Y := []:
  Y_eq := []:
  for i from 1 to m do
    Y := [op(Y), []]:
    poly_d := numer(lhs(y_eqs[i]) - rhs(y_eqs[i])):
    for j from 0 to s + 1 do
      leader := MakeDerivative(y_vars[i], j):
      separant := diff(poly_d, leader):
      Y[i] := [op(Y[i]), poly_d]:
      Y_eq := [op(Y_eq), leader = -(poly_d - separant * leader) / separant]:
      poly_d := Differentiate(poly_d, all_vars):
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
  sample := SamplePoint(D1, x_vars, y_vars, u_vars, mu, X_eq, Y_eq, Q):
  all_subs := sample[4]:
  u_hat := sample[2]:
  y_hat := sample[1]:
 
  # (e) ------------------
  alpha := [seq(1, i = 1..n)]:
  beta := [seq(0, i = 1..m)]:
  Et := [];
  # TODO: improve for arbitrary derivatives
  x_theta_vars := all_params:
  prolongation_possible := [seq(1, i = 1..m)]:

  # (f) ------------------
  while foldl(`+`, op(prolongation_possible)) > 0 do
    for i from 1 to m do
      if prolongation_possible[i] = 1 then
        eqs_i := [op(Et), Y[i][beta[i] + 1]]:
        JacX := VectorCalculus[Jacobian](subs({op(u_hat), op(y_hat)}, eqs_i), x_theta_vars = subs(all_subs, x_theta_vars)):
        if LinearAlgebra[Rank](JacX) = nops(eqs_i) then
          Et := [op(Et), Y[i][beta[i] + 1]]:
          beta[i] := beta[i] + 1:
          # adding necessary X-equations
          polys_to_process := [op(Et), seq(Y[k][beta[k] + 1], k=1..m)]:
          while nops(polys_to_process) <> 0 do
            new_to_process := []:
            vars := {};
            for poly in polys_to_process do
              vars := vars union { op(GetVars(poly, x_vars)) }:
            end do:
            vars_to_add := { op(remove(v -> evalb(v in x_theta_vars), vars)) };
            for v in vars_to_add do
              x_theta_vars := [op(x_theta_vars), v];
              ord_var := GetOrderVar(v);
              var_index := ListTools[Search](ord_var[1], x_vars):
              poly := X[ var_index ][ ord_var[2] ]:
              Et := [op(Et), poly]:
              new_to_process := [op(new_to_process), poly]:
              alpha[ var_index ] := max(alpha[ var_index ], ord_var[2] + 1):
            end do:
            polys_to_process := new_to_process:
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
      for v in GetVars(Y[i][j], x_vars) do
        if not (v in x_theta_vars) then
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
 
  if infolevel > 0 then
    printf("%s %a\n", `Locally identifiable paramters: `, map(x -> ParamToOuter(x, all_vars), theta_l));
    printf("%s %a\n", `Nonidentifiable parameter: `, map(x -> ParamToOuter(x, all_vars), [op({op(theta)} minus {op(theta_l)})]));
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
  sample := SamplePoint(D2, x_vars, y_vars, u_vars, mu, X_eq, Y_eq, Q):
  y_hat := sample[1]:
  u_hat := sample[2]:
  # all_subs := sample[4]:
  theta_hat := sample[3]:
  if infolevel > 1 then
    printf("%s %a\n", `Random sample for the outputs and inputs is generated from `, theta_hat):
  end if:

  # (d) ------------
  Et_hat := map(e -> subs([op(y_hat), op(u_hat)], e), Et):
  Et_x_vars := {}:
  for poly in Et_hat do
    Et_x_vars := Et_x_vars union { op(GetVars(poly, x_vars)) }:
  end do:
  if infolevel > 1 then
    printf("%s %a %s %a %s\n", `The polynomial system \widehat{E^t} contains `, nops(Et_hat), `equations in `, nops(Et_x_vars) + nops(mu), ` variables`);
  end if:
  Q_hat := subs(u_hat, Q):

  vars := [
    op(sort([op(Et_x_vars)], (a, b) -> CompareDiffVar(a, b, x_vars))),
    z_aux, w_aux,
    op(sort(mu))
  ]:

  ###########
  non_id := map(x -> ParamToOuter(x, all_vars), [op({op(theta)} minus {op(theta_l)})]):
  if use_weights then
  	if infolevel > 0 then
    		PrintHeader("Applying Weighted Ordering", output_targets[log]):
    		LogText(sprintf("\t=> Applying Weighted Ordering"), ProgressBar):
  	end if:
    input_table := table(
      [
        "sigma"=system_ODEs,
        "poly_system" = [op(Et_hat), z_aux * Q_hat - 1], "poly_vars" = vars, "non_id" = non_id, 
        "s"=s, "m"=m, "x_vars"=x_vars, "y_vars"=y_vars,
        "mu"=mu, "x_eqs"=x_eqs, "y_eqs"=y_eqs, "all_vars"=all_vars
      ]
    );
    weight_subs, poly_system := SubsByDepth(input_table):
    Et_hat := poly_system;
    weights_table := table(weight_subs);
  else
    weights_table := table([seq(_var_=_var_, _var_ in vars)]);
  	Et_hat := [op(Et_hat), z_aux * Q_hat - 1]:
  end if;
  ###########
  
  if infolevel > 1 then
    printf("Variable ordering to be used for Groebner basis computation %a\n", vars);
    printf("%s %a\n", `Weight assignment:`, [entries(weights_table, `pairs`)]);
  end if:

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
            cat("ordering := tdeg(op(", convert(args_node[2], string), ")):"),
            "Basis(polys, ordering);"
          ], 
          string):
      end if:
      gb_loc := Groebner[Basis](op(args_node)):
      gb_loc;
    end proc:

    if nops(theta_l) > 1 and num_nodes > 1 then
      Grid[Setup]("local", numnodes = num_nodes):
      Grid[Set](at_node):
      gb := Grid[Seq](
        at_node(theta_l[i], [
          [op(Et_hat), z_aux * Q_hat - 1, (theta_l[i] - subs(theta_hat, theta_l[i])) * w_aux - 1],
          tdeg(op(vars))
        ]),
        i = 1..nops(theta_l)
      ):
    else
      gb := []:
      for i from 1 to nops(theta_l) do
        gb := [
          op(gb), 
          at_node(
            theta_l[i], 
            [[op(Et_hat), z_aux * Q_hat - 1, (theta_l[i] - subs(theta_hat, theta_l[i])) * w_aux - 1], tdeg(op(vars))]
           ) 
        ]:
      end do:
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
    gb := Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=char);
    for i from 1 to nops(theta_l) do
      if char>0 then
        check := subs(theta_hat, theta_l[i]) mod char:
      else
        check := subs(theta_hat, theta_l[i]):
      end if:
      if Groebner[NormalForm](theta_l[i]^degree(weights_table[theta_l[i]]), gb, tdeg(op(vars)), characteristic=char) = check then
        theta_g := [op(theta_g), theta_l[i]]:
      end if:
    end do:

    if count_solutions then 
      solutions_table := table([]):
      for var in theta_g do
        if infolevel > 0 then
          printf("%s %a %s %a\n",`The number of solutions for`, var, `is`, 1):
        end if:
        solutions_table[var] := 1:
      end do:
	
      for var in select(p -> not p in theta_g, theta_l) do
        G := Groebner[Walk](gb, tdeg(op(vars)), lexdeg([op({op(vars)} minus {var})], [var])):
	      P := select(x->evalb(indets(x)={var}), G):
	      solutions_table[var]:=degree(P[1], [op(indets(P))])/weights_table[var]: 
        if infolevel > 1 then
          printf("%s %a %s %a\n",`The number of solutions for`, var, `is`, degree(P[1], [op(indets(P))])):
        end if:
      end do:
    end if:  
  else
    print(`No such method`):
  end if:

  if infolevel > 0 then
    printf("\n=== Summary ===\n"):
    printf("%s %a\n", `Globally identifiable parameters:                 `, map(x -> ParamToOuter(x, all_vars), theta_g)):
    printf("%s %a\n", `Locally but not globally identifiable parameters: `, map(x -> ParamToOuter(x, all_vars), select(p -> not p in theta_g, theta_l))):
    printf("%s %a\n", `Not identifiable parameters:                      `, map(x -> ParamToOuter(x, all_vars), select(p -> not p in theta_l, theta))):
    printf("===============\n\n"):
  end if:

  output := table([
    globally = {op(map(x -> ParamToOuter(x, all_vars), theta_g))},
    locally_not_globally = {op(map(x -> ParamToOuter(x, all_vars), select(p -> not p in theta_g, theta_l)))},
    non_identifiable = {op(map(x -> ParamToOuter(x, all_vars), select(p -> not p in theta_l, theta)))}
  ]):

  if count_solutions then 
     PrintHeader("WARNING: The result of solution counting is guaranteed with high probability, however it NOT the same probability 'p' as provided in the input."):
     output[num_solutions] := eval(solutions_table):
  end if:

  return output;
end proc:

#===============================================================================
PrintHeader := proc(text):
#===============================================================================
  printf("\n=======================================================\n"):
  printf(text):
  printf("\n=======================================================\n"):
end proc:

#===============================================================================
GetParameters := proc(system_ODEs, {initial_conditions := true}) local initial_values, all_symbols_rhs, mu:
#===============================================================================
  initial_values := map(f -> subs({t = 0}, int(f, t)), select( f -> type(int(f, t), function(name)), map(lhs, system_ODEs) )):
  all_symbols_rhs := foldl(`union`, op( map(e -> indets(rhs(e)), system_ODEs) )) minus {t}:
  mu := select(s -> not type(s, function), all_symbols_rhs):
  if initial_conditions then
    return [op(mu), op(initial_values)]:
  else
    return [op(mu)]:
  end if:
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
ParamToOuter := proc(p, varnames) local s:
#===============================================================================
  s := convert(p, string):
  if length(s) > 2 and s[-2..-1] = "_0" and parse(s[1..-2] )in varnames then
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
Differentiate := proc(diff_poly, var_list) 
#===============================================================================
  local result, aux, v, h, diff_v:
  result := 0:
  for diff_v in indets(diff_poly) do
    aux := GetOrderVar(diff_v):
    # seems that Maple does not have unpacking
    v := aux[1]:
    h := aux[2]:
    if v in var_list then
      result := result + diff(diff_poly, MakeDerivative(v, h)) * MakeDerivative(v, h + 1):
    end if:
  end do:
  simplify(result):
end proc:

#===============================================================================
GetVars := proc(diff_poly, var_list)
#===============================================================================
  local result;
  result := select(v -> evalb(GetOrderVar(v)[1] in var_list), indets(diff_poly)):
  return result:
end proc:


#===============================================================================
GetOrderVar := proc(diff_var)
#===============================================================================
  local s, v, h;
  if not StringTools[RegMatch]("^(.*_)([0-9]+)$", diff_var, s, v, h) then
    return ["", ""]:
  end if:
  return [parse(v), parse(h)]:
end proc:


#===============================================================================
SamplePoint := proc(bound, x_vars, y_vars, u_vars, mu, X_eq, Y_eq, Q)
#===============================================================================
  local n, m, s, all_params, all_vars, theta_hat, u_variables, 
        u_hat, x_hat, y_hat, eq, eq_num, eq_denom, 
        v, poly, i, j, all_subs, roll, to_compute;
  n := nops(x_vars):
  m := nops(y_vars):
  s := nops(mu) + n:
  all_params := [op(mu), op(map(x -> MakeDerivative(x, 0), x_vars ))]:
  all_vars := [ op(x_vars), op(y_vars), op(u_vars) ]:

  roll := rand(0 .. bound):
  while true do
    theta_hat := map(p -> p = roll(), all_params): 
    u_variables := [];
    for i from 1 to nops(u_vars) do
      u_variables := [ op(u_variables), seq(MakeDerivative(u_vars[i], j), j = 0..s + 1) ]:
    end do:
    u_hat := map(p -> p = roll(), u_variables) :   
  
    all_subs := [op(theta_hat), op(u_hat)]:
    if subs(all_subs, Q) = 0 then
      next
    end if:
    to_compute := [op(X_eq), op(Y_eq)]:
    while nops(to_compute) <> 0 do
      to_compute := map(e -> lhs(e) = subs(all_subs, rhs(e)), to_compute);
      all_subs := [ op(all_subs), op(select(e -> type(rhs(e), numeric), to_compute)) ]:
      to_compute := remove(e -> type(rhs(e), numeric), to_compute):
    end do:
    y_hat := map(e -> lhs(e) = subs(all_subs, rhs(e)), Y_eq):
    x_hat := map(e -> lhs(e) = subs(all_subs, rhs(e)), X_eq):
    break:
  end do:

  return [y_hat, u_hat, theta_hat, all_subs];
end proc:

#===============================================================================
GenerateReplica := proc(equations, r)
#===============================================================================
  # generates a system of equations corresponding to r independent trajectories of
  # the original system. Time-dependent variabes are replicated, parameters are not
  local all_functions, zero_order, first_order, funcs, without_t, result, i, subst:
  all_functions := select(f -> type(f, function), foldl(`union`, op( map(indets, equations) ))):
  zero_order := select(f -> not type(int(f, t), function(name)), all_functions):
  first_order := map(f -> int(f, t), select(f -> type(int(f, t), function(name)), all_functions)):
  funcs := {op(zero_order), op(first_order)}:
  without_t := map(f -> convert(convert(f, string)[1..-4], symbol), funcs):
  result := []:
  for i from 1 to r do
    subst := map(f -> f = convert(cat(convert(f, string), "_r", i), symbol), without_t):
    result := [op(result), op(map(e -> subs(subst, e), equations))]:
  end do: 
  return result:
end proc:

#===============================================================================
CompareDiffVar := proc(dvl, dvr, var_list)
#===============================================================================
  local vl, vr, hl, hr;
  vl, hl := op(GetOrderVar(dvl, var_list)):
  vr, hr := op(GetOrderVar(dvr, var_list)):
  if evalb(hl <> hr) then
    return evalb(hl > hr):
  end if:
  if evalb(length(vl) <> length(vr)) then
    return evalb(length(vl) > length(vr)):
  end if:
  return StringTools[Compare](vr, vl):
end proc:


#===============================================================================
# Weights
#===============================================================================

# get rid of unerscore, e.g. x1(t) -> x1_ -> x1
get_function_name := f -> parse(convert(FunctionToVariable(f), string)[..-2]):

# check if depends on (t)
is_function:= f->StringTools[Has](convert(f, string), "(t)"):
idtfm := x->x:

# check if is derivative
is_diff := f->type(int(f, t), function(name)):

lhs_name := ff -> if convert(ff, string)[-1] = "_" then parse(convert(ff, string)[..-2]) else ff; end if:

#===============================================================================
GetStateName := proc(state, x_vars, mu)
#===============================================================================
  local state_;
  if state in mu then
    return state;
  end if;
  state_ := parse(cat(StringTools[Join](StringTools[Split](convert(state, string), "_")[..-2], "_"), "_")):
  if state_ in x_vars then
    return state_;
  end if:
end proc:

#===============================================================================
GetMinLevelBFS := proc(s, m, x_vars, y_vars, mu, x_eqs, y_eqs_in, all_vars)
# s: number of states and parameters
# m: number of outputs
#===============================================================================
  # this part is copied from original SIAN code
  local current_level, visible_states, visibility_table, i, j, continue, poly_d,
  leader, separant, candidates, each, differentiate_, k, y_eqs:
  y_eqs := y_eqs_in:
  current_level := 0:
  # get functions on level 0, we consider parameters and states indistinguishable
  # i.e. parameters are states with d/dt = 0
  visible_states :=  foldl(`union`, op(map(x->indets(rhs(x)) minus {t}, y_eqs))); #select(f -> f in x_zero_vars, ); # map(x->parse(convert(x, string)[..-2]), select(f -> f in x_zero_vars, foldl(`union`, op(map(x->indets(rhs(x)), y_eqs))))):# cat(StringTools[Split](convert(x, string), "_")[1], "_")

  # construct a hash table of "visibility"
  visibility_table := table([seq(GetStateName(each, x_vars, mu)=current_level, each in visible_states)]):
  # this is a flag array: if i-th position == 1 then we must differentiat i-th y(t) function 
  differentiate_ := [seq(1, i=1..nops(y_eqs))]: 

  for j from 1 to s + 1 do
    # begin differentiation
    current_level := current_level + 1:
    continue:=true:

    for i from 1 to m do
      if differentiate_[i]=1 then 
        poly_d := numer(lhs(y_eqs[i]) - rhs(y_eqs[i])):
        poly_d := Differentiate(poly_d, all_vars):
        leader := MakeDerivative(y_vars[i], j):
        separant := diff(poly_d, leader):
        poly_d := simplify(leader - subs(x_eqs, -(poly_d - separant * leader) / separant)):
        y_eqs[i]:= leader = simplify(poly_d - leader);
        candidates := select(x-> not (GetOrderVar(x)[1]  in y_vars), indets(y_eqs[i])):
        if op(map(x->not assigned(visibility_table[GetStateName(x, x_vars, mu)]), candidates)) <> NULL then
          continue := foldl(`or`, op(map(x->not assigned(visibility_table[GetStateName(x, x_vars, mu)]), candidates))):
        else
          continue := false;
        fi:
        if continue then
          differentiate_[i]:=1:
        else
          differentiate_[i]:=0:
        fi:
        for each in candidates do
          if not assigned(visibility_table[GetStateName(each, x_vars, mu)]) then 
              visibility_table[GetStateName(each, x_vars, mu)] := current_level:
          fi:
        od;
      fi:
    od:
    if add(k, k in differentiate_)=0 then
        break:
    fi:
  od:
  return visibility_table;
end proc:

#===============================================================================
SubsByDepth := proc(input_table, {trdegsub:=false})
# input_table is a Maple table with key value pairs:
# sigma: input ODE system
# poly_system: polynomial system Et_hat
# poly_vars: variables of Et_hat in the proper order
# non_id: non-identifiable parameters
# s: number of parameters + number of states in sigma
# m: number of outputs
# x_vars: state variables
# y_vars: output variables (y-functions)
# mu: list of parameters from sigma
# x_eqs: ODEs from sigma in polynomial form
# y_eqs: output equations from sigma in polynomial form
# all_vars: all variables in sigma
#===============================================================================
  local counting_table_states, min_count, vts, rhs_terms, max_possible,
        rhs_term, indets_, term, substitutions, each, alg_indep,
        all_subs, names, selection, other, all_odes, each_ode;

  vts := GetMinLevelBFS(
    input_table["s"],
    input_table["m"],
    input_table["x_vars"],
    input_table["y_vars"],
    input_table["mu"],
    input_table["x_eqs"],
    input_table["y_eqs"],
    input_table["all_vars"]
  ):
  substitutions := table([]);
  all_odes := map(x->expand(rhs(x)), select(f->is_diff(lhs(f)), input_table["sigma"]));
  rhs_monoms := []:
  for each_ode in all_odes do
    if whattype(each_ode) in [`+`,`*`,`^`] then
      rhs_monoms := [op(rhs_monoms), op(each_ode)];
    end if;
    if whattype(each_ode) in [function] then
     rhs_monoms := [op(rhs_monoms), (each_ode)]; 
    end if;
  end do;
  
  max_possible := max(map(rhs, [entries(vts, `pairs`)]));
  for rhs_term in rhs_monoms do
    indets_ := convert(indets(rhs_term) minus {t}, list):
    for term in indets_ do
      if is_function(term) then
        if assigned(vts[FunctionToVariable(term)]) then
          substitutions[FunctionToVariable(term)] := vts[FunctionToVariable(term)]+1:
        end if;
      else
        if not term in input_table["non_id"] and vts[term]=max_possible and assigned(vts[term]) then # 
          substitutions[term] := vts[term]+1:
        end if;
      end if:
    end do;
  end do:
  substitutions[z_aux]:=min(3, max_possible):

  new_et_hat := input_table["poly_system"]:
  all_subs := {}:
  names := [indices(substitutions, `nolist`)];
  for each in names do #system_vars[2] do
    selection := select(sys_var->StringTools[IsPrefix](convert(each, string), sys_var), input_table["poly_vars"]);
    for other in selection do
        new_et_hat := subs({other = other^substitutions[each]}, new_et_hat):
        all_subs := all_subs union {other = other^substitutions[each]}:
    end do;
  od:
  return all_subs, new_et_hat;
end proc:

