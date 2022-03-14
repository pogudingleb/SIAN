#==============================================================================
IdentifiabilityODE := proc(system_ODEs, params_to_assess, system_name, {sub_transc:=true, p := 0.99, max_comb:=1000, count_solutions:=true, infolevel := 2, method := 2, num_nodes := 6}) 
#===============================================================================
 local i, j, k, n, m, s, all_params, all_vars, eqs, Q, X, Y, poly, d0, D1, 
        sample, all_subs,alpha, beta, Et, x_theta_vars, prolongation_possible, 
        eqs_i, JacX, vars, vars_to_add, ord_var, var_index, deg_variety, D2, 
        y_hat, u_hat, theta_hat, Et_hat, Q_hat, theta_l, theta_g, gb, v, X_eq, Y_eq, 
        poly_d, separant, leader,vars_local, x_functions, y_functions, u_functions,
        all_symbols_rhs, mu, x_vars, y_vars, u_vars, theta, subst_first_order,
        subst_zero_order, x_eqs, y_eqs, param, other_params, to_add, at_node,
        prime, max_rank, R, tr, e, p_local, xy_ders, polys_to_process, new_to_process, solutions_table,
        Et_x_vars, var, G, P, output, alg_indep, rrefJacX, pivots, row_idx, row, pivot_idx, non_id, faux_equations,
        y_faux, alg_indep_derivs, alg_indep_params, faux_odes, faux_outputs, derivs, x_theta_vars_, sigma_new, identifiable_states,
        x_theta_vars_filtered, x_theta_vars_to_be_removed, number_of_choices, choices, current_choice, choice_idx,
        vars_old, Et_hat_old, et_hat_monomials, sum_degrees, each, entropies, entr, rhs_cols, idxs, parameter, idx, A, in_span_count,
        solution, sum_degrees_new, set_of_mnm, counter, count_greater_than_ones, entr_new, Et_hat_subs:
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
  alg_indep := {}:
  # x_theta_vars := ListTools[Reverse](x_theta_vars):
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
              x_theta_vars := [v, op(x_theta_vars)];
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
  pivots := {}:
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
  derivs := {op(x_theta_vars)} minus {op(mu)};
  non_id := [op({op(theta)} minus {op(theta_l)})]:

  x_theta_vars_ := ListTools[Reverse]([op({op(x_theta_vars)} minus {op(theta_l)})]);
  x_theta_vars := [op(theta_l), op(x_theta_vars_)];
  if sub_transc then
    JacX := VectorCalculus[Jacobian](subs({op(u_hat), op(y_hat)}, Et), x_theta_vars = subs(all_subs, x_theta_vars)); 
    rrefJacX := LinearAlgebra[ReducedRowEchelonForm](JacX):
    pivots := {}:
    for row_idx from 1 to nops(Et) do #nops(theta) do
        row := rrefJacX[row_idx]:
        pivot_idx := 1:
        while row[pivot_idx]=0 and add(row)<>0 do
          pivot_idx := pivot_idx + 1:
        end do:
        if pivot_idx <= numelems(row) then
          pivots := {op(pivots), x_theta_vars[pivot_idx]}:
        end if:
    end do:
    alg_indep := {op(x_theta_vars)} minus pivots:
  end if:
  sigma_new := system_ODEs:
  if infolevel > 0 then
    printf("%s %a\n", `Locally identifiable paramters: `, map(x -> ParamToOuter(x, all_vars), theta_l));
    printf("%s %a\n", `Nonidentifiable parameter: `, map(x -> ParamToOuter(x, all_vars), [op({op(theta)} minus {op(theta_l)})]));
  end if:

  identifiable_states := map(each->GetOrderVar(each)[1], {op(select(each->GetOrderVar(each)[1] <>"", theta_l))}); # pick names of states whose IC is identifiable
  x_theta_vars_to_be_removed := select(each-> GetOrderVar(each)[1] in identifiable_states, x_theta_vars_);
  x_theta_vars_filtered := [op({op(x_theta_vars_)} minus {op(x_theta_vars_to_be_removed)})];
  number_of_choices:= binomial(numelems(x_theta_vars_filtered), numelems(alg_indep));
  if infolevel > 0 then
    printf("%s %a\n", `Number of transcendence basis choices before \"filtering\":`, binomial(numelems(x_theta_vars_), numelems(alg_indep)));
    printf("%s %a\n", `Number of transcendence basis choices after \"filtering\":`, number_of_choices);
  end if:
  if number_of_choices < max_comb then
    choices := combinat[choose](x_theta_vars_filtered, numelems(alg_indep));
  else
    current_choice := [op(alg_indep)];
    choices := [current_choice];
    while numelems(choices) < max_comb do
      choices := {op(choices), combinat[randperm](x_theta_vars_filtered)[..numelems(alg_indep)]};
    end do;
  end if;
  choice_idx := 1;
  if infolevel > 0 then
    PrintHeader("Creating Truncated system (old)"):
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

  vars_old := [
    op(sort([op(Et_x_vars)], (a, b) -> CompareDiffVar(a, b, x_vars))),
    z_aux, w_aux,
    op(sort(mu))
  ]:

  Et_hat_old := [op(Et_hat), z_aux*Q_hat - 1]:
  et_hat_monomials := map(each->op(expand(each)), Et_hat_old);
  sum_degrees := 0;
  degrees := [];
  seen := []:
  for each in et_hat_monomials do
    for param in alg_indep do # {x, y}
        if param in {op(each)} then # {y}, seen = [x*y^2]
            sum_degrees := sum_degrees+degree(each);
            degrees := [op(degrees), degree(each)];
            # if each in seen then
            #   # emphasize the double-count
            # end if:
        end if:
    end do:
  end do:
  entropies := [];
  entr := 0:
  writeto(cat(system_name, "/bad_pairs.out"));
  printf("Bad pairs for %s system\n", system_name);
  writeto(terminal);
  for alg_indep in choices do  
    if sub_transc and numelems(alg_indep)<>0 then
      printf("%s %a\n", `Algebraically independent parameters (candidates)`, map(x-> ParamToOuter(x, all_vars), alg_indep)):
      if sub_transc then 
        PrintHeader("Substituting transcendence basis."):
      end if:
      rhs_cols := []:
      idxs := []:
      for parameter in alg_indep do
        idx := ListTools[Search](parameter, x_theta_vars):
        idxs := [op(idxs), idx]:
        rhs_cols := [op(rhs_cols), JacX[..,idx]]:
      end do:
      rhs_cols := Matrix(rhs_cols):
      A := LinearAlgebra[DeleteColumn](JacX, idxs):

      # test column space membership:
      in_span_count := 0:
      #for col_idx from 1 to numelems(alg_indep) do
      try
        solution := LinearAlgebra[LinearSolve](A, rhs_cols):  # Ax = b, b=rhs_cols
      catch :
        appendto(cat(system_name, "/bad_pairs.out")):
        printf("%s %a %s\n", "Pair", alg_indep, "is not transcendence basis" ):
        writeto(terminal):
        choice_idx:=choice_idx+1:
        next
      end try;

      alg_indep_derivs := {op(alg_indep)} intersect derivs:
      alg_indep_params := ({op(alg_indep)} intersect {op(non_id)}) minus {op(alg_indep_derivs)}:
      
      faux_outputs := []: # seq(parse(cat("y_faux", idx, "(t)"))=alg_indep_params[idx](t), idx in 1..numelems(alg_indep_params))
      faux_odes := []:
      idx := 1:
      sigma_new := system_ODEs:
      for each in alg_indep_params do
        if not (each in x_vars) then
          sigma_new := subs({each=each(t)}, sigma_new):
          faux_outputs := [op(faux_outputs), parse(cat("y_faux", idx, "(t)"))=each(t)]:
          faux_odes := [op(faux_odes), diff(each(t), t)=0]:
        else
          faux_outputs := [op(faux_outputs), parse(cat("y_faux", idx, "(t)"))=parse(convert(each, string)[..-2])(t)]:
        end if:
        idx := idx+1:
      end do:

      sigma_new := [op(faux_odes), op(sigma_new), op(faux_outputs)]:
      
      if infolevel>0 then
        printf("%s %a\n", `Algebraically independent parameters among nonidentifiable:`, map(x-> ParamToOuter(x, all_vars), alg_indep_params)):
        printf("%s %a\n", `Algebraically independent parameters among derivatives:`, map(x-> ParamToOuter(x, all_vars), alg_indep_derivs)):
      end if:
      
      if infolevel>1 then
        printf("\t%s %a\n", `Adding ODEs:`, faux_odes):
        printf("\t%s %a\n", `Adding output functions:`, faux_outputs):
        printf("\t%s %a\n", `New system:`, sigma_new):
      end if:
      X_eq, Y_eq, Et, theta_l, x_vars, y_vars, mu, beta, Q, d0 := PreprocessODE(sigma_new, GetParameters(sigma_new)):
      if numelems(alg_indep_derivs)>0 then
        faux_equations := [seq(parse(cat("y_faux", idx+numelems(alg_indep_params), "_0"))=alg_indep_derivs[idx], idx in 1..numelems(alg_indep_derivs))]:
        y_faux := [seq(parse(cat("y_faux", idx+numelems(alg_indep_params), "_")), idx=1..numelems(alg_indep_derivs))]:
        Et := [op(Et), op(map(x->lhs(x)-rhs(x), faux_equations))]:
        Y_eq := [op(Y_eq), op(faux_equations)]:
        if infolevel>1 then
          printf("\t%s %a\n", `Adding new y-equations:`, faux_equations):
          printf("\t%s %a\n", `New system:`, Et):
          printf("\t%s %a\n", `New system:`, Y_eq):
        end if:
      end if; 
    else
      printf("%s\n", `No algebraically independent parameters found.`);
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

    sum_degrees_new := 0;
    degree_table := table([seq(param = [], param in alg_indep)]);
    occurrence_table := table([seq(param = 0, param in alg_indep)]):
    for each in et_hat_monomials do
      for param in alg_indep do # {beta_, q}
          if param in {op(each)} then
              sum_degrees_new := sum_degrees_new+degree(each);
              occurrence_table[param] += 1;
              degree_table[param] := [op(degree_table[param]), degree(each)];
          end if:
      end do:
    end do:

    Et_hat_subs := subs({op(map(each -> each=1, alg_indep))}, map(expand, Et_hat_old));
    set_of_mnm:=[op(map(x->x/max(1, Groebner[LeadingCoefficient](x, tdeg(op(vars_old)))), map(x->op(x), Et_hat_subs)))];
    counter := table([]):
    for m in set_of_mnm do
        if assigned(counter[m]) then
            counter[m] := counter[m]+1:
        else
            counter[m] := 1:
        end if:
    end do:
    counter := sort([entries(counter, `pairs`)], (a,b)->rhs(a)>rhs(b)):
    count_greater_than_ones := numelems(select(x->rhs(x)>1, counter));
    entr_new:=convert(add(map(x->rhs(x)/add(map(rhs, counter))*log((rhs(x)/add(map(rhs, counter)))), counter)), float);
    entropies := [op(entropies), [alg_indep, numelems(counter), entr_new]];
    if sub_transc then
      # if entr_new<entr then
        entr := entr_new;
        system(cat("mkdir -p ", system_name));
        print(cat(system_name, "/", "transcendence_", convert(choice_idx, string), "_entropy_heuristic.mpl"));
        if numelems(alg_indep) = 0 then
          writeto(cat(system_name, "/", "transcendence_", convert(choice_idx, string), ".mpl")):
          printf("%s\n", `No Algebraically independent parameters`):
          writeto(cat(system_name, "/", "transcendence.mpl")):
          printf("%s\n", `quit;`):
        else
          writeto(cat(system_name, "/", "transcendence_dc_", convert(choice_idx, string), "_entropy_heuristic.mpl")):
          printf("infolevel[Groebner]:=10;\nEt_hat := %a;\nvars:=%a;\n", [op(Et_hat), z_aux * Q_hat - 1], vars);
          printf("gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');\n");
          printf("# %a %a\n# %a, %a\n# %a, %a\n# %a\n", alg_indep, sum_degrees_new, numelems(counter), counter, add(map(rhs, counter)), entr, count_greater_than_ones);
          printf("# %a\n# %a", [entries(occurrence_table, `pairs`)], [entries(degree_table, `pairs`)]);
        end if:
    else
      writeto(cat(system_name, "/", "no_transcendence.mpl")):
      printf("infolevel[Groebner]:=10;\nEt_hat := %a;\nvars:=%a;\n", [op(Et_hat), z_aux * Q_hat - 1], vars);
      printf("gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');\n");
      # writeto(cat(system_name, "/", "no_transcendence.m")):
      # printf("SetNthreads(64);\nQ:= GF(11863279);\nSetVerbose(\"Faugere\", 2);\nP<%s>:= PolynomialRing(Q, %d, \"grevlex\");\nG := [%s];\n", convert(vars, string)[2..-2], nops(vars), convert([op(Et_hat), z_aux * Q_hat - 1], string)[2..-2]);
      # printf("time GroebnerBasis(G);\nquit;");
    end if:
    writeto(terminal);
    choice_idx := choice_idx+1;
  end do;
  # Sorting :=proc(a, b)
  #   if a[2] < b[2] then
  #     return true;
  #   elif a[2] > b[2] then
  #     return false;
  #   elif a[3] < b[3] then
  #     return true;
  #   end if;
  #   return false;
  # end proc;
  # entropies := sort(entropies, (a,b)-> Sorting(a,b));
  appendto(cat(system_name, "/", "entropies.mpl"));
  for each in entropies do
    printf("%a\n",each);
  end do;
  writeto(terminal);
  return;

  if infolevel > 1 then
    printf("Variable ordering to be used for Groebner basis computation %a\n", vars);
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
    gb := Groebner[Basis]([op(Et_hat), z_aux * Q_hat - 1], tdeg(op(vars)));
    for i from 1 to nops(theta_l) do
      if Groebner[NormalForm](theta_l[i], gb, tdeg(op(vars))) = subs(theta_hat, theta_l[i]) then
        theta_g := [ op(theta_g), theta_l[i] ]:
      end if:
    end do:
    theta_g := select(x-> not (x in alg_indep), theta_g):
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
	solutions_table[var]:=degree(P[1], [op(indets(P))]): 
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

get_function_name := f -> parse(convert(FunctionToVariable(f), string)[..-2]): 
GetStateName := proc(state, x_vars, mu) # i.e. Xij_q => Xij
  local state_;
  if state in mu then
    return state;
  end if;
  state_ := parse(cat(StringTools[Join](StringTools[Split](convert(state, string), "_")[..-2], "_"), "_")):
  if state_ in x_vars then
    return state_; #parse(convert(state_, string)[..-2]);
  end if:
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

#==============================================================================
PreprocessODE := proc(system_ODEs, params_to_assess, {p := 0.99, infolevel := 1, method := 2}) 
#===============================================================================
 local i, j, k, n, m, s, all_params, all_vars, eqs, Q, X, Y, poly, d0, D1, 
        sample, all_subs,alpha, beta, Et, x_theta_vars, prolongation_possible, 
        eqs_i, JacX, vars, vars_to_add, ord_var, var_index, deg_variety, D2, 
        y_hat, u_hat, theta_hat, Et_hat, Q_hat, theta_l, theta_g, gb, v, X_eq, Y_eq, 
        poly_d, separant, leader,vars_local, x_functions, y_functions, u_functions,
        all_symbols_rhs, mu, x_vars, y_vars, u_vars, theta, subst_first_order,
        subst_zero_order, x_eqs, y_eqs, param, other_params, to_add, at_node,
        prime, max_rank, R, tr, e, p_local, xy_ders, polys_to_process, new_to_process, solutions_table,
        Et_x_vars, var, G, P, output, alg_indep, rrefJacX, pivots, row_idx, row, pivot_idx, non_id, faux_equations,
        y_faux:

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

  # if infolevel > 0 then
  #   PrintHeader("0. Extracting states, inputs, outputs, and parameters from the system"):
  # end if:

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

  # if infolevel > 0 then
  #   PrintHeader("1. Constructing the maximal polynomial system"):
  # end if:

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

  # if infolevel > 0 then
  #   PrintHeader("2. Truncating the polynomial system based on the Jacobian condition"):
  # end if:

  # (a) ---------------
  d0 := max(op( map(f -> degree( simplify(Q * rhs(f)) ), eqs) ), degree(Q)):

  # (b) ---------------
  # extra factor nops(theta) + 1 compared to the formula in the paper is to
  # provide probability gaurantee to the local identifiability test
  D1 := floor( (nops(theta) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p_local) ):
  # prime := nextprime(D1):
  # if infolevel > 1 then
  #   printf("%s %a\n", `Bound D_1 for testing the rank of the Jacobian probabilistically: `, D1);
  # end if:

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
  alg_indep := {}:
  x_theta_vars := ListTools[Reverse](x_theta_vars):
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
              x_theta_vars := [v, op(x_theta_vars)];
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
 
  # if infolevel > 1 then
  #   printf("%s %a\n", `Orders of prolongations of the outputs (beta) = `, beta):
  #   printf("%s %a\n", `Orders of prolongations of the state variables (alpha) = `, alpha):
  # end if:
 
  ##############################

  # if infolevel > 0 then
  #   PrintHeader("3. Assessing local identifiability"):
  # end if:
 
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
 return X_eq, Y_eq, Et, theta_l, x_vars, y_vars, mu, beta, Q, d0:
 end proc: