%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTSOLVERINTERFACE
% DATE: 20/05/2023
% AUTHOR: BRUNO BARBOSA DO AMARAL
%
% DESCRIPTION: BUILDS AN OPTIMIZATION PROBLEM IN MEMORY AND PROVIDES
% ITS REPRESENTATION TO SPECIFIC OPTIMIZATION SOFTWARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Macros
VAR_TYPE_CONT  = 0;
CTR_TYPE_EQUAL = 0;

function [upt_struct, var_pointer] = add_var(lower_bound, upper_bound, obj_val, var_type, struct)

  % Copy struct content
  upt_struct = struct;

  % Increment structure size
  var_pointer = upt_struct.size;
  upt_struct.size = upt_struct.size + 1;

  upt_struct.lb(upt_struct.size)   = lower_bound;
  upt_struct.ub(upt_struct.size)   = upper_bound;
  upt_struct.obj(upt_struct.size)  = obj_val;
  upt_struct.type(upt_struct.size) = var_type;

end

function [upt_struct, ctr_pointer] = add_ctr(lower_bound, upper_bound, ctr_type, struct)

  % Copy struct content
  upt_struct = struct;

  % Increment structure size
  ctr_pointer = upt_struct.size;
  upt_struct.size = upt_struct.size + 1;

  upt_struct.lb(upt_struct.size)   = lower_bound;
  upt_struct.ub(upt_struct.size)   = upper_bound;
  upt_struct.type(upt_struct.size) = ctr_type;

end

function upt_struct = add_var_to_ctr(coef, var_pointer, ctr_pointer, struct)

  % Copy struct content
  upt_struct = struct;

  upt_struct.var(ctr_pointer,var_pointer) = coef;
end

function [opt_sol, opt_obj, error, info] = solve_glpk_prob(sense, optvar, optctr)

  % Temp
  VAR_TYPE_CONT  = 0;
  CTR_TYPE_EQUAL = 0;

  % 1 for min/-1 for max
  prob_type = sense;

  % Matrices
  A = optctr.var;
  b = [];

  ctr_type = "";
  for i = 1:optctr.size
    % Configuring constraint type (= -> S,>= -> L,<= -> U, <= >= -> D, Free -> F)
    if (optctr.type == CTR_TYPE_EQUAL)
      ctr_type = strcat(ctr_type,"S");
      b(i) = optctr.lb(i);
    end
  end

  var_type = "";
  for i = 1:optvar.size
    % Configuring variable type
    if (optvar.type == VAR_TYPE_CONT)
      var_type = strcat(var_type,"C");
    end
  end
  var_lower_bounds = optvar.lb;
  var_upper_bounds = optvar.ub;

  % Other parameters
  param.msglev=1;
  param.itlim=100;

  % Solving optimization problem
  [opt_sol, opt_obj, error, info] = glpk(optvar.obj, A, b,...
   var_lower_bounds, var_upper_bounds, ctr_type, var_type, prob_type, param);

end

% Add variables
optvar.size = 0;
optvar.lb   = [];
optvar.lb   = [];
optvar.type = [];

optctr.size = 0;
optctr.lb   = [];
optctr.up   = [];
optctr.type = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example - Thermal dispatch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nthermal = 3;
nthermal = 4;

% Capacities [MW]
max_cap = zeros(nthermal,1);
max_cap(1) = 5.0;
max_cap(2) = 10.0;
max_cap(3) = 18.0;
max_cap(4) = 11.0;

% Operation costs (k$/MW)
thermal_ope_costs = zeros(nthermal,1);
thermal_ope_costs(1) = 10.0;
thermal_ope_costs(2) = 15.0;
thermal_ope_costs(3) = 30.0;
thermal_ope_costs(4) = 8.0;

% Demand [MW]
nblocks = 1;
electrical_demand = zeros(nblocks,1);
electrical_demand(1) = 12.0;

var_thermal_pointer = optvar.size;
for i = 1:nthermal
 optvar = add_var(0,max_cap(i),thermal_ope_costs(i),VAR_TYPE_CONT,optvar);
end

ctr_demand_balance_pointer = optctr.size;
for i = 1:nblocks
  optctr = add_ctr(electrical_demand(i), ...
                   electrical_demand(i), ...
                   CTR_TYPE_EQUAL,       ...
                   optctr);
end
optctr.var = zeros(optctr.size,optvar.size);

fprintf("=======================================================================\n");
fprintf("Variables' structure attributes: \n");
fprintf("=======================================================================\n");
fprintf("Size: %d\n",optvar.size);
fprintf("Types:\n")
disp(optvar.type);
fprintf("Lower bounds:\n")
disp(optvar.lb);
fprintf("Upper bounds:\n")
disp(optvar.ub);
fprintf("Objective function coefficients:\n")
disp(optvar.obj);
disp("");

fprintf("=======================================================================\n");
fprintf("Constraints' structure attributes: \n");
fprintf("=======================================================================\n");
fprintf("Size: %d\n",optctr.size);
fprintf("Lower bounds:\n")
disp(optctr.lb);
fprintf("Upper bounds:\n")
disp(optctr.ub);

fprintf("Types:\n")
disp(optctr.type);
disp("");

% Associates variables to each constraint
for i = 1:nblocks
  demand_ctr_index = ctr_demand_balance_pointer + i;
  for j = 1:nthermal
    coef = 1.0;
    thermal_var_index = var_thermal_pointer + j;
    optctr = add_var_to_ctr(coef,thermal_var_index,demand_ctr_index,optctr);
  end
end

fprintf("=======================================================================\n");
fprintf("Constraint parameters: \n")
fprintf("=======================================================================\n");
for i=1:size(optctr.var, 1)
  for j=1:size(optctr.var, 2)
    fprintf("%3.3f  ", optctr.var(i,j))
  end
  fprintf("  | %3.3f \n", optctr.ub(i))
end


[opt_sol, opt_obj, error, info] = solve_glpk_prob(1,optvar,optctr);

% Display results
disp("");
fprintf("=======================================================================\n");
fprintf("Results\n");
fprintf("=======================================================================\n");
if (error==0)
  resp="No";
  fprintf("Error in solution? %s\n", resp)
end
fprintf("\n");

fprintf("Optimal obj. function value: k$ %3.3f\n", opt_obj)
fprintf("Optimal solution:\n")
for i=1:nthermal
  fprintf("Output of plant %d: %3.3f MWh\n", i, opt_sol(i))
end
fprintf("\n");

fprintf("--------------------------------------------------------------------------\n")
fprintf("Sensibility analysis\n")
fprintf("--------------------------------------------------------------------------\n")
disp("");

fprintf("Marginal costs: \n")
for i=1:length(info.lambda)
  fprintf("Restrição %d: %3.3f\n", [i, info.lambda(i)])
end
disp("");

fprintf("Reduced costs (variable marginal costs): \n")
for i=1:length(info.redcosts)
  fprintf("Variable %d: %3.3f\n", [i, info.redcosts(i)])
end
