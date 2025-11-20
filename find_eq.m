% This function finds the equilibrium state of an unforced mass-spring
% system

function V_eq = find_eq(rate_func, V_0)
    
    % solver params indicate diff equation unknown
    solver_params = struct();
    solver_params.numerical_diff = 0;

    % Find root = V_eq
    [V_eq, ~] = multi_newton_solver(rate_func,V_0,solver_params);

    
end