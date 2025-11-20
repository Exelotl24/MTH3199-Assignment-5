% MAIN FUNCTION
% RUN THIS PROGRAM TO SEE ASSIGNMENT 

main

function main() 
    
    % clear workspace & close figures
    close all
    clear



    % define system parameters
    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2+LW^2);

    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];

    l0 = 1.5*LG;

    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;

    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];

    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];

    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];

    P_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];
    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];

    % package system parameters in a struct
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(P_world,2));
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;

        % Define Dormand Prince BT Struct
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0;...
        3/40, 9/40, 0, 0, 0, 0,0;...
        44/45, -56/15, 32/9, 0, 0, 0,0;...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
        35/384, 0, 500/1113, 125/192, -2187/6784,11/84,0];



    % --------------- INITIAL SIMULAITON ------------

    x0 = 0.1; y0 = 0.1; theta0 = 0.1; 
    vx0 = 0; vy0 = 0; vtheta0 = 0;
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];

    tspan = [0, 5];
    h_ref = 0.01; 

    % Run animation
    % simulate_box(V0, tspan, h_ref, box_params, DormandPrince);

    % -------------- EQUILIBRIUM SIMULATION ------------------

    rate_func = @(t, V) box_rate_func(t,V,box_params);
    V_eq_F = find_eq((@(V)rate_func(0, V)), V0);

    % test equilibrium function
    % simulate_box(V_eq_F, tspan, h_ref, box_params, DormandPrince);

    % ---------- LINEARIZATION -----------

    jacobian_eq = approximate_jacobian((@(V)rate_func(0, V)), V_eq_F);
    my_linear_rate = @(t_in, V_in) jacobian_eq*(V_in-V_eq_F);
    V_eq_J = find_eq(@(V_in) my_linear_rate(0, V_in), V0);

    % test equilibrium function
    % simulate_box(V_eq_J, tspan, h_ref, box_params, DormandPrince);



    % -------------- PERTURBATION ------------------


    x0 = 0.1; y0 = 0.1; theta0 = 0.1; 
    vx0 = 0; vy0 = 0; vtheta0 = 0;


    dx0 = 1;
    dy0 = 1;
    dtheta0 = 1;
    vx0 = 1;
    vy0 = 1;
    vtheta0 = 1;

    %small number used to scale initial perturbation
    epsilon = 0.01;
    V0_perturbation = V_eq_F + epsilon*[dx0;dy0;dtheta0;vx0;vy0;vtheta0];


    [t_list_nonlinear,V_list_nonlinear,~, ~] = explicit_RK_fixed_step_integration(rate_func,tspan,V0_perturbation,h_ref,DormandPrince);
    [t_list_linear,V_list_linear,~, ~] = explicit_RK_fixed_step_integration(my_linear_rate,tspan,V0_perturbation,h_ref,DormandPrince);

    figure()
    plot(t_list_nonlinear, vecnorm(V_list_nonlinear'), 'DisplayName', 'nonlinear')
    hold on
    plot(t_list_linear, vecnorm(V_list_linear'), 'DisplayName', 'linear')
    xlabel('time')
    ylabel('V (norm)')
    legend()


    %run the integration of nonlinear system
    % [tlist_nonlinear,Vlist_nonlinear] =...
    % your_integrator(my_rate_func,tspan,V0,...);

    %run the integration of linear system
    % [tlist_linear,Vlist_linear] =...
    % your_integrator(my_linear_rate,tspan,V0,...);

end