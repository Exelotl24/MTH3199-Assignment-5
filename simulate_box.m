function simulate_box()

    %define system parameters
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(P_world,2));
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;
    
    %load the system parameters into the rate function
    %via an anonymous function
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    
    x0 = 0.5;
    y0 = 0.5;
    theta0 = 0.1; %radians
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;
    
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    
    tspan = [0, 10];
    dt = 0.01;
    
    %run the integration
    [tlist, Vlist, h_avg, num_evals] = explicit_RK_fixed_step_integration(box_rate_func, tspan, V0, h_ref, BT_struct);

end