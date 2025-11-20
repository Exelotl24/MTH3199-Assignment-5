function simulate_box(V0, tspan, h_ref, box_params, BT_Struct)

    

    % Start Figure
    figure()
    hold on
    axis equal
    axis([-10, 10, -10, 10])
    
    % Box plot
    h_box_plot = plot(0, 0, 'b-', 'linewidth', 2);
    
    % Spring plots
    num_zigs = 5;
    w = 0.1;
    num_springs = size(box_params.P_world, 2);
    spring_plots = cell(1, num_springs);

    for i = 1:num_springs
        spring_plots{i} = initialize_spring_plot(num_zigs, w);
    end

    %load the system parameters into the rate function
    %via an anonymous function
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

    
    %run the integration
    [~, Vlist, ~, ~] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, DormandPrince);

    % Animation loop
   
    for k = 1:size(Vlist, 1)
        V_i = Vlist(k,:)';
        x = V_i(1);
        y = V_i(2);
        theta = V_i(3);

        % Update box
        box_world = compute_rbt(x, y, theta, box_params.boundary_pts);
        set(h_box_plot, 'XData', box_world(1,:), 'YData', box_world(2,:));

        % Update spring
        for i = 1:num_springs
            PA = box_params.P_world(:, i);
            PB = compute_rbt(x, y, theta, box_params.P_box(:, i));
            update_spring_plot(spring_plots{i}, PA, PB);
        end
  
        drawnow;
    end

end