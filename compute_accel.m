%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of static mounting
% points for the spring (in the world frame)
%box_params.P_box: 2 x n list of mounting points
% for the spring (in the box frame)
%
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)

    m = box_params.m;
    I = box_params.I;
    g = box_params.g;
    
    % Centroid position
    PC = [x; y];

    num_springs = size(box_params.P_world, 2); 
    
    % Net forces
    F_net = [0; -m * g]; 
    tau_net = 0;   

    Plist_world = compute_rbt(x,y,theta,box_params.P_box);

    for i = 1:num_springs

        % Spring params
        k_i = box_params.k_list(i);
        l0_i = box_params.l0_list(i);

        % Static mounting point
        PA = box_params.P_world(:,i);
        % Attachment point on box
        PB = Plist_world(:,i);

        % Calculate force by spring on point B
        F_i = compute_spring_force(k_i,l0_i,PA,PB);

        % Update F_net
        F_net = F_net + F_i;

        r_i = PB - PC;

        r_x = r_i(1);
        r_y = r_i(2);

        F_x = F_i(1);
        F_y = F_i(2);

        tau_i = (r_x * F_y) - (r_y * F_x);

        % Update tau_net
        tau_net = tau_net + tau_i;
    end

    % Linear acceleration
    ax = F_net(1) / m;
    ay = F_net(2) / m;

    % Angular acceleration
    atheta = tau_net / I;

end