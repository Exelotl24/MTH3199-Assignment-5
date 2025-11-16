clear
close all

main()

function main()

    % Create BT struct for Forward Euler
    BT_struct_FE = struct();
    BT_struct_FE.C = [0];
    BT_struct_FE.B = [1];
    BT_struct_FE.A = [0];

    % Create BT struct for Explicit Midpoint
    BT_struct_EM = struct();
    BT_struct_EM.C = [0; 0.5];
    BT_struct_EM.B = [0, 1];
    BT_struct_EM.A = [0, 0; 0.5, 0];

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

    % Run animation
    spring_plotting_example()


end