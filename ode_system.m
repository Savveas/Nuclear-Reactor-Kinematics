function dY = ode_system(t, Y)
global rho rhoi a_f a_c T_f T_c T_i M_f M_c R_f W l_p tau c_f c_p beta lambda b_d T_f0 T_c0
    dY=zeros(4,1);
    P = Y(1);
    c = Y(2);
    T_f = Y(3);
    T_c = Y(4);
    rhoi0=rhoi;
    rho=rhoi0+a_f*(Y(3)-T_f0)+a_c*(Y(4)-T_c0);
    dY(1) = (rho-beta)*Y(1)/l_p+lambda*Y(2);
    dY(2) = beta * Y(1) / l_p - lambda*Y(2);
    dY(3) = 1/(M_f * c_f)*Y(1) - 1/tau*(Y(3)-Y(4));
    dY(4) = 1/(M_c*c_p*R_f)*(Y(3)-Y(4))-2*W/M_c*(Y(4)-T_i);

return
