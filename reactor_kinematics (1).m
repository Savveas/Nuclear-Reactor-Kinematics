% Initial contitions
clear;
clc;
global rho rhoi a_f a_c T_i M_f M_c R_f W l_p tau c_f c_p beta lambda b_d T_f0 T_c0
k_f=2.5;
k_c=13;
c_p=5.42*10^(3);
h=30*10^(3);
N=50952;
c_f=0.247*(10^3);
a_f=-2*10^(-5);
a_c=-4*10^(-5);
H=3.7;
niSigma_f=15.7;
u=2200;
l_p=1/(niSigma_f*u);
T_i=293.7+273.15;
W=17400;
K_c=14;
a_UO2=8.2*10^(-3);
b_d=0.57*10^(-3);
M_f=101000;
M_c=750*8.79*10^(-5)*50952*3.7;
beta=0.0065;
a=4.1*10^(-3);
A=8.79*10^(-5);
lambda=1/13;
l_p=1/[niSigma_f*u];
R_f =(1/(3.14)*N*H*4.0*k_f)+(log((a+b_d)/b_d)/(3.14*N*H*2.0*K_c))+(1/3.14*N*(H^2)*2*(a+b_d));
tau=R_f*M_f*c_f;


    % Initial conditions
    P0 = 3*10^(6);
    c0 = (beta * P0) / (lambda * l_p);
    rhoi=0.1*beta;
    T_c0 = (P0 / (2 * W * c_p)) + T_i;
    T_f0=T_c0+R_f*P0;
    Y0=[P0,c0,T_f0,T_c0];
    rho=rhoi;
    % Time span
    tspan = [0 100]; % Adjust the end time as needed

    % Solve the system of ODEs
    [t,Y] = ode15s(@ode_system_fun, tspan, Y0);

    rho=rhoi+a_f*(Y(3)-T_f0)+a_c*(Y(4)-T_c0);

    % Extracting the results
    P = Y(:,1);
    c = Y(:,2);
    T_f = Y(:,3);
    T_c = Y(:,4);

 % Plot the results if needed
    figure;

    subplot(2, 2, 1);
    plot(t, P);
    xlabel('Time');
    ylabel('P');
    title('P vs Time');

    subplot(2, 2, 2);
    plot(t, c);
    xlabel('Time');
    ylabel('c');
    title('c vs Time');

    subplot(2, 2, 3);
    plot(t, T_f);
    xlabel('Time');
    ylabel('T_f');
    title('T_f vs Time');

    subplot(2, 2, 4);
    plot(t, T_c);
    xlabel('Time');
    ylabel('T_c');
    title('T_c vs Time');




