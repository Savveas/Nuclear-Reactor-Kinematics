function dY = ode_system(t, Y)

global ro ro0 af ac Tf0 Tc0 Rf tau1 tau2 W Ti Mc Mf cf cp b bd lp lamda
    dY = zeros(4,1);


    ro=ro0+af*(Y(3)-Tf0)+ac*(Y(4)-Tc0);         ### reactivity
    dY(1) = (ro - bd) * Y(1)/ lp + lamda * Y(2);##### isxus
    dY(2) = bd * Y(1)/ lp - lamda * Y(2);       ### c precurssor consentration
    dY(3) = Y(1) / (Mf * cp) - (Y(3) - Y(4)) / tau1;    ####  Tf
    dY(4) = (Y(3) - Y(4)) / tau2 - 2 * W * (Y(4) - Ti) / Mc;  ### Tc


return

