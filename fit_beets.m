function yr = fit_beets(p,t)   
% Function containing two coupled ODEs. Solves the ODES at timepoints t
% with parameter vector p
    x0=[abs(p(3)), 0]; 
    [t,m]=ode45(@beet,t,x0);

    % Function defining the ODE
    function dydt = beet(t,m)
        dydt = zeros(2,1);
        dydt(1) = p(1)*(800-m(1))*m(1);
        dydt(2)= p(2)*(150-m(2))*m(1);
    end
    yr = m(:,2);
end