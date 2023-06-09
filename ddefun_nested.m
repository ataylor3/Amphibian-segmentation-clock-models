function dydt = ddefun_nested(Tm, Tp, HL_m, HL_p,a,k,p_crit) % DDE SOLVER
    
    tspan = [0 3100];
    
    delays = [Tm Tp];
    
    b = log(2)/HL_p; % rate degradation of protein
    c = log(2)/HL_m; %rate degradation of mRNA
    
    dydt = ddesd(@ddefun, delays, @history, tspan);
    
    function dydt = ddefun(t,y,Z)
        dydt = [a*Z(2,2) - b*y(1); k/(1+(Z(1,1)/p_crit)^2) - c*y(2)];
    end
end