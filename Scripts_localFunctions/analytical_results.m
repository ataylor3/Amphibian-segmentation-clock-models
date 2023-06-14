%% ANALYTICAL RESULTS: PLOT TCRIT AGAINST PROTEIN HALF-LIFE 

% start by clearing workspace and command terminal
clear
clc

% Model specific parameters
p_crit_XL = 161;
p_crit_AM = 420;

initial_hl_m = 3;

Texp_amex_BM = 11.97;
Texp_amex_fBM = 26.27;
%%

% Global parameters (some based on model specific parameters defined above)

k = 33;
a = 4.5;
hl_p = 3:0.5:23; % Fixed range of protein half-life

hl_m = Texp_amex_fBM/4; % mRNA half-life (model specific, adjust as needed)
p_crit = p_crit_AM; % pcrit (model specific, adjust as needed) 

c = log(2)/hl_m; % degradation of mRNA
blist = log(2)./hl_p; % degradation of protein

y_o = p_crit/(k*a);

%% 

%Steady state y* (y_s below)

Nos = length(blist);
for i=1:Nos
    b=blist(i);
first_term = -((3*2^(1/3)*y_o^2)/ ...
    ((((27*y_o^2)/(b*c)) + ...
    3*sqrt(((81*y_o^4)/(b^2*c^2)) + 12*y_o^6))^(1/3)));

second_term = ((((27*y_o^2)/(b*c)) + ...
    3*sqrt(((81*y_o^4)/(b^2*c^2)) + 12*y_o^6))^(1/3))/(2^(1/3));

y_s = (1/3)*(first_term + second_term);

beta = (y_s/y_o)^2;

K = (2*beta)/(y_s*(1+beta)^2);
if K > b*c

    w = sqrt((-(c^2 + b^2) + sqrt((c^2 + b^2)^2 + 4*K^2 - 4*b^2*c^2))/2);

    T_crit(i) = (1/w)*asin((w*(c+b))/K);
end

end

%% 

% Plot results

plot(hl_p, T_crit)
xlabel("Protein half-life, hl_p (min)",'FontSize',15)
ylabel("T_{crit}",'FontSize',15)
title("T_{crit} across increasing protein stability","{\it A. mexicanum}, h_m = 1/4 T_{exp} (fBM)", 'FontSize',18)
xlim([hl_p(1) hl_p(end)])
ax = gca;
ax.FontSize = 15;