%% TEST THE IMPACT OF FURTHER INCREASING GENE STABILITY

% A. mexicanum Brownian Motion model

% Can a period of 154 minutes be achieved by further increasing gene
% product stability (without any increase in total delay time)?

% Take A. mexicanum BM model that is closest to achieving a 154 minute
% period of gene expression: 
% A. mexicanum Brownian Motion model, mRNA half-life hl_m = Texp (BM)

% Start by clearing workspace and command window
clc
clear

% Set global parameters 
tfinal = 3100; 
a = 4.5; % protein synthesis rate
k = 33; % mRNA synthesis rate in absence of inhibition
p_crit = 420; % critical protein threshold (species-specific)

%% 

% Increase protein stability: start by changing the range of protein stability 
% to a higher range (i.e. change from 3 to 23 minutes to 15 to 35 minutes

% Delay parameters stay the same
Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 11.97; % export delay (BM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5; % take off 5 min to get +/- 5 minute range
Tp = Ttl;

HL_m = Texp; % mRNA half life equal to export delay (BM model)
HL_p = 15; % starting protein half life 


HL_m_vector = [];
P_mRNA = [];

for i = 1:20

    Tm = Tm_base;

    for j = 1:20

        soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

        t = soln.x;
        counts = soln.y;
        protein = counts(1,:);
        mRNA = counts(2,:);

        osc = osc_behavior(t,tfinal,protein,mRNA); 

        P_mRNA(i,j) = osc(1,1);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5; % add on ~10 min to get +/- 5 min range

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot resulting heat map
figure()
imAlpha=ones(size(P_mRNA));
imAlpha(isinf(P_mRNA))=0;
imagesc(Total_delay,HL_p_vector, P_mRNA,'AlphaData',imAlpha);
set(gca,'color',[0.15, 0.07, 0.25]);
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;
%%

% Additional increase in mRNA stability (increase by 25%)

HL_m = Texp * 1.25; % 25% increase, all other parameters are the same
HL_p = 15; % starting protein half life (still looking at 15 to 35 minute range)

for i = 1:20

    Tm = Tm_base;

    for j = 1:20

        soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

        t = soln.x;
        counts = soln.y;
        protein = counts(1,:);
        mRNA = counts(2,:);

        osc = osc_behavior(t,tfinal,protein,mRNA); 

        P_mRNA(i,j) = osc(1,1);
        A_mRNA(i,j) = osc(1,3);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5; % add on ~10 min to get +/- 5 min range

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot resulting heat map
figure()
imAlpha=ones(size(P_mRNA));
imAlpha(isinf(P_mRNA))=0;
imagesc(Total_delay,HL_p_vector, P_mRNA,'AlphaData',imAlpha);
set(gca,'color',[0.15, 0.07, 0.25]);
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = 14.96}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% Additional increase in mRNA stability (increase by 50%, relative to Texp)

HL_m = Texp * 1.50; % 50% increase, all other parameters are the same
HL_p = 15; % starting protein half life

for i = 1:20

    Tm = Tm_base;

    for j = 1:20

        soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

        t = soln.x;
        counts = soln.y;
        protein = counts(1,:);
        mRNA = counts(2,:);

        osc = osc_behavior(t,tfinal,protein,mRNA); 

        P_mRNA(i,j) = osc(1,1);
        A_mRNA(i,j) = osc(1,3);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5; % add on ~10 min to get +/- 5 min range

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot resulting heat map
figure()
imAlpha=ones(size(P_mRNA));
imAlpha(isinf(P_mRNA))=0;
imagesc(Total_delay,HL_p_vector, P_mRNA,'AlphaData',imAlpha);
set(gca,'color',[0.15, 0.07, 0.25]);
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = 17.96}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;



