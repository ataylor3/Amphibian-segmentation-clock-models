%% FIGURES S6 and S7
clear 
clc

% global parameters
tfinal = 3100; % total time span (min)

a = 4.5; % protein synthesis rate
k = 33; % mRNA synthesis rate in absence of inhibition
p_crit = 420; % critical protein threshold (species-specific)

Ttx = 6.89; % transcriptional delay
Tin = 12.87; % intron splicing delay
Tp = 2.18; % translation delay


%% S6A and S7A: A. mexicanum, Brownian Motion (normal diffusion) model
% mRNA half-life = export time

Texp = 11.97; % export delay (BM simulations, radius = 5.5)
Tm_base = Ttx + Tin + Texp - 5; 

HL_m = Texp; % mRNA half life equal to export delay
HL_p = 3; % base protein half-life

for i = 1:20

    Tm = Tm_base;

    for j = 1:20

        soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

        t = soln.x;
        counts = soln.y;
        protein = counts(1,:);
        mRNA = counts(2,:);

        osc = osc_behavior(t,tfinal,protein,mRNA); 

        P_mRNA(i,j) = osc(1,1); % store period of gene expression
        A_mRNA(i,j) = osc(1,3); % store amplitude of mRNA expression
        A_protein(i,j) = osc(1,4); % store amplitude of protein expression

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay; % store total delay

        Tm = Tm + 0.5; % add on ~10 min to get +/- 5 min range

    end


    HL_p_vector(i,1) = HL_p; % store protein half-life
    HL_p = HL_p + 1;


end

% plot amplitude of mRNA expression 
A_mRNA(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped 

figure()
imagesc(Total_delay,HL_p_vector,A_mRNA)
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of mRNA expression (molecules)'
c.Label.FontSize = 15;

% plot amplitude of protein expression 
A_protein(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped 

figure()
imagesc(Total_delay,HL_p_vector,A_protein)
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of protein expression (molecules)'
c.Label.FontSize = 15;

%% S6C and S7C: A. mexicanum, Brownian Motion (normal diffusion) model
% mRNA half-life = 1/2 export time

HL_m = Texp/2; % mRNA half life equal to 1/2 export delay
HL_p = 3;% base protein half-life

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
        A_protein(i,j) = osc(1,4);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5; % add on ~10 min to get +/- 5 min range

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot amplitude of mRNA expression 
A_mRNA(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_mRNA)
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = 1/2 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of mRNA expression (molecules)'
c.Label.FontSize = 15;

% plot amplitude of protein expression 
A_protein(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_protein)
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = 1/2 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of protein expression (molecules)'
c.Label.FontSize = 15;


%% S6E and S7E: A. mexicanum, Brownian Motion (normal diffusion) model
% mRNA half-life = 1/4 export time

HL_m = Texp/4; % mRNA half life equal to 1/4 export delay
HL_p = 3; 

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
        A_protein(i,j) = osc(1,4);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5; % add on ~10 min to get +/- 5 min range

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot amplitude of mRNA expression 
A_mRNA(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_mRNA)
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = 1/4 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of mRNA expression (molecules)'
c.Label.FontSize = 15;

% plot amplitude of protein expression 
A_protein(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_protein)
c = colorbar;
title('{\it A. mexicanum} (BM model)', '{h_m = 1/4 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of protein expression (molecules)'
c.Label.FontSize = 15;

%% S6B and S7B: A. mexicanum, fractional Brownian Motion (obstructed diffusion) model
% mRNA half-life = export time

Texp = 26.27; % export delay (fBM simulations, radius = 5.5)

Tm_base = Ttx + Tin + Texp - 5;

HL_m = Texp; % mRNA half-life equal to export delay
HL_p = 3; 

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
        A_protein(i,j) = osc(1,4);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5;

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot amplitude of mRNA expression 
A_mRNA(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_mRNA)
c = colorbar;
title('{\it A. mexicanum} (fBM model)', '{h_m = T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of mRNA expression (molecules)'
c.Label.FontSize = 15;


% plot amplitude of protein expression 
A_protein(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_protein)
c = colorbar;
title('{\it A. mexicanum} (fBM model)', '{h_m = T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of protein expression (molecules)'
c.Label.FontSize = 15;

%% S6D and S7D: A. mexicanum, fractional Brownian Motion (obstructed diffusion) model
% mRNA half-life = 1/2 export time

HL_m = Texp/2; % mRNA half-life equal to 1/2 export delay
HL_p = 3; 

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
        A_protein(i,j) = osc(1,4);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5;

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot amplitude of mRNA expression 
A_mRNA(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_mRNA)
c = colorbar;
title('{\it A. mexicanum} (fBM model)', '{h_m = 1/2 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of mRNA expression (molecules)'
c.Label.FontSize = 15;

% plot amplitude of protein expression 
A_protein(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_protein)
c = colorbar;
title('{\it A. mexicanum} (fBM model)', '{h_m = 1/2 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of protein expression (molecules)'
c.Label.FontSize = 15;

%% S6F and S7F: A. mexicanum, fractional Brownian Motion (obstructed diffusion) model
% mRNA half-life = 1/4 export time

HL_m = Texp/4; % mRNA half-life equal to 1/4 export delay
HL_p = 3; 

p_crit = 420;

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
        A_protein(i,j) = osc(1,4);

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay;

        Tm = Tm + 0.5;

    end


    HL_p_vector(i,1) = HL_p;
    HL_p = HL_p + 1;


end


% plot amplitude of mRNA expression 
A_mRNA(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_mRNA)
c = colorbar;
title('{\it A. mexicanum} (fBM model)', '{h_m = 1/4 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of mRNA expression (molecules)'
c.Label.FontSize = 15;

% plot amplitude of protein expression 
A_protein(isinf(P_mRNA)) = 0; % set amplitude to 0 if oscillations are damped

figure()
imagesc(Total_delay,HL_p_vector,A_protein)
c = colorbar;
title('{\it A. mexicanum} (fBM model)', '{h_m = 1/4 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Amplitude of protein expression (molecules)'
c.Label.FontSize = 15;
