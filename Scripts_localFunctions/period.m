%% ASSESSING PERIOD OF GENE EXPRESSION (Figure 1)

% start by clearing workspace and command window
clear
clc

% Global parameters
tfinal = 3100; % global parameter
a = 4.5; % rate of protein synthesis (protein/mRNA/min)
k = 33; % mRNA synthesis rate in absence of inhibition (mRNA/cell/min)

% Initial set of species- and diffusion-specific models:

%%

% Xenopus laevis, Brownian Motion (normal diffusion) model

Ttx = 1.34; % transcription delay
Tin = 4.65; % intron delay
Texp = 6.39; % export delay (BM simulations, radius = 4)
Ttl = 1.29; % translation delay

Tm_base = Ttx + Tin + Texp - 5; % mRNA synthesis delay is a sum of 
% transcription, intron, and export delay
% start 5 min under estimated total delay to gt

Tp = Ttl; % protein synthesis delay is equal to translation delay

HL_m = 3; % mRNA half life
HL_p = 3; % protein half life lower bound

p_crit = 161; % critical protein threshold (species-specific value)

% set up vectors and matrices for storing data
HL_p_vector = [];
P_mRNA = [];


for i = 1:20 % use outer loop for increasing protein half-life

    Tm = Tm_base;

    for j = 1:20 % use inner loop for increasing total delay time

        % solve the DDE system associated with each protein stability/
        % delay parameter combination
        soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

        t = soln.x;
        counts = soln.y;
        protein = counts(1,:);
        mRNA = counts(2,:);

        % use local function to assess oscillatory behavior of solution
        % (i.e. period and amplitude of expression)
        osc = osc_behavior(t,tfinal,protein,mRNA); 

        P_mRNA(i,j) = osc(1,1); % store period for each combination of protein half-life and total delay

        tot_delay = Tm + Tp;
        Total_delay(j,1) = tot_delay; % store total delay time

        Tm = Tm + 0.5; % increase by 10 min to get range +/- 5 min within estimated total delay time
        % Note: we increase only Tm because we are summing over Tm and Tp
        % to get total delay; individual delays contribute to oscillatory
        % behavior through their impacts on total delay 

    end


    HL_p_vector(i,1) = HL_p; % store protein half-life
    HL_p = HL_p + 1; % add 1 minute onto protein half-life


end


% plot resulting heat map
figure()
imAlpha=ones(size(P_mRNA));
imAlpha(isinf(P_mRNA))=0;
imagesc(Total_delay,HL_p_vector, P_mRNA,'AlphaData',imAlpha);
set(gca,'color',[0.15, 0.07, 0.25]);
c = colorbar;
title('{\it X. laevis}', '(BM model)','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%% 

% Ambystoma mexicanum, Brownian Motion (normal diffusion) model

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron delay
Texp = 11.97; % export delay (BM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5;
Tp = Ttl;

HL_m = 3; 
HL_p = 3; 

p_crit = 420;

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

        Tm = Tm + 0.5; 

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
title('{\it A. mexicanum}', '(BM model)','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% Xenopus laevis, fractional Brownian Motion (obstructed diffusion) model

Ttx = 1.34; % transcriptional delay
Tin = 4.65; % intron delay
Texp = 8.36; % export delay (fBM simulations, radius = 4)
Ttl = 1.29; % translation delay

Tm_base = Ttx + Tin + Texp - 5;
Tp = Ttl;

HL_m = 3; 
HL_p = 3; 

p_crit = 161;

HL_p_vector = [];
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

        Tm = Tm + 0.5; 

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
title('{\it X. laevis}', '(fBM model)','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% Ambystoma mexicanum, fractional Brownian Motion (obstructed diffusion)
% model

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron delay
Texp = 26.27; % export delay (fBM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5;
Tp = Ttl;

HL_m = 3; 
HL_p = 3; 

p_crit = 420;

HL_p_vector = [];
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

        Tm = Tm + 0.5; 

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
title('{\it A. mexicanum}', '(fBM model)','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%% TESTING OUT MRNA STABILITY SCENARIOS (Figure 2)

% A. mexicanum, Brownian Motion (normal diffusion) model
% mRNA half-life = export time

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 11.97; % export delay (BM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5; 
Tp = Ttl;

HL_m = Texp; % mRNA half life equal to export delay
HL_p = 3; 

p_crit = 420; 

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

% A. mexicanum, Brownian Motion (normal diffusion) model
% mRNA half-life = 1/2 export time

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 11.97; % export delay (BM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5; 
Tp = Ttl;

HL_m = Texp/2; % mRNA half life equal to 1/2 export delay
HL_p = 3; 

p_crit = 420; 

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
title('{\it A. mexicanum} (BM model)', '{h_m = 1/2 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% A. mexicanum, Brownian Motion (normal diffusion) model
% mRNA half-life = 1/4 export time

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 11.97; % export delay (BM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5; 
Tp = Ttl;

HL_m = Texp/4; % mRNA half life equal to 1/2 export delay
HL_p = 3; 

p_crit = 420; 

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
title('{\it A. mexicanum} (BM model)', '{h_m = 1/4 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% A. mexicanum, fractional Brownian Motion (obstructed diffusion) model
% mRNA half-life = export time

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 26.27; % export delay (fBM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5;
Tp = Ttl;

HL_m = Texp; % mRNA half-life equal to export delay
HL_p = 3; 

p_crit = 420;

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

        Tm = Tm + 0.5;

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
title('{\it A. mexicanum} (fBM model)', '{h_m = T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% A. mexicanum, fractional Brownian Motion (obstructed diffusion) model
% mRNA half-life = 1/2 export time

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 26.27; % export delay (fBM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5;
Tp = Ttl;

HL_m = Texp/2; % mRNA half-life equal to export delay
HL_p = 3; 

p_crit = 420;

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

        Tm = Tm + 0.5;

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
title('{\it A. mexicanum} (fBM model)', '{h_m = 1/2 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

%%

% A. mexicanum, fractional Brownian Motion (obstructed diffusion) model
% mRNA half-life = 1/4 export time

Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 26.27; % export delay (fBM simulations, radius = 5.5)
Ttl = 2.18; % translation delay

Tm_base = Ttx + Tin + Texp - 5;
Tp = Ttl;

HL_m = Texp/4; % mRNA half-life equal to export delay
HL_p = 3; 

p_crit = 420;

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

        Tm = Tm + 0.5;

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
title('{\it A. mexicanum} (fBM model)', '{h_m = 1/4 T_{exp}}','FontSize',18)
xlabel('Total delay (min), {T_m + T_p}','FontSize',15)
ylabel('Protein half-life (min), {h_p}','FontSize',15)
ax = gca;
ax.FontSize = 15;
c.Label.String = 'Period of gene expression (min)'
c.Label.FontSize = 15;

