%% Sensitivity of period and amplitude of mRNA expression to various parameter changes 


% Start by clearing workspace and command window
clear
clc

% Set global parameter
tfinal = 3100; % time span

%% % A. mexicanum BM model (hlm = 3), shown in main text

% set delay parameters 
Ttx = 6.89; % transcriptional delay
Tin = 12.78; % intron splicing delay
Texp = 11.97; % export delay (simulat,ions)
Ttl = 2.18; % translation delay

% Testing the impact of a 50% increase and 50% decrease in each parameter
delta = [1.5, 0.5];

for i = 1:2 % outer loop: 1st run --> 50% inc.; 2nd run --> 50% dec.
    for j = 1:7 % inner loop: changes each parameter individually
        
        % original parameter values corresponding to model
        a = 4.5; 
        k = 33;
        p_crit = 420;
        HL_m = 3; 
        HL_p = 15; % choose a protein half-life value in the middle of the 
        % range tested in model (3 to 23 minutes)
        Tm = Ttx + Tin + Texp;
        Tp = Ttl;
        Tm + Tp
    
        V = [a,k,p_crit,HL_m,HL_p,Tm,Tp]; % vector of parameters

        if j == 1 % on the first run, solve system with original set of parameters
            soln = ddefun_nested(Tm,Tp,HL_m,HL_p,a,k,p_crit);

        elseif j == 7 % on the last run, address Tm and Tp at the same time (change total delay)
            V(6:7) = V(6:7) * delta(i); % 50% inc/dec in total delay
            soln = ddefun_nested(V(6),V(7),V(4),V(5),V(1),V(2),V(3)); % solve with updated set of parameters 

        else % otherwise, change each parameter in vector V individually 
            V(j-1) = V(j-1) * delta(i); % 50% inc/dec in each parameter
            soln = ddefun_nested(V(6),V(7),V(4),V(5),V(1),V(2),V(3)); % solve with updated set of parameters
        end 
    
        % system solution
        t = soln.x;
        counts = soln.y;
        protein = counts(1,:);
        mRNA = counts(2,:);
    
        osc = osc_behavior(t,tfinal,protein,mRNA);
    
        P_mRNA(j,i) = osc(1,1); % store period of gene expression
        A_mRNA(j,i) = osc(1,3); % store amplitude of gene expression
    
            
    end
end

% Create tables to track the impact of parameter changes on period and
% amplitude of gene expression

parameter_changes = ["No change","a","k","p_crit","h_m","h_p","T"];
column_names = ["Parameter changes","Period of mRNA expression","Amplitude of mRNA expression"];

inc_50 = table(parameter_changes', P_mRNA(:,1), A_mRNA(:,1), 'VariableNames', column_names)
dec_50 = table(parameter_changes', P_mRNA(:,2), A_mRNA(:,2), 'VariableNames', column_names)
