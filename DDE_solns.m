%% GENERAL METHOD OF SOLVING DDE SYSTEMS 

% USE THIS TO SHOW MINIMAL DIFFERENCES BETWEEN MRNA AND PROTEIN PERIODS OF GENE EXPRESSION
% (Supplemental Material 1)

% Start by clearing workspace and command window
clear
clc

% Global parameters 
tfinal = 3100; % always set the time span to 3,100 minutes
a = 4.5; % rate of protein synthesis (protein/mRNA/min)
k = 33; % mRNA synthesized in absence of inhibition (mRNA/cell/min)

HL_m = 3; % mRNA half life
HL_p = 3; % protein half life

%%

% X. laevis parameters (Brownian Motion model hm, hp = 3 min)

Tp = 1.29; % protein production delay (translation delay)

Ttx = 1.34; % transcriptional delay
Tin = 4.15; % intron delay
Texp = 6.39; % export delay

Tm = Ttx + Tin + Texp; % mRNA production delay

p_crit = 161; % critical protein threshold

% Solution to DDE system

% use local function ddefun_nested to generate a solution to the DDE system
% given the model parameters plugged in 
soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

t = soln.x; % gives a time vector
counts = soln.y; % gives a matrix of mRNA and protein levels
% split counts matrix into specific mRNA and protein vectors 
protein = counts(1,:); % corresponding vector of mRNA levels over time
mRNA = counts(2,:); % corresponding vector of protein levels over time

% plug solution into local function osc_behavior that assesses period and
% amplitude of gene expression
osc = osc_behavior(t,tfinal,protein,mRNA); 

% we can retrieve the period of mRNA and protein expression to compare
P_mRNA = osc(1,1) % period of mRNA expression
P_protein = osc(1,2) % period of protein expression

%%

% A. mexicanum parameters (Brownian Motion model, hm, hp = 3 min)

Tp = 2.18; % protein production delay (translation delay)

Ttx = 6.89; % transcriptional delay
Tin = 12.87; % intron delay
Texp = 11.97; % export delay

Tm = Ttx + Tin + Texp; % mRNA production delay

p_crit = 420; % critical protein threshold

% Solution to DDE system

soln = ddefun_nested(Tm, Tp, HL_m, HL_p, a, k, p_crit);

t = soln.x;
counts = soln.y;
protein = counts(1,:);
mRNA = counts(2,:);

osc = osc_behavior(t,tfinal,protein,mRNA); 
P_mRNA = osc(1,1) % period of mRNA expression
P_protein = osc(1,2) % period of protein expression

%% FUNCTIONS CALLED

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


% USES LOCAL EXTREMA TO FIND PERIOD AND AMPLITUDE OF EXPRESSION
function osc_behavior = f(t,tfinal,protein,mRNA)

    % look for local extrema
    % islocal functions give logical array, 1s at extrema, 0s otherwise

    p_max = protein(islocalmax(protein)); % vector of local protein max
    p_min = protein(islocalmin(protein)); % vector of local protein min
    m_max = mRNA(islocalmax(mRNA)); % vector of local mRNA max
    m_min = mRNA(islocalmin(mRNA)); % vector of local mRNA min
    
    t_p_min = t(islocalmin(protein)); % vector of time stamps, local p min
    t_m_min = t(islocalmin(mRNA)); % vector of time stamps, local m min
    
    % cut out first 5 cycles of oscillation to avoid skewing data
    p_max = p_max(1,6:end);
    p_min = p_min(1,5:end);
    m_max = m_max(1,6:end);
    m_min = m_min(1,5:end);
    
    t_p_min = t_p_min(1,5:end);
    t_m_min = t_m_min(1,5:end);

    % find period of expression
    periodicity_mRNA = P(t_m_min,tfinal,m_min,m_max);
    periodicity_protein = P(t_p_min,tfinal,p_min,p_max);
    
    % find amplitude of expression
    amplitude_mRNA = Amp(m_max,m_min);
    amplitude_protein = Amp(p_max,p_min);
    
    % store results
    osc_behavior = [periodicity_mRNA, periodicity_protein, 
        amplitude_mRNA, amplitdue_protein];
end


function P = p(min_time_vector,tfinal,min_mol,max_mol) % PERIOD OF EXPRESSION

    P = mean(diff(min_time_vector)); % period of expression is the average number of minutes 
    % between successive local minima 

    % if the amplitude of expression < 10 for any cycle of expression, set
    % period to infinity (i.e. damped oscillations)

    % set number of complete cycles of oscillation to s
    if size(max_mol,2) == size(min_mol,2) | size(min_mol,2) > size(max_mol,2)
        s = size(max_mol,2);
    elseif size(max_mol,2) > size(min_mol,2)
        s = size(min_mol,2);
    end
    
    % go through each cycle of oscillation to check for amplitude < 10 
    for i = 1:s
        if abs(max_mol(1,i) - min_mol(1,i)) < 10
            P = inf;
            break
        end
    end
end


function Amp = amp(max_mol,min_mol) % AMPLITUDE OF EXPRESSION

    % make sure corresponding maxima and minima are matched up
    if size(max_mol,2) == size(min_mol,2)
        amp = max_mol - min_mol;
    elseif size(max_mol,2) > size(min_mol,2)
        amp = max_mol(:,1:size(min_mol,2)) - min_mol;
    elseif size(min_mol,2) > size(max_mol,2)
        amp = max_mol - min_mol(:,1:size(max_mol,2));
    end
    
    
    Amp = mean(amp);

end


