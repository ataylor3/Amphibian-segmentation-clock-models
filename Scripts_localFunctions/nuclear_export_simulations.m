% NUCLEAR EXPORT SIMULATIONS

% start by clearing workspace and command window
clear
clc

%%

% Simulations: Brownian Motion, initial position at the origin
N = 500000; % length of position vector/number of total steps
r = 0.5:0.5:13; % nucleas radius ranging from 0.5 to 13 to reflect what has been observed across tree of life

for i = 1:10000 % 10,000 iterations

    % x, y, and z position vectors are cumulative sums of step sizes/increments drawn from a normal distribution
    % scaling parameter of 0.12 chosen such that a radius of 3 corresponds
    % to 202 steps (each step is one second, 202 steps --> 3.36 min)
    % by default, the first entry of the x, y, and z position vectors is
    % equal to the first randomly generated number from randn(1,N)
    % this implies an initial radial position at the origin (the first step
    % is equal to the increment between a randomly generated number and 0)

    x = cumsum(randn(1,N))*0.12; 
    y = cumsum(randn(1,N))*0.12;
    z = cumsum(randn(1,N))*0.12;

    position_R = sqrt(x.^2 + y.^2 + z.^2); % Euclidean formula for distance gives radial position of particle over N steps

    for j = 1:size(r,2) % for every sphere size

        fpt = find(position_R>r(j),1); % find the first step for which the radial position exceeds the sphere radius
        % the number of steps gives the time required for the particle to
        % first exit a nucleus of that size

        if isempty(fpt) == 1 % if the particle never exits the sphere for an iteration
            fpt = N; % set the particle exit time to the maximum number of steps (N = 50,000 steps/seconds --> 833.33 minutes)
        end
    
        FPT_BM_o(i,j) = fpt; % store first exit times in a matrix
        % rows are iterations, columns are the different radii lengths 

    end


end


print = "done with BM, initial position at the origin" % print statement to keep track of code

%%

% Simulations: fraction Brownian Motion, initial position at the origin

H = 0.25; % set Hurst parameter, this corresponds to a sub-scaling exponent alpha = 0.5 
% Sub-diffusion (obstructed diffusion): mean squared displacement scales
% sub-linearlly with time 

for i = 1:10000
    
    % x, y, and z position vectors are generated directly using the
    % fractional Brownian Motion function wfbm
    % scaling parameter of 0.93 chosen such that a radius of 3 corresponds
    % to 202 steps (each step is one second, 202 steps --> 3.36 min)
    % by default, the first entry of the x, y, and z position vectors is 0
    % which gives an initial position at the origin

    x = wfbm(H,N,10)*0.93; 
    y = wfbm(H,N,10)*0.93;
    z = wfbm(H,N,10)*0.93;

    position_R = sqrt(x.^2 + y.^2 + z.^2);
    position_R = position_R(2:N); % we do not want to count the initial 
    % position at the origin as a step, so we remove the first entry

    for j = 1:size(r,2)

        fpt = find(position_R>r(j),1);

        if isempty(fpt) == 1
            fpt = N;
        end
    
        FPT_fBM_o(i,j) = fpt; 
    end


end

print = "done with fBM, initial position at the origin"

%%

% Simulations: Brownian Motion, initial position drawn from a uniform
% distribution within 3/4 of sphere radius of the origin (i.e. exlcuding
% nuclear periphery)

for i = 1:10000

    x = cumsum(randn(1,N))*0.12;
    y = cumsum(randn(1,N))*0.12;
    z = cumsum(randn(1,N))*0.12;

    % generate 3 initial positions (one for each x, y, z position vector)
    % from a uniform distribution across 3/4 of the total sphere radius
    initial_pos = (0.75*r+0.75*r).*rand(3,1) - 0.75*r; 

    for j = 1:size(r,2)
        
        % update position vectors to account for non-zero 0 initial positions
        % the increments/step sizes stay the same, but the position
        % reflects a starting position that is not at the origin
        x_u = x + initial_pos(1,j);
        y_u = y + initial_pos(2,j);
        z_u = z + initial_pos(3,j);
        
        position_R = sqrt(x_u.^2 + y_u.^2 + z_u.^2); % because the first entry
        % is the position after the first step, we do not need to cut it 

        fpt = find(position_R>r(j),1);

        if isempty(fpt) == 1
            fpt = N;
        end
    
        FPT_BM_ud(i,j) = fpt; 

    end


end

print = "done with BM, initial position drawn from uniform distribution"

%%

% Simulations: fractional Brownian Motion, initial position drawn from a uniform
% distribution within 3/4 of sphere radius of the origin (i.e. exlcuding
% nuclear periphery

for i = 1:10000
    
    x = wfbm(H,N,10)*0.93; 
    y = wfbm(H,N,10)*0.93;
    z = wfbm(H,N,10)*0.93;

    initial_pos = (0.75*r+0.75*r).*rand(3,1) - 0.75*r;

    for j = 1:size(r,2)

        x_u = x + initial_pos(1,j);
        y_u = y + initial_pos(2,j);
        z_u = z + initial_pos(3,j);
        
        position_R = sqrt(x_u.^2 + y_u.^2 + z_u.^2); 
        position_R = position_R(2:N); % cut first entry, 
        % initial position should not count as a step

        fpt = find(position_R>r(j),1);

        if isempty(fpt) == 1
            fpt = L;
        end
        
        FPT_fBM_ud(i,j) = fpt;

    end


end

print = "done with fBM, initial position drawn from uniform distribution"

%%

% VISUALIZING SIMULATION RESULTS 

% Generate statistics 
% We would like mean and standard deviation to be in terms of minutes
% each step is one second, divide by 60 to get in terms of minutes

mean_BM_o = mean(FPT_BM_o,1)/60;
mean_BM_ud = mean(FPT_BM_ud,1)/60;
mean_fBM_o = mean(FPT_fBM_o,1)/60;
mean_fBM_ud = mean(FPT_fBM_ud,1)/60;

sd_BM_o = std(FPT_BM_o,0,1)/60;
sd_BM_ud = std(FPT_BM_ud,0,1)/60;
sd_fBM_o = std(FPT_fBM_o,0,1)/60;
sd_fBM_ud = std(FPT_fBM_ud,0,1)/60;

cv_BM_o = sd_BM_o./mean_BM_o;
cv_BM_ud = sd_BM_ud./mean_BM_ud;
cv_fBM_o = sd_fBM_o./mean_fBM_o;
cv_fBM_ud = sd_fBM_ud./mean_fBM_ud;

%%

% Tables
var_names = ["Radius (\mum)","Mean nuclear export time (min)", "Standard deviation (min)","Coefficient of variation"];

results_BM_o = table(r', mean_BM_o', sd_BM_o', cv_BM_o', 'VariableNames',var_names)
results_BM_ud = table(r', mean_BM_ud', sd_BM_ud', cv_BM_ud', 'VariableNames',var_names)
results_fBM_o = table(r', mean_fBM_o', sd_fBM_o', cv_fBM_o', 'VariableNames',var_names)
results_fBM_ud = table(r', mean_fBM_ud', sd_fBM_ud', cv_fBM_ud', 'VariableNames',var_names)

%%

% Plots
figure()
p1 = plot(r, mean_fBM_o, ':', r, mean_fBM_ud, ':', r, mean_BM_o, r, mean_BM_ud, 'LineWidth',1.0)
p1(1).Color = [0.6350 0.0780 0.1840];
p1(2).Color = [0 0.4470 0.7410];
p1(3).Color = [0.6350 0.0780 0.1840];
p1(4).Color = [0 0.4470 0.7410];
xlim([0.5 13])
title("Nuclear export estimates across the tree of life",'FontSize',18)
xlabel("Radius (\mum)",'FontSize',15)
ylabel("Mean first exit time (min)",'FontSize',15)
legend({"Diffusion: obstructed, initial position: origin", ...
    "Diffusion: obstructed, initial position: from uniform dist.", ...
    "Diffusion: normal, initial position: origin", ...
    "Diffusion: normal, initial position: from uniform dist."}, 'FontSize', 12)
ax = gca;
ax.FontSize = 15;

figure()
p2 = plot(r, mean_fBM_o,':', r, mean_fBM_ud,':', r, mean_BM_o, r, mean_BM_ud, 'LineWidth',1.0)
p2(1).Color = [0.6350 0.0780 0.1840];
p2(2).Color = [0 0.4470 0.7410];
p2(3).Color = [0.6350 0.0780 0.1840];
p2(4).Color = [0 0.4470 0.7410];
xlim([3 6])
title("Nuclear export estimates across increasing radii", 'FontSize',18)
xlabel("Radius (\mum)",'FontSize',15)
ylabel("Mean first exit time (min)", 'FontSize',15)
legend({"Diffusion: obstructed, initial position: origin", ...
    "Diffusion: obstructed, initial position: from uniform dist.", ...
    "Diffusion: normal, initial position: origin", ...
    "Diffusion: normal, initial position: from uniform dist."},'FontSize',12)
ax = gca;
ax.FontSize = 15;

figure()
plot(r, mean_fBM_o, r, mean_BM_o, 'LineWidth', 1.0) % plot simulation results
xlim([3 6])
hold on
xline(4, '--', 'Color', "red")
xline(5.5, '--', 'Color', "red")
hold on
scatter(3, 3.36, '>','r','LineWidth',3.5)
title("Nuclear export estimates across increasing radii", 'FontSize',18)
xlabel("Radius (\mum)", 'FontSize',15)
ylabel("Mean export time (min)", 'FontSize',15)
hold on
plot(r, (r.^2)/(6 * 0.4464), 'LineStyle','--') % plot analytical results (
legend({"Simulation results: obstructed diffusion from origin", ...
    "Simulation results: normal diffusion from origin","","","",...
    "Analytical results: normal diffusion from origin"}, 'FontSize', 12)
ax = gca;
ax.FontSize = 15;
hold off

%%

% Distributions

% Distributions for small radius
% r = 3.5, corresponds to the 7th column in data matrices 

% Brownian Motion, initial position at the origin
figure()
histogram(FPT_BM_o(:,7)/60)
title('Distribution of nuclear export times', ...
    'Radius r = 3.5 \mum, BM model, initial position at orgin', ...
    'FontSize',18)
xlabel('Nuclear export time (min)', 'FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% Brownian Motion, initial position drawn from uniform distribution
figure()
histogram(FPT_BM_ud(:,7)/60, 'FaceColor',[0.3010, 0.7450, 0.9330])
title('Distribution of nuclear export times', ...
    'Radius r = 3.5 \mum, BM model, initial position drawn from UD', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% fractional Brownian Motion, initial position at the origin
figure()
histogram(FPT_fBM_o(:,7)/60, 'FaceColor',[0.6350 0.0780 0.1840])
title('Distribution of nuclear export times', ...
    'Radius r = 3.5 \mum, fBM model, initial position at orgin', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% fractional Brownian Motion, initial position drawn from uniform
% distribution
figure()
histogram(FPT_fBM_ud(:,7)/60, 'FaceColor', [0.8500, 0.3250, 0.0980])
title('Distribution of nuclear export times', ...
    'Radius r = 3.5 \mum, fBM model, initial position drawn from UD', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;


% Distribution for "middle" radius
% r = 6, corresponds to the 12th column in data matrices 

% Brownian Motion, initial position at the origin
figure()
histogram(FPT_BM_o(:,12)/60)
title('Distribution of nuclear export times', ...
    'Radius r = 6 \mum, BM model, initial position at orgin', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% Brownian Motion, initial position drawn from uniform distribution
figure()
histogram(FPT_BM_ud(:,12)/60, 'FaceColor',[0.3010, 0.7450, 0.9330])
title('Distribution of nuclear export times', ...
    'Radius r = 6 \mum, BM model, initial position drawn from UD', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% fractional Brownian Motion, initial position at the origin
figure()
histogram(FPT_fBM_o(:,12)/60, 'FaceColor',[0.6350 0.0780 0.1840])
title('Distribution of nuclear export times', ...
    'Radius r = 6 \mum, fBM model, initial position at orgin', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% fractional Brownian Motion, initial position drawn from uniform
% distribution
figure()
histogram(FPT_fBM_ud(:,12)/60, 'FaceColor', [0.8500, 0.3250, 0.0980])
title('Distribution of nuclear export times', ...
    'Radius r = 6 \mum, fBM model, initial position drawn from UD', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;


% Distribution for big radius
% r = 13, corresponds to the 26th column in data matrices

% Brownian Motion, initial position at the origin
figure()
histogram(FPT_BM_o(:,26)/60)
title('Distribution of nuclear export times', ...
    'Radius r = 13 \mum, BM model, initial position at orgin', ...
    'FontSize', 18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% Brownian Motion, initial position drawn from uniform distribution
figure()
histogram(FPT_BM_ud(:,26)/60, 'FaceColor',[0.3010, 0.7450, 0.9330])
title('Distribution of nuclear export times', ...
    'Radius r = 13 \mum, BM model, initial position drawn from UD', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% fractional Brownian Motion, initial position at the origin
figure()
histogram(FPT_fBM_o(:,26)/60, 'FaceColor',[0.6350 0.0780 0.1840])
title('Distribution of nuclear export times', ...
    'Radius r = 13 \mum, fBM model, initial position at orgin', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

% fractional Brownian Motion, initial position drawn from uniform
% distribution
figure()
histogram(FPT_fBM_ud(:,26)/60, 'FaceColor', [0.8500, 0.3250, 0.0980])
title('Distribution of nuclear export times', ...
    'Radius r = 13 \mum, fBM model, initial position drawn from UD', ...
    'FontSize',18)
xlabel('Nuclear export time (min)','FontSize',15)
ylabel('Absolute frequency','FontSize',15)
ax = gca;
ax.FontSize = 15;

