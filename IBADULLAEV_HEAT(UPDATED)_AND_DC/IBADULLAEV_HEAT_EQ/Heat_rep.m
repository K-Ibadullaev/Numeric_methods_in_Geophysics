% The Heat equation

clear;
clc;
%% Daily temperature oscillations

% initial set-up
omega = 2*pi; % frequency
kappa = 2e-7 * 3600 * 24; % heat conduction constant
dt = 0.001; % time step
t = 0:dt:1; % time span
dz = 0.01; % spatial spacing
z = (0:dz:10).'; % depth
T0 = 12.*ones(length(z),1); % initial value
dT = 20; % temperature oscillation, K
D = assemble_D(z, dz); % operator for central differences


% Calculate analytical solution
T_an = T0 + dT * exp(-z*sqrt(omega/(2*kappa))) .* ...
            sin(omega.*t - z.*sqrt(omega/(2*kappa)));

% Allocate size of vectors and set starting values
deviation_implicit = zeros(size(t));
deviation_explicit = zeros(size(t));
deviation_cn = zeros(size(t));

T_expl = T_an(2:end-1,:);
T_impl = T_an(2:end-1,:);
T_crni = T_an(2:end-1,:);

% Set Boundary Conditions
Ta = T0(1) + dT * sin(omega*dt);% at the surface the temperature changes periodicaly
Tb = T0(1); % at Zmax the temperature is constant
[D, T_bc] = applyDirichletBC(D, Ta, Tb, dz); % implement Dirichlet BC


%% Computation and visualization of daily oscillations

figure(1);

if dt<dz^2/(2 *kappa)
   sprintf("Stable!")
   
   for i = 2:numel(t)

%     % Timedependent Boundary Condition
    Ta = T0(1) + dT * sin(omega*i*dt);
    T_bc(1) = Ta/(dz^2);
%     
%     % Explicit Euler Method 
    T_expl(:,i) = (eye(length(D)) + kappa * dt * D) ...
        * T_expl(:,i-1) + kappa * dt * T_bc;  
    
%     % Implicit Euler Method
    T_impl(:,i) = (eye(length(D)) - kappa * dt * D) ...
                \ (T_impl(:,i-1) + kappa * dt * T_bc);  
            
%     % Crank Nicolson Method
    T_crni(:,i) = (eye(length(D)) - kappa * dt/2 * D) ...
              \ ((eye(length(D)) + kappa * dt/2 * D) ...
              * T_crni(:,i-1) + kappa * dt * T_bc);

%     % Calculate deviations from analytical solution
    deviation_explicit(i) = norm(T_an(2:end-1,i) - T_expl(:,i), 'inf');
    deviation_implicit(i) = norm(T_an(2:end-1,i) - T_impl(:,i), 'inf');
    deviation_cn(i)       = norm(T_an(2:end-1,i) - T_crni(:,i), 'inf');

    % Plot iterative solutions ontop of analytical solution
    if i < numel(t)
        plot(z,T_an(:,i), 'r-')
        hold on
        plot(z(2:end-1),T_expl(:,i), 'green*')    
        plot(z(2:end-1),T_impl(:,i), 'yellow+')
        plot(z(2:end-1),T_crni(:,i), 'mo')
        ylim([-dT, dT] + T0(1));
        xlim([0,0.6]); % doesn't changes after 0.5
        title(sprintf('Hour %.1f', i/(1000/24)));
        legend("Analytical Solution","Explicit Euler","Implicit Euler","Crank-Nicolson")
        xlabel 'Depth in [m]'
        ylabel 'Temperature in [°C]'
         
              
    end
    hold off
    pause(0.01)
        
    
    
  end

else
    sprintf("Unstable!")
end

%%
%listing depth

% analytical
[row_an,col_an] = find((T_an(1,:)==-8));
[z_a, c_a] = find(T_an(:,col_an) >=0,1);
z_an = z_a * dz;
sprintf( 'The first positive temperature value of %.3f°C occurs at %.2f meters',T_an(z_a,col_an),z_an)

% explicit
% [row_ex,col_ex] = find((T_expl(1,:)==-8));
% [z_ex, c_ex] = find(T_expl(:,col_ex) >=0,1);
% z_ex = z_ex * dz

% implicit
% [row_impl,col_impl] = find((T_impl(1,:)==-8));
% [z_impl, c_impl] = find(T_impl(:,col_impl) >=0,1);
% z_impl = z_impl * dz
% 
% crank nicolson
% [row_crni,col_crni] = find((T_crni(1,:)==-8));
% [z_crni, c_crni] = find(T_an(:,col_crni) >=0,1);
% z_crni = z_crni * dz

%%

figure(2)
% Plot Deviations 
plot(1:length(deviation_explicit), deviation_explicit, 'g', ...
     1:length(deviation_explicit), deviation_implicit, 'y', ...
     1:length(deviation_explicit), deviation_cn,       'm')
 legend("Explicit Euler","Implicit Euler","Crank-Nicolson Method")
 ylabel 'Deviation in [K]'
 xlabel 'Iteration Number'
 title 'Deviation from Analytical Solution'
 mean(deviation_explicit(:,:))
 mean(deviation_implicit(:,:))
 mean(deviation_cn(:,:))
 
 %% Yearly oscillations
 
clear;
clc;
%% set-up
 
omega = (2*pi);
kappa = 2e-7 * 3600 * 24 * 365;
dt = 0.001;
t = 0:dt:1;
dz = 0.25;
z = (0:dz:10).';
T0 = 12.*ones(length(z),1);
dT = 20;
D = assemble_D(z, dz);
 
% Calculate Analytical Solution
T_an = T0 + dT * exp(-z * sqrt(omega / (2 * kappa))) .* ...
            sin(omega .* t - z .* sqrt(omega / (2 * kappa)));

% allocate size of vectors and set starting values
T_expl = T_an(2:end-1,:);
T_impl = T_an(2:end-1,:);
T_crni = T_an(2:end-1,:);

deviation_implicit = zeros(size(t));
deviation_explicit = zeros(size(t));
deviation_cn = zeros(size(t));

%% Computation and visualization of yearly oscillations
if dt<dz^2/(2 *kappa)
   sprintf("Stable!")
   
for i = 1:numel(t)
   
    if i == 1        
         % Boundary Conditions
        Ta = T0(1) + dT * sin(omega*i*dt);
        Tb = T0(1);
        [D, T_bc] = applyDirichletBC(D, Ta, Tb, dz);
        
    else
        % Timedependent Boundary Condition
        Ta = T0(1) + dT * sin(omega*i*dt);
        T_bc(1) = Ta/(dz^2);
       
        % Explicit Euler Method
        T_expl(:,i) = (eye(length(D)) + kappa * dt * D) ...
            * T_expl(:,i-1) + kappa * dt * T_bc;  
        
        % Implicit Euler Method
        T_impl(:,i) = (eye(length(D)) - kappa * dt * D) ...
                    \ (T_impl(:,i-1) + kappa * dt * T_bc);  
                
        % Crank Nicolson Method
        T_crni(:,i) = (eye(length(D)) - kappa * dt/2 * D) ...
                  \ ((eye(length(D)) + kappa * dt/2 * D) ...
                  * T_crni(:,i-1) + kappa * dt * T_bc);
    end
    
     % calculate deviations from analytical solution
    deviation_explicit(i) = norm(T_an(2:end-1,i) - T_expl(:,i), 'inf');
    deviation_implicit(i) = norm(T_an(2:end-1,i) - T_impl(:,i), 'inf');
    deviation_cn(i)       = norm(T_an(2:end-1,i) - T_crni(:,i), 'inf');

    % plot iterative methods ontop of analytical solution
    if i < numel(t)
    plot(z,T_an(:,i), 'r-')
    hold on
    plot(z(2:end-1),T_expl(:,i), 'green*')    
    plot(z(2:end-1),T_impl(:,i), 'yellow+')
    plot(z(2:end-1),T_crni(:,i), 'mo')
    ylim([-dT, dT] + T0(1));
    xlim([0,10]);
    title(sprintf('Day %.1f', i/(1000/365)));
    legend("Analytical Solution","Explicit Euler","Implicit Euler",...
    "Crank-Nicolson Method")
    xlabel 'Depth in [m]'
    ylabel 'Temperature in [°C]'
%     anim = getframe(h);
%         im = frame2im(anim);
%         [imind,cm] = rgb2ind(im,256);
%         count = 1;
%         Write to the GIF File
%         if count == 1
%             imwrite(imind,cm,"yearly_o.gif", 'GIF', 'Loopcount',inf);
%         else
%             imwrite(imind,cm,"yearly_o.gif", 'GIF', 'WriteMode','append');
%         end
    end
    hold off
    pause(0.01)
end
else
    sprintf("Unstable!")
end

% Plot Deviations 
plot(1:length(deviation_explicit), deviation_explicit, 'g', ...
     1:length(deviation_explicit), deviation_implicit, 'y', ...
     1:length(deviation_explicit), deviation_cn,       'm')
 legend("Explicit Euler","Implicit Euler","Crank-Nicolson Method")
 ylabel 'Deviation in [K]'
 xlabel 'Iteration Number'
 title 'Deviation from Analytical Solution'
 mean(deviation_explicit(:,:))
 mean(deviation_implicit(:,:))
 mean(deviation_cn(:,:))
 
 %% listing depth where first positive temperature occurs
     %when temperature on the surface is -8°C

[row_an,col_an] = find((T_an(1,:)==-8),5);
[z_a, c_a] = find(T_an(:,col_an) >=0,1);
z_an = z_a * dz;
sprintf( 'The first positive temperature value of %.3f°C occurs at %.2f meters',T_an(z_a,col_an),z_an)
 %% Local functions
%function for assembling the operator D
 function D = assemble_D(z, dz)
    % Build sparse 2nd-order central difference operator for uniform mesh.
    %
    % z  ... boundary value at z_min 
    % dz ... boundary value at z_max
    
    % Initialize.
    n = length(z);    
    D_kern = repmat([1, -2, 1], n, 1);
    
    % Assemble.
    D = 1/dz^2 * spdiags(D_kern, [-1, 0, 1], n, n);    
 end

%Dirichlet boundary conditions
function [D, b] = applyDirichletBC(D, a, e, dz)
    % Map Dirichlet values from D to b and reduce system size.
    %
    % a ... boundary value at z_min 
    % e ... boundary value at z_max

    b = zeros(size(D, 1), 1);
    if ~isempty(a)
        % Remove entries related to boundary z_min;
        D(:, 1) = [];
        D(1, :) = [];
        b(1) = [];
        
        % Map Dirichlet values.
        b(1) = a/dz^2;
    end
    if ~isempty(e)
        % Remove entries related to boundary z_max;
        D(:, end) = [];
        D(end, :) = [];
        b(end) = [];
        
        % Map Dirichlet values.
        b(end) = e/dz^2;
    end
end