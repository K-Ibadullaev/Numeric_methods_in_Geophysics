clear 
clc

%%
% Define the domain
dx = 1; % Spacing for x axis 
dy = 1; % Spacing for x axis 
x = -500:dx:500; % Domain
y = -500:dy:0;   % Domain

% Create the grid
[X, Y] = ndgrid(x,y);
nx = length(x);
ny = length(y);

sig_h = 0.1;                % Conductivity of the earth
sig_d = 10;                 % Conductivity of the dike

k = sig_h*ones(length(x)-1, length(y)-1); % Conductivity matrix

%%
% FOR THE INITIAL ANOMALY:
% k(504: 505, : ) = sig_d;
% or 
k(find(X == 3,1): find(X==4,1),:)= sig_d;

% FOR THE RECTANGULAR ANOMALY:
%k(504: 505, 497:498 ) = sig_d;
% or 
%k(find(X == 3,1): find(X==4,1),find(X == -4,1): find(X==-3,1))= sig_d;

%%

%Defining electrodes / dipoles (Position, Seperation, Currents, Offset)
I_A = -1;                   % Induced current A
I_B = 1;                    % Induced current B
d_ex = -10:1:10;            % Placing of electrodes on x axis
d_ey = zeros(length(d_ex)); % Placing of electrodes in depth
d_p = 1:1:18;               % Distance between emitting and recieving electrodes
d_ab = 1;                   % Distance between electrodes A and B
gN_ = 0;                    % Neumann value.


% Simulating electrical potential distribution 
% Search the boundary indices
[idx_x, idx_y] = ndgrid(1:nx, 1:ny);
idx_x_min = idx_x == 1;
idx_x_max = idx_x == nx;
idx_y_min = idx_y == 1;
idx_y_max = idx_y == ny;
bnd_idx = find(idx_x_min | idx_x_max | idx_y_min | idx_y_max);

%%
% Creating BC
% Homogeneous Dirichlet BC inside the earth
bnd_D = find(idx_y_min | idx_x_max | idx_x_min);
bnd_gD = zeros(length(bnd_D), 1);
% Neumann BC on the surface
bnd_N = find(idx_y_max);
bnd_gN = gN_ + zeros(length(bnd_N), 1);

%%
% Set the source location
for i = 1:1:length(d_p)
    y_0 = d_ey(i);    
    
    for j = 1:1:length(d_ex)
            % Position of electrodes A B M N on x axis
            x_a = d_ex(j);
            x_b = x_a + 1;
            x_m = x_b + d_p(i);
            x_n = x_m + 1;
        
        if(x_b + d_p(i) + 1 <= max(d_ex)) 
            % Indices of electrodes A B M N within X, Y and Phi
            a_idx = find(X == x_a & Y == y_0);
            b_idx = find(X == x_b & Y == y_0);
            m_idx = find(X == x_m & Y == y_0);
            n_idx = find(X == x_n & Y == y_0);
            
            % Potential from the electrode A using Dey & Morrison method
            [A, b_a] = assemble_system('DM2', x, y, k, I_A, a_idx);
            [A, b_a] = apply_NBC(A, b_a, k, x, y, bnd_N, bnd_gN, 'DM2');% First apply Neumann BC
            [A, b_a] = apply_DBC(A, b_a, k, x, y, bnd_D, bnd_gD, 'DM2');% Then apply Dirichlet BC
            Phi_a = A\b_a;
            
            % Potential from the electrode B using Dey & Morrison method
            [A, b_b] = assemble_system('DM2', x, y, k, I_B, b_idx);
            [A, b_b] = apply_NBC(A, b_b, k, x, y, bnd_N, bnd_gN, 'DM2');% First apply Neumann BC
            [A, b_b] = apply_DBC(A, b_b, k, x, y, bnd_D, bnd_gD, 'DM2');% Then apply Dirichlet BC
            Phi_b = A\b_b;
            
            % Calculating potential between M and N 
            U_MN = Phi_a(m_idx) + Phi_b(m_idx) - Phi_a(n_idx) - Phi_b(n_idx);
			
			% Distances between electordes
            r_AM = d_p(i) + 1;
            r_BM = d_p(i);
            r_AN = d_p(i) + 2;
            r_BN = d_p(i) + 1; 
            lns = log(r_AM) - log(r_BM) - log(r_AN) + log(r_BN);
            
			% Apparent resistivity calculation  using the Neumann formula
            rho_a(i,j) = -(U_MN/abs(I_A))* (pi/lns);            
        end
    end
end

%%
% Plot the pseudosection
plot_pseudosection(rho_a, d_ex)






