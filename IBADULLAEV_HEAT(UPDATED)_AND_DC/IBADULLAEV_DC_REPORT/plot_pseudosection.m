function plot_pseudosection(rho, d_ex)
% Plots rho defined at cells dependent on positioning d_ex of electrodes
%
% INPUT PARAMETER
%   rho     ...     resistivity ,upper left matrix 
%   d_ex    ...     vector for positions of electrode A on X axis

% We move electrodes MN along the survey line stepwise
% Define coordinates of "blocks" for the pseudosection
X = zeros(4, 0);
Y = zeros(4, 0);
    for i = 1:numel(d_ex)
        yref = i;
        for j = 1:numel(d_ex)
            xref = d_ex(j) + (i+2)/2;
            if(xref <= d_ex(end) - (i+2)/2)
            X(:,end+1) = [xref - 0.5; xref - 0.5; xref + 0.5; xref + 0.5];
            Y(:,end+1) = [yref + 0.5; yref - 0.5; yref - 0.5; yref + 0.5];
            end
        end
    end
    
    patch(X,Y,nonzeros(rho'))
	
    xlabel('x_{ref} in [m]');
    ylabel('d_{x}^p in [m]');
    axis('equal');
    set(gca, 'YDir','reverse') %colormap
    hcb = colorbar();
    hcb.Title.String = "rho_{a}";
    ylim([0, 20]);
    xlim([d_ex(1),d_ex(end)]);
    title('Pseudosection')
end
