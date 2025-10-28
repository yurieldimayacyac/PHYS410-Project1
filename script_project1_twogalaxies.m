% This script simulates the interaction of two colliding galaxies, represented 
% using the Toomre Model

% Initialize Galaxy 1
[m_core1, rcore_initial1, vcore_initial1, n_stars1, rstar_initial1, vstar_initial1] = initializegalaxy(1);
N_cores = 1;

% Initialize Galaxy 2 (if testing one galaxy at rest, comment this out)
[m_core2, rcore_initial2, vcore_initial2, n_stars2, rstar_initial2, vstar_initial2] = initializegalaxy(2);
N_cores = N_cores + 1;

if (N_cores ~= 2)
    m_core2 = 0;
    rcore_initial2 = zeros(1,3);
    vcore_initial2 = zeros(1,3);
    n_stars2 = 0;
    rstar_initial2 = zeros(n_stars2,3);
    vstar_initial2 = zeros(n_stars2,3);
end  

% Initialize tmax & time steps of finite difference grid
l = 8;
tmax = 200;
nt = 2^l + 1; % number of time steps
delta_t = tmax/(2^l); % timestep size
t_plot = (0:nt-1)*delta_t; % time array used for plotting

% Store position vector r and mass as an array
% r(Body Index, Coordinate, Number of Time Steps)
r_c = zeros(2,3,nt);    
mass_c = [m_core1, m_core2];

r_star1 = zeros(n_stars1, 3, nt); % Note, stars have vanishing mass
r_star2 = zeros(n_stars2, 3, nt);

% Initial Values of cores 

capture_radius = 5; % initialized 'radius' of core

r_c(1,:,1) = rcore_initial1;
r_c(2,:,1) = rcore_initial2;

v0_c1 = vcore_initial1;
v0_c2 = vcore_initial2;

a0_c1 = nbodyaccn(r_c(1,:,1),r_c(2,:,1),m_core2);
a0_c2 = nbodyaccn(r_c(2,:,1),r_c(1,:,1),m_core1);

% Find core position at n = 2 using Taylor Series
r_c(1,:,2) = r_c(1,:,1) + v0_c1*delta_t + 0.5*a0_c1*delta_t^2;
r_c(2,:,2) = r_c(2,:,1) + v0_c2*delta_t+ 0.5*a0_c2*delta_t^2;

% Initial Values of stars 
r_star1(:,:,1) = rstar_initial1;
r_star2(:,:,1) = rstar_initial2; 

v0_star1 = zeros(n_stars1,3);
v0_star2 = zeros(n_stars2,3);

v0_star1(:,:) = vstar_initial1;
v0_star2(:,:) = vstar_initial2;

a0_star1 = zeros(n_stars1, 3);
a0_star2 = zeros(n_stars2, 3);

% Initialize stars from galaxy 1

if n_stars1 > 0

    for i = 1:n_stars1
       
        acc = [0 0 0];
        for j = 1:N_cores
            if mass_c(j) ~= 0
                acc = acc + nbodyaccn(r_star1(i,:,1), r_c(j,:,1), mass_c(j));
            end
        end
        a0_star1(i,:) = acc;

        r_star1(i,:,2) = r_star1(i,:,1) + delta_t * v0_star1(i,:) + delta_t^2 * 0.5 * a0_star1(i,:);

    end 

end 

% Initialize stars from galaxy 2 

if n_stars2 > 0

    for i = 1:n_stars2
       
        acc = [0 0 0];
        for j = 1:N_cores
            if mass_c(j) ~= 0
                acc = acc + nbodyaccn(r_star2(i,:,1), r_c(j,:,1), mass_c(j));
            end
        end
        a0_star2(i,:) = acc;
        r_star2(i,:,2) = r_star2(i,:,1) + delta_t * v0_star2(i,:) + delta_t^2 * 0.5 * a0_star2(i,:);
    end 
end 


% Finite Difference Solver

for n = 2 : nt-1

    % Core - Core reaction
    for i = 1 : N_cores

        total_acceleration = zeros(1,3);
    
        for j = 1 : N_cores
    
            % Do not count core's interaction with itself
            if i == j 
                continue 
            end 

            % Calculate acceleration based on gravity
            accel = nbodyaccn(squeeze(r_c(i,:,n)), squeeze(r_c(j,:,n)), mass_c(j));
            total_acceleration = total_acceleration + accel;

        end 

        % Add advanced position into position vector
        r_c(i,:,n+1) = 2*r_c(i,:,n) - r_c(i,:,n-1) + total_acceleration * delta_t^2;

    end 
    
    % Each core's stars interactions with the other core
    for m = 1 : N_cores

        % Set number of stars according to which core's stars we are
        % focusing on
        if (m == 1)
            N_star = n_stars1;
            r_s = r_star1;
        elseif (m == 2)
            N_star = n_stars2;
            r_s = r_star2;
        end 
    
        % Star - Core reaction
        for i = 1 : N_star
    
            total_acceleration = zeros(1,3);

            capture_core = 0;
        
            for j = 1: N_cores
        
                % Calculate acceleration based on gravity
                accel = nbodyaccn(r_s(i,:,n), r_c(j,:,n), mass_c(j));
    
                % Superpose acceleration of i-th star due to j-th core
                total_acceleration = total_acceleration + accel;

                if dot(squeeze(r_s(i,:,n))-squeeze(r_c(j,:,n)), squeeze(r_s(i,:,n))-squeeze(r_c(j,:,n))) < capture_radius^2
                    capture_core = j; 
                end 

            end 
    
            % Add advanced position into position vector
            if capture_core ~= 0
                r_s(i,:,n) = r_c(capture_core,:,n);
                r_s(i,:,n+1) = r_c(capture_core,:,n+1);
            else 
                r_s(i,:,n+1) = 2*r_s(i,:,n) - r_s(i,:,n-1) + total_acceleration * delta_t^2;
            end 
        end 

        % Save position vectors
        if (m == 1)
            r_star1 = r_s;
        elseif (m == 2)
            r_star2 = r_s;
        end 

    end 
    
end 

%% Plot the trajectories of the cores

x_core1 = squeeze(r_c(1, 1, :));
y_core1 = squeeze(r_c(1,2,:));
x_core2 = squeeze(r_c(2, 1, :));
y_core2 = squeeze(r_c(2,2,:));

x_star1 = squeeze(r_star1(:,1,:));
y_star1 = squeeze(r_star1(:,2,:));
x_star2 = squeeze(r_star2(:,1,:));
y_star2 = squeeze(r_star2(:,2,:));

figure; 
hold on;
grid on;

% Plot the cores
for m = 1:N_cores
    plot (r_c(m,1,n),r_c(m,2,n),'o','MarkerSize',18,'MarkerFaceColor', ...
        [0.1 0.5 1], 'MarkerEdgeColor','k','DisplayName','Core 1');
end 

for i = 1:n_stars1

    color = rand(1,3);

    x_star = squeeze(r_star1(i,1,:));
    y_star = squeeze(r_star1(i,2,:));
    plot(x_star, y_star, 'o', 'Color', color, 'DisplayName', sprintf('Star %d', i));

end 

% Plot the stars from galaxy 2
for i = 1:n_stars2
    color = rand(1,3);
    x_star = squeeze(r_star2(i,1,:));
    y_star = squeeze(r_star2(i,2,:));
    plot(x_star, y_star, 'o', 'Color', color, 'DisplayName', sprintf('Star %d', i + n_stars1));
end 

%% Animation

pausesecs = 0;
avienable = 0; % toggle recording
avifilename = 'galaxymotion.avi';
aviframerate = 25;

if avienable
   aviobj = VideoWriter(avifilename);
   open(aviobj);
end

% Initialize range for plotting
x_array = [x_core1(:); x_core2(:); x_star1(:); x_star2(:)];
y_array = [y_core1(:); y_core2(:); y_star1(:); y_star2(:)];

xmin = min(x_array);
xmax = max(x_array);
ymin = min(y_array);
ymax = max(y_array);

% Change
cx0 = mean([r_c(1,1,1), r_c(2,1,1)]);  % x-center of both cores
cy0 = mean([r_c(1,2,1), r_c(2,2,1)]);  % y-center of both cores

% Fix the half-width of the camera window
half = 150;  % adjust to zoom level you like

if xmin == xmax
    xmin = xmin - 1;
    xmax = xmax + 1;
end
if ymin == ymax
    ymin = ymin - 1;
    ymax = ymax + 1;
end

maxrange = max(xmax-xmin,ymax-ymin);

core1_color      = [0.90 0.20 0.20];
core2_color      = [0.20 0.35 0.95];
g1_star_color    = [0.40 0.70 1.00];
g2_star_color    = [0.90 0.85 0.30];

for k = 1:nt

    t = (k-1)*delta_t;

    % Graphics
    clf;
    hold on;

    % Fixed camera centered at the initial cores (compute cx0, cy0, half once)
    ax = gca;
    xlim(ax, [cx0 - half, cx0 + half]);
    ylim(ax, [cy0 - half, cy0 + half]);
    
    axis(ax, 'equal');   % keep aspect ratio 1:1
    axis(ax, 'manual');  % <-- FREEZE limits; plotting won't change them

    axis equal;


    titlestr = sprintf('Step: %d   Time: %.1f', k, t),;
    title(titlestr, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'g');

    % Plot trails of galaxy 1
    plot(x_core1(1:k),y_core1(1:k),'-','LineWidth',1,'Color', core1_color);

    for s = 1 : n_stars1 
        plot(squeeze(r_star1(s,1,1:k)),squeeze(r_star1(s,2,1:k)), '-','LineWidth',1,'Color', g1_star_color);

    end 
    
    % Plot bodies of galaxy 1

    plot(x_core1(k),y_core1(k), 'o', 'MarkerSize', 16, ...
             'MarkerFaceColor', core1_color, 'MarkerEdgeColor', 'k');

    for s = 1 : n_stars1 
        plot(squeeze(r_star1(s,1,k)),squeeze(r_star1(s,2,k)), 'o', 'MarkerSize', 5, ...
             'MarkerFaceColor', g1_star_color, 'MarkerEdgeColor', 'k');
    end 

    if m_core2 ~= 0 % if galaxy 2 exists in code

        % Plot trails of galaxy 2

        plot(x_core2(1:k), y_core2(1:k), '-', 'LineWidth', 1, 'Color', core2_color);

        for s = 1:n_stars2
            plot(squeeze(r_star2(s,1,1:k)),squeeze(r_star2(s,2,1:k)), '-', 'LineWidth', 1, 'Color', g2_star_color);
        end

        % Plot bodies of galaxy 2

        plot(x_core2(k), y_core2(k), 'o', 'MarkerSize', 16, ...
         'MarkerFaceColor', core2_color, 'MarkerEdgeColor', 'k');

        for s = 1:n_stars2
            plot(squeeze(r_star2(s,1,k)),squeeze(r_star2(s,2,k)), 'o', 'MarkerSize', 4, ...
                 'MarkerFaceColor', g2_star_color, 'MarkerEdgeColor', g2_star_color);
        end

    end 

    drawnow; 

    if avienable
         if t == 0
            framecount = 1 * aviframerate;
         else
            framecount = 1;
         end
         for iframe = 1 : framecount
            writeVideo(aviobj, getframe(gcf));
         end

    end

    pause(pausesecs);
    
end 

if avienable
    close(aviobj);
end 