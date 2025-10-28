function [m_core,r_core, v_core, n_stars, r_star, v_star] = initializegalaxy(galaxytype)
% Initialize galaxy
% Input: None
% Output: Mass of Core, Position of Core, Velocity of Core, 
% # Stars, Position of Stars, Velocity of Stars

arguments (Output)
    m_core;
    r_core;
    v_core;
    n_stars;
    r_star;
    v_star;
end

G = 1;

circular_orbit = 1;

% Assign order (and value) of mass
m_core = 50; 

% Generate galaxy at rest conditions
r_core = zeros(1,3); 
v_core = zeros(1,3);

% Give core a velocity and position;
if galaxytype == 1 
    v_core = [0.5, 0.1, 0];
    r_core = [-100, 0, 0];
else 
    v_core = [-0.5, -0.1, 0];
    r_core = [100, 0, 0];
end 

% Number of stars
n_stars = randi([5,10]);

% Assign min val and max val according to order of mass
minVal = -m_core;
maxVal = m_core;

% Generate random initial 3D position vector from min val to max val
r_star = r_core + [minVal + (maxVal - minVal)*rand(n_stars,2),zeros(n_stars,1)];
size(r_star);

% Generate random initial 3D velocity vector from min val to max val
if (circular_orbit == 1) 
   
    v_star = zeros(n_stars, 3);

    for i = 1:n_stars

        r_vec = r_star(i,:) - r_core;
        r_xy = sqrt(r_vec(1)^2 + r_vec(2)^2); 

        v_mag = sqrt(G*m_core/(r_xy + 1e-6));

        t_direction = [-r_vec(2),r_vec(1),0]/r_xy;

        v_star(i,:) = v_mag*t_direction + v_core; % initial velocity of star should be relative to the core's velocity 

    end 

else 
    v_star = [minVal + (maxVal - minVal)*rand(n_stars,2),zeros(n_stars,1)];
end 

