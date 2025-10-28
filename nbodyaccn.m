function accel = nbodyaccn(r_ci, r_cj, mass_j)

G = 1;

d = r_cj - r_ci;
softening_parameter = 1e-6;
distance = sqrt((norm(d))^2 + softening_parameter^2);

accel = G * (mass_j/distance^3)*(d);
end 