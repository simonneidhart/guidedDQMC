load eigenvalues
load eigenvalues_mass

m = zeros(18,1);
for i=1:18
    m(i) = eigenvalues(i)/eigenvalues_mass(i);
end

save('n_masses','m','-ascii')
