load hessian
%load C2H6.x
load geom_opt.x
C2H6 = geom_opt;

Bohr_Ang = 0.529;

C2H6 = C2H6./Bohr_Ang;

hess = reshape(hessian',[24,24]);
hess1 = hess;

%mass weighted hessian
masses = [12.011, 12.011, 12.011, 12.011, 12.011, 12.011, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007];
masses = masses.*1836.1;

M = diag(masses);
hess = M^(-1/2)*hess*M^(-1/2);

[V,D] = eig(hess);

V = real(V);

do = zeros(24);
for i=1:24
    for j=1:24
        %hess(i,j)=hess(i,j)/sqrt(masses(i)*masses(j));
        do(i,j) = dot(V(:,i),V(:,j));
        %do(i,j) = abs(hess(i,j)-hess(j,i));
    end
end

e = real(eig(hess))

ener = 0;
for i=1:18
    ener = ener + 1/2*sqrt(e(i));
end 
ener

e = e(1:18);
filename = 'eigenvalues';
save(filename,'e','-ascii')

sum_e = sum(e);
e = e./sum_e;
new_masses = e.*sum(masses);

filename = 'new_masses';
save(filename,'new_masses','-ascii')

V = V(:,1:18);

rx = V\C2H6;

r = V*rx;
rxyz = reshape(r,3,8)';

filename = 'inital';
save(filename,'rx','-ascii') 

filename = 'eigenvectors';
save(filename,'V','-ascii')




    




