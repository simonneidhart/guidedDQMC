load hessian.out
%load C2H6.x
load geom_opt.x
C2H6 = geom_opt;

Bohr_Ang = 0.529;

C2H6 = C2H6./Bohr_Ang;

hess = reshape(hessian',[24,24]);

%mass weighted hessian
masses = [12.011, 12.011, 12.011, 12.011, 12.011, 12.011, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007];
masses = masses.*1836.1;

% for i=1:24
%     for j=1:24
%         hess(i,j)=hess(i,j)/sqrt(masses(i)*masses(j));
%     end
% end

M = diag(masses);

hess = M^(-1/2)*hess*M^(-1/2);

[V,D,W] = eig(hess);

V = real(V);

e = real(eig(hess));

ener = 0;
for i=1:18
    ener = ener + 1/2*sqrt(e(i));
end 
ener

filename = 'eigenvalues';

save(filename,'e','-ascii')

% r = zeros(24,1);
% for i=1:24
%     r = r + e(i)*V(:,i);
% end 
% 
% rxyz = reshape(r,3,8)';
% 
% scatter3(rxyz(:,1),rxyz(:,2),rxyz(:,3))

% A = rref([V(:,1:18) C2H6]);

rx = V\C2H6;

mx = V\masses';
mxyz = V(:,1:18)*mx(1:18);

filename = 'masses';
save(filename,'mx','-ascii')


%rx(3) = 0; %identical eigenvectors
%rx(4) = 0;

r = V(:,1:18)*rx(1:18);

rxyz = reshape(r,3,8)';

c = [2,2,1,1,1,1,1,1];
s = 20.0;
scatter3(rxyz(:,1),rxyz(:,2),rxyz(:,3),s,c,'filled')

filename = 'inital';

rx_out = rx(1:18);

save(filename,'rx_out','-ascii')

filename = 'eigenvectors';

V_out = V(:,1:18);

save(filename,'V_out','-ascii')




    




