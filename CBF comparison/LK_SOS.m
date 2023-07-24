% clear
K = [0.3162    2.6124];
vdes = 5;
W = 1.8;
alpha = 1;
% R = sdpvar(1,1);

Y = sdpvar(1,1);
v = sdpvar(1,1);
theta = sdpvar(1,1);
f = [v*(theta-theta^3/3);vdes-v;-K*[Y;theta]];

[h,ch] = polynomial([Y v theta],2,0);
x = [Y,v,theta]';
xbar = [Y,v-vdes,theta]';
dhdt = jacobian(h,x)*f;

[s1,c1] = polynomial([Y v theta],2);
[s2,c2] = polynomial([Y v theta],2);
[s3,c3] = polynomial([Y v theta],4);

Q = diag([1,1,5]);

R=0.1;
dR = 0.1;
feasble_R = R;
while dR>1e-3
    
    Constraints = [sos(h-s1*(R^2-x'*Q*x)),sos(s1),sos(dhdt+alpha*h-s3*((2*R)^2-x'*Q*x)),sos(-h-s2*(Y^2-W^2)),sos(s2),sos(s3),ch(1)==1];
    % Constraints = [sos(h)];
    ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
    
    sol = optimize(Constraints,[],ops,[ch;c1;c2;c3]);
    ch_val=value(ch);
    if sol.problem==0||sol.problem==3||sol.problem==4
        R = R+dR
        feasble_R = R;
    else
        dR = dR*0.7;
    end
end
%%
R = feasble_R;
% R = 1.1;
Constraints = [sos(h-s1*(R^2-xbar'*Q*xbar)),sos(s1),sos(dhdt+alpha*h-s3*((3*R)^2-xbar'*Q*xbar)),sos(-h-s2*(Y^2-W^2)),sos(s2),sos(s3),ch(1)==1];
% Constraints = [sos(h)];
ops = sdpsettings('solver','mosek','verbose',1,'debug',1);

sol = optimize(Constraints,[],ops,[ch;c1;c2;c3]);
ch_val=value(ch);

%%

% ch(1)+Y*ch(2)+v*ch(3)+theta*ch(4)+Y^2*ch(5)+Y*v*ch(6)+v^2*ch(7)+Y*theta*ch(8)+v*theta*ch(9)+theta^2*ch(10)

clear g
g.dim = 3;
g.min = [ -2; 0; -pi/3 ];
g.max = [ 2; 10; pi/3 ];

g.bdry = @addGhostExtrapolate;
g.N = [50,30,30]';
g = processGrid(g);
Ygrid = linspace(g.min(1),g.max(1),g.N(1));
vgrid = linspace(g.min(2),g.max(2),g.N(2));
psigrid = linspace(g.min(3),g.max(3),g.N(3));

data2 = zeros(g.N');

for i=1:g.N(1)
    i
    for j=1:g.N(2)
        for k=1:g.N(3)
            x0 = [Ygrid(i);vgrid(j);psigrid(k)];
            data2(i,j,k) = ch_val(1)+x0(1)*ch_val(2)+x0(2)*ch_val(3)+x0(3)*ch_val(4)+x0(1)^2*ch_val(5)+x0(1)*x0(2)*ch_val(6)+x0(2)^2*ch_val(7)+x0(1)*x0(3)*ch_val(8)+x0(2)*x0(3)*ch_val(9)+x0(3)^2*ch_val(10);
        end
    end
end
[ mesh_xs, mesh_data ] = gridnd2mesh(g, data2);
    
h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
set(h, 'facealpha',0.3,'facecolor','r', 'EdgeColor', 'none');






