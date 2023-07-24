clear g
g.dim = 3;
g.min = [ -5; -5; -pi ];
g.max = [ 5; 5; pi ];

g.bdry = @addGhostExtrapolate;
g.N = [41,41,41]';
g = processGrid(g);
Xgrid = linspace(g.min(1),g.max(1),g.N(1));
Ygrid = linspace(g.min(2),g.max(2),g.N(2));
psigrid = linspace(g.min(3),g.max(3),g.N(3));
data1 = zeros(g.N');
um = 0.2;
R = 1;
va = 1;
vb = 1;
Ts = 0.04;

con1 = @(x)-sign(x(2))*um;

f = @(t,x)Air3D(t,x,con1(x),va,vb);
Opt    = odeset('Events', @stop_crit);

for i=1:g.N(1)
    i
    for j=1:g.N(2)
        for k=1:g.N(3)
            x0 = [Xgrid(i);Ygrid(j);psigrid(k)];
            x = x0;
            Y = [];
%             [T, Y] = ode45(f, [0,2], x0, Opt);
            for t=0:Ts:2
                x = x+f(t,x)*Ts;
                Y = [Y;x'];
            end
            data1(i,j,k) = min(vecnorm(Y(:,1:2)')-R);
        end
    end
end

%%
data = Air_HJ();
data_compressed = zeros(size(data,1),size(data,2),size(data,3));
for i=1:size(data,1)
    for j=1:size(data,2)
        for k=1:size(data,3)
            data_compressed(i,j,k)=min(data(i,j,k,:));
        end
    end
end

%% plotting
figure(1)
clf
subplot(121)
title('Hamilton Jacobi Isaac result','fontsize',15)
h = visSetIm(g, data_compressed, 'b', 0);
xlabel('$\Delta X$','interpreter','latex','fontsize',15)
ylabel('$\Delta Y$','interpreter','latex','fontsize',15)
zlabel('$\Delta \psi$','interpreter','latex','fontsize',15)
% [ mesh_xs, mesh_data ] = gridnd2mesh(g, data);
    
% h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
% isonormals(mesh_xs{:}, mesh_data, h);
% set(h, 'facealpha',0.3,'facecolor','g', 'EdgeColor', 'none');
subplot(122)
title('Backup CBF result','fontsize',15)
hold on
h = visSetIm(g, data1, 'r', 0);
xlabel('$\Delta X$','interpreter','latex','fontsize',15)
ylabel('$\Delta Y$','interpreter','latex','fontsize',15)
zlabel('$\Delta \psi$','interpreter','latex','fontsize',15)

% xlabel('$Y$','interpreter','latex','fontsize',20)
% ylabel('$v$','interpreter','latex','fontsize',20)
% zlabel('$\psi$','interpreter','latex','fontsize',20)
% m = legend('HJ PDE','Backup CBF','SOS');
% set(m,'fontsize',15,'edgecolor','none')

% tightfig


%%


% function con_des(x,x_des)
% des_heading = atan2(x_des(2)-x(2),x_des(1)-x(1));
% heading_error = x(4)-des_heading;
% while heading_error>pi
%     heading_error = heading_error-2*pi;
% end
% while heading_error<-pi
%     heading_error = heading_error+2*pi;
% end
%
% end
function res = clip(x,xm)
res = x/(norm(x)+1e-6)*min(norm(x),xm);
end
function dxdt = dubin(t,x,con)
u = con(x);
dxdt = [x(2)*sin(x(3));u(1);u(2)];
end
function [value, isterminal, direction] = stop_crit(t, x)
value      = double(abs(x(1))>2||t>5||(abs(x(1))<0.1&&abs(x(3))<0.05));
isterminal = 1;   % Stop the integration
direction  = 0;
end
% function [b,db] = obs(x,center,r)
% b = (x(1:2)-center)'*(x(1:2)-center)-r^2;
% db = [2*x(1:2)-center;zeros(2,1)];
% end
function dx = Air3D(t,x,u,va,vb)
dx =zeros(3,1);
dx(1) = -va + vb * cos(x(3)) + u*x(2);
dx(2) = vb * sin(x(3)) - u*x(1);
dx(3) = - u;
end
