clear g
g.dim = 3;
g.min = [ -2; 0; -1.5 ];
g.max = [ 2; 10; 1.5 ];

g.bdry = @addGhostExtrapolate;
g.N = [50,30,30]';
g = processGrid(g);
Ygrid = linspace(g.min(1),g.max(1),g.N(1));
vgrid = linspace(g.min(2),g.max(2),g.N(2));
psigrid = linspace(g.min(3),g.max(3),g.N(3));
data1 = zeros(g.N');
amax = 5;
rmax = 0.5;
W = 1.8;

x0 = [1;10;0.1];
x = x0;
A = [0 5;0 0];
B = [0;1];
[K,S,e] = lqr(A,B,diag([1;5]),1,zeros(2,1));

con1 = @(x)[clip(5-x(2),amax);clip(-K*[x(1);x(3)],rmax)];
% con1 = @(x)[max(-amax,-10*x(2));clip(-9*K*[x(1);x(3)],rmax)];
f = @(t,x)dubin(t,x,con1);
Opt    = odeset('Events', @stop_crit);

for i=1:g.N(1)
    i
    for j=1:g.N(2)
        for k=1:g.N(3)
            x0 = [Ygrid(i);vgrid(j);psigrid(k)];
            [T, Y] = ode45(f, [0,5], x0, Opt);
            data1(i,j,k) = min(W - max(abs(Y(:,1))),pi/3-max(abs(Y(:,3))));
        end
    end
end

%%
tic
[ data, g, data0 ] = Lane_keeping_HJ('high');
toc

% f = @(t,x)dubin(t,x,con1);
% x0 = [0;5.1;-0.7];
% [T, Y] = ode45(f, [0,5], x0, Opt);
%% SOS
tic
vdes = 5;
W = 1.8;
alpha = 1;
% R = sdpvar(1,1);

Y = sdpvar(1,1);
v = sdpvar(1,1);
theta = sdpvar(1,1);
f = [v*(theta-theta^3/3);vdes-v;-K*[Y;theta]];

[h,ch] = polynomial([Y v theta],4,0);
x = [Y,v,theta]';
xbar = [Y,v-vdes,theta]';
dhdt = jacobian(h,x)*f;

[s1,c1] = polynomial([Y v theta],4);
[s2,c2] = polynomial([Y v theta],4);
[s3,c3] = polynomial([Y v theta],4);
[s4,c4] = polynomial([Y v theta],6);

Q = diag([1,1,5]);

R = 1.1;
Constraints = [sos(h-s1*(R^2-xbar'*Q*xbar)),sos(s1),sos(dhdt+alpha*h-s4*((3*R)^2-xbar'*Q*xbar)),sos(-h-s2*(Y^2-W^2)),sos(-h-s3*(theta^2-(pi/6)^2)),sos(s2),sos(s3),sos(s4),ch(1)==1];
% Constraints = [sos(h)];
ops = sdpsettings('solver','mosek','verbose',1,'debug',1);

sol = optimize(Constraints,[],ops,[ch;c1;c2;c3;c4]);
ch=value(ch);

data2 = zeros(g.N');

for i=1:g.N(1)
    for j=1:g.N(2)
        for k=1:g.N(3)
            x0 = [Ygrid(i);vgrid(j);psigrid(k)];
            Y = Ygrid(i);
            v = vgrid(j);
            theta = psigrid(k);
            data2(i,j,k) = ch(1)+Y*ch(2)+v*ch(3)+theta*ch(4)+Y^2*ch(5)+Y*v*ch(6)+v^2*ch(7)+Y*theta*ch(8)+v*theta*ch(9)+theta^2*ch(10)...
                +Y^3*ch(11)+Y^2*v*ch(12)+Y*v^2*ch(13)+v^3*ch(14)+Y^2*theta*ch(15)+Y*v*theta*ch(16)+v^2*theta*ch(17)+Y*theta^2*ch(18)+...
                v*theta^2*ch(19)+theta^3*ch(20)+Y^4*ch(21)+Y^3*v*ch(22)+Y^2*v^2*ch(23)+Y*v^3*ch(24)+v^4*ch(25)+Y^3*theta*ch(26)+Y^2*v*theta*ch(27)...
                +Y*v^2*theta*ch(28)+v^3*theta*ch(29)+Y^2*theta^2*ch(30)+Y*v*theta^2*ch(31)+v^2*theta^2*ch(32)+Y*theta^3*ch(33)+v*theta^3*ch(34)+theta^4*ch(35);
        end
    end
end
toc


%% plotting
figure(1)
clf
subplot(211)
[ mesh_xs, mesh_data ] = gridnd2mesh(g, data);
    
h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
% isonormals(mesh_xs{:}, mesh_data, h);
set(h, 'facealpha',0.3,'facecolor','g', 'EdgeColor', 'none');
[ mesh_xs, mesh_data ] = gridnd2mesh(g, data1);
    
h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
set(h, 'facealpha',0.3,'facecolor','r', 'EdgeColor', 'none');

[ mesh_xs, mesh_data ] = gridnd2mesh(g, data2);
    
h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
set(h, 'facealpha',0.3,'facecolor','b', 'EdgeColor', 'none');
isonormals(mesh_xs{:}, mesh_data, h);
axis([-2,2,0.5,10,-1.8,1.8])
xlabel('$Y$','interpreter','latex','fontsize',20)
ylabel('$v$','interpreter','latex','fontsize',20)
zlabel('$\psi$','interpreter','latex','fontsize',20)
m = legend('HJ PDE','Backup CBF','SOS');
set(m,'fontsize',15,'edgecolor','none')
subplot(212)
hold on
[X,Y]=meshgrid(Ygrid,psigrid);

contour(X',Y',squeeze(data(:,16,:)),[0,0],'color','g','linewidth',2)

contour(X',Y',squeeze(data1(:,16,:)),[0,0],'color','r','linewidth',2)

contour(X',Y',squeeze(data2(:,16,:)),[0,0],'color','b','linewidth',2)
xlabel('$Y$','interpreter','latex','fontsize',20)

ylabel('$\psi$','interpreter','latex','fontsize',20)
title('$v=5.1m/s$','interpreter','latex','fontsize',20)
grid on
m = legend('HJ PDE','Backup CBF','SOS');
set(m,'fontsize',15,'edgecolor','none')
tightfig

%%
figure
hold on
% [X,Y]=meshgrid(Ygrid,psigrid);
% contour(X',Y',squeeze(data(:,16,:)),[0,0],'color','g','linewidth',2)
% contour(X',Y',squeeze(data1(:,16,:)),[0,0],'color','r','linewidth',2)
% 
% xlabel('$Y$','interpreter','latex','fontsize',20)
% ylabel('$\psi$','interpreter','latex','fontsize',20)
% grid on
% m = legend('HJ PDE','Backup CBF');
% set(m,'fontsize',15,'edgecolor','none')
[ mesh_xs, mesh_data ] = gridnd2mesh(g, data);
    
h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
% isonormals(mesh_xs{:}, mesh_data, h);
set(h, 'facealpha',0.3,'facecolor','g', 'EdgeColor', 'none');
[ mesh_xs, mesh_data ] = gridnd2mesh(g, data1);
    
h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
set(h, 'facealpha',0.3,'facecolor','r', 'EdgeColor', 'none');
isonormals(mesh_xs{:}, mesh_data, h);
axis([-2,2,0.5,10,-1.8,1.8])
xlabel('$Y$','interpreter','latex','fontsize',20)
ylabel('$v$','interpreter','latex','fontsize',20)
zlabel('$\psi$','interpreter','latex','fontsize',20)
m = legend('HJ PDE','Backup CBF','SOS');
set(m,'fontsize',15,'edgecolor','none')
tightfig

%%
deriv={};
for i=1:3
    deriv{i}=0*data;
end
deriv{1}(1:g.N(1)-1,:,:)=(data(2:g.N(1),:,:)-data(1:g.N(1)-1,:,:))/(g.xs{1}(2,1,1)-g.xs{1}(1,1,1));
deriv{2}(:,1:g.N(2)-1,:)=(data(:,2:g.N(2),:)-data(:,1:g.N(2)-1,:))/(g.xs{2}(1,2,1)-g.xs{2}(1,1,1));
deriv{3}(:,:,1:g.N(3)-1)=(data(:,:,2:g.N(3))-data(:,:,1:g.N(3)-1))/(g.xs{3}(1,1,3)-g.xs{3}(1,1,1));
options = optimoptions('quadprog','Display','none'); 
H = diag([1,1,0]);


%%
tic
for m = 1:100
    i=randsample(g.N(1),1);
    j=randsample(g.N(2),1);
    k=randsample(g.N(3),1);
    x = [g.xs{1}(i,j,k);g.xs{2}(i,j,k);g.xs{3}(i,j,k)];
    h = data(i,j,k);
    dh = [deriv{1}(i,j,k);deriv{2}(i,j,k);deriv{3}(i,j,k)];
    u0 = [0.1;0.1];
    Lfh = dh'*[x(2)*sin(x(3));0;0];
    Lgh = dh(2:3)';
    Aineq = [-Lgh -1];
    bineq = h+Lfh;

    [res,~,~] = quadprog(H,[-u0;1e6],Aineq,bineq,[],[],[-amax;-rmax;0],[amax;rmax;inf],[u0;0],options);
%         settings.verbose = 0;
%         params.H = diag([1,1,0]);
%         params.A = Aineq;
%         params.b = bineq;
%         params.lb = [-amax;-rmax;0];
%         params.ub = [amax;rmax;1e6];
%         params.f = 2*[-u0; 1e6];
%         [vars, status] = csolve(params, settings);
%         u = vars.u(1:2);
end
toc

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

