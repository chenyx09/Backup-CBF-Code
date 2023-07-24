function [u,backup_state] = backup_CBF_con_v1(t,x,u0,bdry,backup_state)
global Ts vm rm am
bs = backup_state;
Ts = 0.05;
vm = 0.18;
rm = 1.2;
am = 0.3;
con1 = @(x)[-am*sign(x(3));0];
con2 = @(x)[-am*sign(x(3)/2);rm];
con3 = @(x)[-am*sign(x(3)/2);-rm];
stop_crit = @(t,x)(abs(x(3))<=0);

% if isempty(tt0)
%     [tt0,xx0,QQ0] = generate_backup_traj([0;0;vm;0],con1,stop_crit);
% end
alpha = 1;

r = 0.17;


%% x = [X,Y,v,theta], u = [a,r]
N_robot = size(x,2);
u = zeros(2,N_robot);


H = diag([1 1 0]);

% options = qpOASES_options('printLevel',0);
options = optimoptions('quadprog','Display','none');
% x_dot_backup = zeros(4,N_robot);
for k=1:N_robot
    [tt{k,1},xx{k,1},QQ{k,1}] = generate_backup_traj(x(:,k),con1,stop_crit);
    [tt{k,2},xx{k,2},QQ{k,2}] = generate_backup_traj(x(:,k),con2,stop_crit);
    [tt{k,3},xx{k,3},QQ{k,3}] = generate_backup_traj(x(:,k),con3,stop_crit);
    
end

for k=1:N_robot
    uu = [con1(x(:,k)) con2(x(:,k)) con3(x(:,k))];
    if mod(t,10)==0
        bsc = bs(k);
    else
        bsc = 1:3;
    end
    for l =1:length(bsc)
        m = bsc(l);
        Aineq = [];
        bineq = [];
        [f,g] = dubin(x(:,k));
        h = 1e-3;
        db = zeros(4,1);
        
        for j = 1:length(tt{k,m})
            for n = 1:N_robot
                if n~=k
                    if j<=length(tt{n,bs(n)})
                        
                        [b,db] = two_robot_XY(xx{k,m}(j,:)',xx{n,bs(n)}(j,:)',r);
                        if b<-0.02
                            disp('')
                        end
                        if b<0.2
                            if x(3,k)>0.1
                                disp('')
                            end
                            
                            
                            Aineq = [Aineq;-db'*QQ{k,m}(4*(j-1)+1:4*j,:)*g];
                            bineq = [bineq;alpha*(b-0.01)+db'*QQ{k,m}(4*(j-1)+1:4*j,:)*f];
                        end
                    else
                        
                        [b,db] = two_robot_XY(xx{k,m}(j,:)',xx{n,bs(n)}(end,:)',r);
                        if b<-0.02
                            disp('')
                        end
                        if b<0.2
                            Aineq = [Aineq;-db'*QQ{k,m}(4*(j-1)+1:4*j,:)*g];
                            bineq = [bineq;alpha*(b-0.01)+db'*QQ{k,m}(4*(j-1)+1:4*j,:)*f];
                            if all(Aineq(end,:)==0) && bineq(end)<0
                                disp('')
                            end
                        end
                    end
                end
            end
            [b,db] = bdry_obs_XY(xx{k,m}(j,:)',bdry);
            if b<0.2
                Aineq = [Aineq;-db'*QQ{k}(4*(j-1)+1:4*j,:)*g];
                bineq = [bineq;alpha*b+db'*QQ{k}(4*(j-1)+1:4*j,:)*f];
                if all(Aineq(end,:)==0) && bineq(end)<0
                    disp('')
                end
            end
            
        end
        if ~isempty(Aineq)
            idx = find(all(Aineq'==0)');
            Aineq(idx,:) = [];
            bineq(idx,:) = [];
            Aineq = [Aineq -ones(length(bineq),1)];
        end
        
        if isempty(Aineq)
            u(:,k) = u0(:,k);
            backup_state(k)=m;
            break
        else
            [res,fval,exitflag] = quadprog(H,[-u0(:,k);1e3],Aineq,bineq,[],[],[-am;-rm;0],[am;rm;inf],[u0(:,k);0],options);
            if length(bsc)==1
                u(:,k) = res(1:2);
                backup_state(k)=m;
            else
                uu(:,m) = res(1:2);
            end
            
        end
        
    end
    if length(bsc)==3&&m==3
        interv_norm = vecnorm(uu-u0(:,k));
        [~,idx]=min(interv_norm);
        backup_state(k)=idx;
        u(:,k)=uu(:,idx);
        if x(3,k)<-vm
            u(1,k)=max(u(1,k),-x(3,k)-vm);
        elseif x(3,k)>vm
            u(1,k)=min(u(1,k),vm-x(3,k));
        end
    end
end

end

function [tt,xx,QQ] = symmetry_backup(x,tt0,xx0,QQ0)
if x(3)<0
    x(4)=round_2pi(x(4)+pi);
    x(3)=-x(3);
    [tt,xx,QQ] = symmetry_backup(x,tt0,xx0,QQ0);
else
    ts = tt0(2)-tt0(1);
    v = x(3);
    theta = x(4);
    rot = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    N = size(xx0,1);
    idx = min(find(xx0(:,3)<=v));
    xx_xy = xx0(idx:end,1:2)*rot'+x(1:2)';
    xx = [xx_xy xx0(idx:end,3) xx0(idx:end,4)+theta];
    QQ = [];
    for i=1:N-idx+1
        dxy_dvtheta = rot*QQ0([4*(idx+i-2)+1,4*(idx+i-2)+2],[3,4]);
        QQ = [QQ;[eye(2) dxy_dvtheta;zeros(2,2),eye(2)]];
    end
    tt = ts*(1:N-idx+1)';
end

end

function b = bdry_obs(x,bdry)
global am
ts = 0.03;

v = x(3);
t = 0:ts:abs(v)/am;
a = -am*sign(v)*[cos(x(4));sin(x(4))];
x_traj =x(1:2)+v*t+0.5*a*t.^2;
b = min([x_traj(1,:)-bdry(1),bdry(2)-x_traj(1,:),x_traj(2,:)-bdry(3),bdry(4)-x_traj(2,:)]);

end

function [b,db] = bdry_obs_XY(x,bdry)
[b,idx] = min([x(1)-bdry(1),bdry(2)-x(1),x(2)-bdry(3),bdry(4)-x(2)]);
switch idx
    case 1
        db = [1;0;0;0];
    case 2
        db = [-1;0;0;0];
    case 3
        db = [0;1;0;0];
    case 4
        db = [0;-1;0;0];
end
end
function [b,db] = two_robot_XY(x1,x2,r)
b0 = norm(x1(1:2)-x2(1:2));
db = [(x1(1:2)-x2(1:2))/b0;0;0];
b = b0-r;
end


function b = two_robot_obs(x1,x2,r)
global rm
R1 = x1(3)/rm;
R2 = x2(3)/rm;
c11 = x1(1:2)+R1*[sin(x1(4));-cos(x1(4))];
c12 = x1(1:2)+R1*[-sin(x1(4));cos(x1(4))];

c21 = x2(1:2)+R2*[sin(x2(4));-cos(x2(4))];
c22 = x2(1:2)+R2*[-sin(x2(4));cos(x2(4))];
b = max([norm(c11-c21),norm(c11-c22),norm(c12-c21),norm(c21-c22)])-R1-R2-r;
end

function b = two_robot_obs_brake(x1,x2,r)
global am
ts = 0.01;
if x1(3)>x2(3)
    temp = x1;
    x1 = x2;
    x2 = temp;
end

v1 = x1(3)*[cos(x1(4));sin(x1(4))];
v2 = x2(3)*[cos(x2(4));sin(x2(4))];
a1 = -am*[cos(x1(4));sin(x1(4))];
a2 = -am*[cos(x2(4));sin(x2(4))];
t1 = 0:ts:max(x1(3)/am,0);
t2 = 0:ts:max(x2(3)/am,0);
x1_traj =x1(1:2)+v1*t1+0.5*a1*t1.^2;
x2_traj =x2(1:2)+v2*t2+0.5*a2*t2.^2;
delta_t = length(t2)-length(t1);
x1_traj = [x1_traj x1_traj(:,end)*ones(1,delta_t)];
[dmin,idx]=min(vecnorm(x1_traj - x2_traj));
if idx==1 || idx ==size(x1_traj,2)
    b = dmin-r;
else
    x1_bar = x1_traj(:,idx-1);
    x2_bar = x2_traj(:,idx-1);
    t = t2(idx-1);
    v2_bar = x2(3)-a2*t;
    if t<t1(end)
        v1_bar = v1-a1*t;
    else
        v1_bar = 0;
    end
    ts_bar = 1e-4;
    t1_bar = 0:ts_bar:min(norm(v1_bar)/am,2*ts);
    t2_bar = 0:ts_bar:2*ts;
    x1_traj_bar =x1_bar+v1_bar*t1_bar+0.5*a1*t1_bar.^2;
    x2_traj_bar =x2_bar+v2_bar*t2_bar+0.5*a2*t2_bar.^2;
    delta_t = length(t2_bar)-length(t1_bar);
    x1_traj_bar = [x1_traj_bar x1_traj_bar(:,end)*ones(1,delta_t)];
    b = min(vecnorm(x1_traj_bar - x2_traj_bar))-r;
end
end



function [f,g] = dubin(x)
global vm rm
f = [x(3)*cos(x(4));x(3)*sin(x(4));0;0];
g = [zeros(2,2);eye(2)];

end

function ja = dubin_f_x(x,con)
h = 1e-6;
dudt = zeros(2,4);
dudt(:,1) = (con(x+[h;0;0;0])-con(x-[h;0;0;0]))/2/h;
dudt(:,2) = (con(x+[0;h;0;0])-con(x-[0;h;0;0]))/2/h;
dudt(:,3) = (con(x+[0;0;h;0])-con(x-[0;0;h;0]))/2/h;
dudt(:,4) = (con(x+[0;0;0;h])-con(x-[0;0;0;h]))/2/h;
ja = [0 0 cos(x(4)) -x(3)*sin(x(4));...
    0 0 sin(x(4)) x(3)*cos(x(4));...
    dudt];
end

function [tt,xx,QQ] = generate_backup_traj(x,con,stop_criterion,ts)
global am
if nargin<4
    ts = 0.05;
end
global vm
t = 0;
tt = 0;
xx = x';
Q = eye(4);
QQ = Q;
while ~stop_criterion(t,x)
    
    
    u = con(x);
    [f,g] = dubin(x);
    xdot = f+g*u;
    if x(3)<-vm
        xdot(3)=max(u(1),-x(3)-vm);
    elseif x(3)>vm
        xdot(3)=min(u(1),vm-x(3));
    end
    
    ja = dubin_f_x(x,con);
    x = x + xdot*ts;
    if abs(x(3))<ts*am/2
        x(3)=0;
    end
    Q = Q + ja*Q*ts;
    t = t + ts;
    tt = [tt;t];
    xx = [xx;x'];
    QQ = [QQ;Q];
    
end
end

