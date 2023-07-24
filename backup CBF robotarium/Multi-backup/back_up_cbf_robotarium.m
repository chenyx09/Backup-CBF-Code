%% Single-integrator Barrier Certificate Algorithm
% by Sean Wilson
% 7/2019

%% Set up Robotarium object
% Before starting the algorithm, we need to initialize the Robotarium
% object so that we can communicate with the agents


N = 3;
v = zeros(1,N);
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);
am = 0.2;
ts = Robotarium.time_step;
vm = Robotarium.max_linear_velocity;
u0_rec = [];
u_rec = [];
x_rec = [];

%This iis how many times the main loop will execute.
iterations = 3000;

%% Experiment constants
% Next, we set up some experiment constants

% Initialize velocity vector for agents.  Each agent expects a 2 x 1
% velocity vector containing the linear and angular velocity, respectively.
dx = zeros(2, N);

% This code ensures that the agents are initially distributed around an
% ellipse.
xybound = [-1, 1, -0.8, 0.8];
succ = 0;
while ~succ
    fail = 0;
    pose = [xybound(1)*ones(1,N)+(xybound(2)-xybound(1))*rand(1,N);xybound(3)*ones(1,N)+(xybound(4)-xybound(3))*rand(1,N)];
    for i =1:N-1
        for j=i+1:N
            if norm(pose(:,i)-pose(:,j))<0.2
                fail=1;
                break
            end
        end
        if fail
            break
        end
    end
    if ~fail
        succ = 1;
    end
end

r.initialize([pose;2*pi*rand(1,N)]);
p_theta = (1:2:2*N)/(2*N)*2*pi;
p_circ = [xybound(2)*cos(p_theta) xybound(2)*cos(p_theta+pi); xybound(4)*sin(p_theta)  xybound(4)*sin(p_theta+pi)];

x_goal = p_circ(:,1:N);

flag = 0; %flag of task completion

%% Retrieve tools for single-integrator -> unicycle mapping

% Let's retrieve some of the tools we'll need.  We would like a
% single-integrator position controller, a single-integrator barrier
% function, and a mapping from single-integrator to unicycle dynamics
position_control = create_si_position_controller();
si_to_uni_dyn = create_si_to_uni_dynamics();
uni_barrier_certificate = create_uni_barrier_certificate2();

%% Begin the experiment
% This section contains the actual implementation of the barrier
% certificate experiment.

%Iterate for the previously specified number of iterations
for t = 1:iterations
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();
    x_e = [x(1:2,:); v; x(3,:)];
    x_rec = [x_rec;x_e];
    x_temp = x(1:2,:);
    
    %% Algorithm
    
    % Let's make sure we're close enough the the goals
    if norm(x_goal-x_temp,1)<0.03
        flag = 1-flag;
    end
    
    % This code makes the robots switch positions on the ellipse
    if flag == 0
        x_goal = p_circ(:,1:N);
    else
        x_goal = p_circ(:,N+1:2*N);
    end
    
    % Use a single-integrator position controller to drive the agents to
    % the circular formation
    togoal = zeros(2,N);
    adj = zeros(N,N);
    for i=1:N
        togoal(:,i) = x_goal(:,i)-x(1:2,i);
        for j=1:N
            toobs  = x(1:2,j)-x(1:2,i);
            if j~=i && norm(toobs)<0.3 && abs(atan2(togoal(1,i)*toobs(2)-toobs(1)*togoal(2,i),togoal(1,i)*toobs(1)+togoal(2,i)*toobs(2)))<pi/3
                adj(i,j)=1;
            end
        end
    end
    x_goal1 = x_goal;
    for i=1:N
        if abs(v(i))<0.05 && sum(adj(i,:))>0 && norm(togoal(:,i))>0.5
            delta = [x_goal(2,i)-x(2,i);x(1,i)-x_goal(1,i)];
            delta = delta/norm(delta)*0.4;
            wp1 = x(1:2,i)+delta;
            wp2 = x(1:2,i)-delta;
            min_dis1 = vecnorm(x(1:2,adj(i,:)==1)-wp1);
            min_dis2 = vecnorm(x(1:2,adj(i,:)==1)-wp2);
            if min_dis1>min_dis2
                x_goal1(:,i)=wp1;
            else
                x_goal1(:,i)=wp2;
            end
        end
    end
    %     for i=1:N
    %         if abs(v(i))<0.01 && sum(adj(i,:))>0 && norm(togoal(:,i))>0.5
    %             collective_delta = x(1:2,i)-sum(x(1:2,[1:i-1,i+1:end]),2)/(N-1);
    %             x_goal1(:,i)= x(1:2,i)+collective_delta/norm(collective_delta)*0.6;
    %         end
    %     end
    
    dx = position_control(x(1:2, :), x_goal1);
    
    %% Apply barrier certs. and map to unicycle dynamics
    % Transform the single-integrator dynamics to unicycle dynamics using a
    % diffeomorphism, which can be found in the utilities
    dx = si_to_uni_dyn(dx, x);
    dx(1,:) = min(max(dx(1,:),v-am*ts),v+am*ts);
    dx(1,:) = min(max(dx(1,:),-vm),vm);
    u0 = [(dx(1,:)-v)/ts;dx(2,:)];
    u0_rec = [u0_rec;u0];
    %Ensure the robots don't collide
    u = backup_CBF_con(x_e,u0,xybound);
    u_rec = [u_rec;u];
    dx = [v+u(1,:)*ts;u(2,:)];
    v = dx(1,:);
    for i=1:N
        if abs(dx(1,i))+r.base_length*abs(dx(2,i))>0.2
            dx(2,i)=sign(dx(2,i))*(0.199-abs(dx(1,i)))/r.base_length;
        end
    end
    
    %% Set the velocities of the agents
    % Set velocities of agents 1,...,N
    r.set_velocities(1:N, dx);
    
    % Send the previously set velocities to the agents.  This function must be called!
    r.step();
end
save('result_log','x_rec','u0_rec','u_rec');
% Print out any simulation problems that will produce implementation
%differences and potential submission rejection.
r.debug()