function [ data, g, data0 ] = Lane_keeping_HJ(accuracy)

tMax = 5.0;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';


W = 1.8;
vmax = 10;
thetamax = pi/6;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the 3D reachable set.
displayType = 'surface';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Create the grid.
clear g
g.dim = 3;
g.min = [ -2; 0; -1.5 ];
g.max = [ 2; 10; 1.5 ];
g.bdry = @addGhostExtrapolate;
g.N = [50,30,30]';
g = processGrid(g);
% rangeTx = [-L,L];
% rangeTy = [-W,W];
%---------------------------------------------------------------------------
% Create initial conditions.
%   A rectangle is the intersection of four halfplanes.
    y_grid=linspace(g.min(1),g.max(1),g.N(1));
    psi_grid=linspace(g.min(3),g.max(3),g.N(3));
    Z=zeros(g.N(1),1);
    for i=1:length(y_grid)
        Z(i)= W-abs(y_grid(i));
    end
    data0=zeros(g.N');
for i = 1:g.N(1)
    for j =1:g.N(2)
        for k=1:g.N(3)
          data0(i,j,k)=min(W-abs(y_grid(i)),pi/3-abs(psi_grid(k)));
        end
    end
end
data = data0;

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @lanekeepingHam;
schemeData.partialFunc = @LKPartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.rmax = 0.5;
schemeData.amax = 5;

% For evaluating the evader's speed, 
%   we might as well precompute the term min(sqrt(x^2 + y^2), S).

%---------------------------------------------------------------------------
% Choose degree of dissipation.

switch(dissType)
 case 'global'
  schemeData.dissFunc = @artificialDissipationGLF;
 case 'local'
  schemeData.dissFunc = @artificialDissipationLLF;
 case 'locallocal'
  schemeData.dissFunc = @artificialDissipationLLLF;
 otherwise
  error('Unknown dissipation function %s', dissFunc);
end

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.
%   The Lax-Friedrichs approximation scheme MUST already be completely set up.
innerFunc = schemeFunc;
innerData = schemeData;
clear schemeFunc schemeData;

% Wrap the true Hamiltonian inside the term approximation restriction routine.
schemeFunc = @termRestrictUpdate;
schemeData.innerFunc = innerFunc;
schemeData.innerData = innerData;
schemeData.positive = 0;

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);
xlabel('$Y$','interpreter','latex')
ylabel('$v$','interpreter','latex')
zlabel('$\psi$','interpreter','latex')

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Get correct figure, and remember its current view.
  figure(f);

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
  xlabel('$Y$','interpreter','latex')
  ylabel('$v$','interpreter','latex')
  zlabel('$\psi$','interpreter','latex')
end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = lanekeepingHam(t, data, deriv, schemeData)


checkStructureFields(schemeData, 'grid', 'rmax', 'amax');

grid = schemeData.grid;

% implements equation at the bottom of p.57 from my thesis term by term
%   with allowances for nonunit \script A and \script B
%   where deriv{i} is p_i
%         x is grid.xs{1}, y is grid.xs{2}
%         \script A is boundA and \script B is boundB
ax = schemeData.amax*sign(deriv{2}).*(grid.xs{2}>0);
r = schemeData.rmax*sign(deriv{3});
hamValue = -(deriv{1}.*grid.xs{2}.*sin(grid.xs{3}) +deriv{2}.*ax + deriv{3}.*r); 



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = ...
              LKPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
% acousticPartialFunc: Hamiltonian partial fcn for acoustic capture game.
%
% alpha = acousticPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype for the acoustic
%   capture game.  The Hamiltonian can be found in my PhD thesis
%   chapter 3.2, bottom of page 57.
%
% It calculates the extrema of the absolute value of the partials of the 
%   analytic Hamiltonian with respect to the costate (gradient).
%
% Parameters:
%   t            Time at beginning of timestep (ignored).
%   data         Data array.
%   derivMin	 Cell vector of minimum values of the costate (\grad \phi).
%   derivMax	 Cell vector of maximum values of the costate (\grad \phi).
%   schemeData	 A structure (see below).
%   dim          Dimension in which the partial derivatives is taken.
%
%   alpha	 Maximum absolute value of the partial of the Hamiltonian
%		   with respect to the costate in dimension dim for the 
%                  specified range of costate values (O&F equation 5.12).
%		   Note that alpha can (and should) be evaluated separately
%		   at each node of the grid.
%
% schemeData is a structure containing data specific to this Hamiltonian
%   For this function it contains the field(s):
%
%
%   .grid	 Grid structure.
%   .We          Speed of evader.
%   .Wp          Speed of pursuer.
%   .R           Turn radius of pursuer.
%   .boundA      Range of norm of input A (evader).
%   .boundB      Range of magnitude of input B (pursuer).
%   .speedBound  precomputed term min(sqrt(x^2 + y^2), S).
%
% Ian Mitchell 4/21/04

checkStructureFields(schemeData, 'grid', 'rmax', 'amax');

grid = schemeData.grid;



switch dim
  case 1
    alpha = abs(grid.xs{2});

  case 2 
    alpha = schemeData.amax;
  case 3
    alpha = schemeData.rmax;
  
  
  otherwise
    error([ 'Partials for the acoustic capture game' ...
            ' only exist in dimensions 1-2' ]);
end
