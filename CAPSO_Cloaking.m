%This is a CAPSO variant of the APSO algorithm, specifically adjusted to a
%specific electromagnetic coating problem. More information regarding the 
%acquisition of the specific optimization problem can be found in the
%README file.

%The code is strongly based upon Yang's work which is referenced below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] Gandomi, A.H.; Yun, G.J.; Yang, X.S.; Talatahari, S. Chaos-enhanced
% accelerated particle swarm optimization.%
% Communications inNonlinear Science and Numerical Simulation2013,18,
% 327–340.%
% Files of the Matlab programs for APSO included in the book:       %
% [2] Xin-She Yang, Nature-Inspired Metaheuristic Algorithms,  %
% Second Edition, Luniver Press, (2010).   www.luniver.com %
% or can be found at:
% XS Yang (2021). Accelerated Particle Swarm Optimization (APSO)           
%(https://www.mathworks.com/matlabcentral/fileexchange/74766-accelerated-particle-swarm-optimization-apso),
% MATLAB Central File Exchange.
% A good point of reference for the basics of PSO.
% [3] Rao, S. S. (2019). Engineering optimization: theory and practice. John Wiley & Sons.
% ======================================================== %
function [fmin,]=capso_ndim_demo(n, Num_iterations)
clear all;
clc;
% n=number of particles
% Num_iterations=total number of iterations
% default = (n=40, Num_iterations = 250)
if nargin<2,   Num_iterations=500;  end
if nargin<1,   n=180;          end

%%Seach domains settings
%N = Number of coating layers.
N = 4;
%nd =number of dimensions - this can be used if the default objective
%function is changed.
%nd for this problem refers to the number of layers, and their relative magnetic permeabilities and dielectric permittivities
nd = 3*(N-1); 

%material boundaries (this regards \epsilon and \mu, permittivities and
%permeabilities. They correspond to physical constraints.
%umaterial = upper material limit
%lmaterial = lower material limit
umaterial = 10;
lmaterial = 1;

%boundary vectors for materials:
lb = lmaterial*ones(1,2*(N-1));
ub = umaterial*ones(1, 2*(N-1));
%experiment cycle settings
cycles = 5;

% final = result table, here we initialize it
final = zeros(cycles,3.*(N-1) +1);
%--------------------------------------------------------------------
% Setting the parameters: alpha, beta
beta=0.7; %for sinusoidal_s we need b= 0.7, for most experiments it remained like this
alpha=1; %alpha= X can be various numbers around 1 to 10 for this experiment in this format of a, in general we need to scale to a problem though
gamma=(10^(-20)/alpha)^(1/Num_iterations); %a number smaller than 1 so it can decrease if raised in a power

% Start Chaotic Accelerated Particle Swarm Optimization ---------------------
% generating the initial locations of n particles
%initialize the swarm

%Start Experiment Loop -----------------------------------------------
for cycle = 1:cycles
    %optional messages
    disp('Cycle number');
    disp(cycle);
    cxo =[]; 
    %It is possible to work with both min/max, it depends on how you handle
    %the obj.function. The min/max default values also depend on the
    %problem.

    %initial swarm of this cycle
    [xn]=initial_pso(n,lb, ub, N);
    %disp("Initial randomized swarm for this cycle");
    %disp(xn);
    i = 1;
    sflag = false;
    while (i <= Num_iterations)&&(sflag == false)
        %find the current g* - global best, be it minimum or maximum
        %[fmin,xo]=findmin(xn,N);
        [fmax,xo] = findmax(xn,N);
        %g = [xo fmin];
        %display the global best if you want
        g = [xo fmax];
        %store it in a best array if you want history
        %disp(g);
        
        %alpha can be handled with many different options [1][2], here it is a simple decreasing function
      
        alpha=gamma.^i;
        %beta chaotic map options. Singer and Sinusoidal(s) were used for
        %experimentation.
        %beta = singer(beta);
        %beta = sinusoidal_s(beta); %note, sinusoidal s is very easily developped and produces good noise. In general the majority of the best designs were produced with sinusoidal. In bibliography both sinusoidal(s) and singer did really well
        beta = sinusoidal(beta);
        
        % Move all the particles to new locations
        [xn]=pso_move(xn,xo,alpha,beta,lb,ub,N);
               
        %optional display messages
        %disp('new positions');
        %disp(xn);
            
        %optional displays
        %display the best fmin every 100 iterations
        %if round(i/100)==i/100
            %best(i,1:end)
        %nd
        
        %optional convergence checks through swarm position standard
        %deviation
        S = std(xn);
        if S(1,:) < 0.000001
            disp("This cycle closes due to reaching standard deviation standards as convergence");
            sflag = true;
            final(cycle,1:nd)= cxo; final(cycle,end)=fobj_expt(cxo,N); 
        end
        %check if reached the end of iterations for this cycle
        if i==Num_iterations
            final(cycle,1:nd)= xo; final(cycle,end)=fobj_expt(xo,N); 
        end
    i = i +1;
    end   %%%%% end of iterations
    %
end
%display simulation/ experiment table, it displays the final global best per experiment
disp(final);
end


% ----- All subfunctions are listed here -----

%Initialize radii seperately, due to their specific conditions and
%accretive sizing
function [rad] = initialize_rad(N)
%Initialize radii according to these boundaries: [2π/20,π/2]
%Double floats as values
rad = zeros(1,N-1);
%mysum refers to the core's radius, we ran multiple experiments.
%mysum = 2.*pi;
mysum = pi;
for i=1:(N-1)
    myrand = (2.*pi)./20 + ((pi)./2 -(2.*pi)./20).*rand(1,1,'double');
    rad(N-i)= mysum + myrand;
    mysum= rad(N-i);
end
end
% Move all the particles toward the  xo
function [xn]=pso_move(xn,xo,alpha,beta,Lb,Ub,N)
nd=size(xn,2);  
for j=1:size(xn,1)
    xn(j,:)=xn(j,:).*(1-beta)+xo.*beta+alpha.*randn(1,nd);
    % Check if the new solution is within limits,boundary checking
    xn(j,:)=simplebounds(xn(j,:),Lb,Ub,N);
end
end
% Application of constraints or bounds
function s =simplebounds(s,Lb,Ub,N)
  %The strategy is to ensure whether all variables of the coating problem keep to
  %their physical boundaries. In case they go out of bounds, they are
  %checked and re-randomized wrt the restrictions.
  %Apply boundaries on radii prorities 
  flag = false;
  ns_tmp=s;
  r = ns_tmp(1,1:(N-1));
  %if (((r(1)<r(2))||(r(2)<r(3)))||(r(3)<(2*pi)))
      %r = initialize_rad(N);
  %end
  for i=1:(N-2)
      if (r(1,i)<=r(1,i+1))
          flag = true;
      end
  end
  if r(1,end)<= (pi) %(2.pi) <=VERY IMPORTANT, this has to refer to the core we chose.
      flag = true;
  end
  if flag
      r = initialize_rad(N);
  end
  
  %Apply the lower bound - materials e, μ
  materials = ns_tmp(1,N:end);
  I=materials<Lb;
  materials(I)=Lb(I);
  % Apply the upper bounds - materials e, μ
  J=materials>Ub;
  materials(J)=Ub(J);
  % Update this new move by forming a new xn restricted to boundaries and
  % physical limitations
  ns_tmp = [r materials];
  s=ns_tmp;
end

% Initial locations of n particles - generic benchmark/ artificial landscape version
function [xn]=init_pso(n,Lb,Ub)
% Uniform sampling of the search space for the initial population
for j=1:n
    xn(j,:)=Lb+(Ub-Lb).*rand(size(Lb));
end
end

% Initial locations of n particles - ELECTROMAGNETIC CLOAKING PROBLEM.
function [xn]=initial_pso(n,Lb,Ub,N)
% Uniform sampling of the search space for the initial population
for j=1:n
    xn(j,:)=(Lb+(Ub-Lb).*rand(size(Lb)));
end
for j = 1:n
    an(j,:) = initialize_rad(N);
end
xn = [an xn];
end

% Find the  solution xo and its fmin
function [fmin,xo]=findmin(xn,N)
fmin=10^10;
for j=1:size(xn,1)
   fnew=fobj_expt(xn(j,:),N);
   if fnew<fmin
       fmin=fnew;
       xo=xn(j,:);
   end
end
end
% Find the  solution xo and its fmax
function [fmax,xo]=findmax(xn,N)
fmax=-10^10;
for j=1:size(xn,1)
   fnew=fobj_exp(xn(j,:),N);
   if fnew>fmax
       fmax=fnew;
       xo=xn(j,:);
   end
end
end
%Chaotic Maps

%alternative sinusoidal, works with specific a, x0 [1]
function [update] = sinusoidal_s(beta)
%x0 = 0.7 (b = 0.7), a = 2.3
update = sin(pi.*beta);
end
%Sinusoidal Map
function [update] = sinusoidal(beta)
%set a as various nums but 2.3 is good
a = 2.3
update = (a).*(beta.^2).*sin(pi.*beta);
end
%Singer Map
function [update] = singer(beta)
%μ is subjected to change, [0.9, 1.08], picked 1.07 due to loosely related bibliography
m = 0.908; %small ms seemed to fare well
update = m.*(7.86.*beta - 23.31.*(beta.^2) + 28.75.*(beta.^3) - 13.302875.*(beta.^4));
end

% Objective Functions You can replace the following by your own functions
% and readjust the code. Please refer to [2] for APSO, or extra clarity.
% A d-dimensional (nd) objective function

%example function, Sphere
function z=fobj(u)
% The D-dimensional sphere function sum_j=1^D (u_j-1)^2. 
%  with the minimum fmin=0 at (1,1, ...., 1); 
z=sum((u - 1).^2);
end

%Electromagnetic cloaking optimization problem, the experiment's
%obj.function. This version is used for evaluation of solutions, while
%running the algorithm.
function z = fobj_exp(u,N)
%optics experiment
a = zeros(1,N);
for i=1:N-1
    a(i)=u(i);
end
%alternative core radii
%a(N)=2.*pi;
a(N)= pi;
%alternative distances of light source
%r0 = (x).*a(1);
r0 = 10.*a(N);
%
%
er = u(N:2*(N-1));
mr = u((2*(N-1)+1):3*(N-1));
%
%core permit.:
er(N)= 2.1;
%er(N)= 18.7;

%core permeab.:
mr(N)= 1;
%ind = 'PEC' or 'dil' depending on the type of core we utilize.
ind = 'dil';
%for access to the files of this optimization problem model, please refer
%to the README file
flag = objective_function_rcs_total_sph_new(ind,r0,a,er,mr);
%This obj.function is certain to be always positive in result, so we use
%this trick to keep treating it as a maximization problem, [3]
z = 1/(1+flag);
end

%this version of the obj.function provides the logarithmized results for
%reasons of clarity. The results matrix stores the results in this form.
function z = fobj_expt(u,N)
%optics experiment
a = zeros(1,N);
for i=1:N-1
    a(i)=u(i);
end
%core radii
%a(N)=2.*pi;
a(N)= pi;
%source of light
%r0 = (x).*a(1);
r0 = 10.*a(N);
%
er = u(N:2*(N-1));
mr = u((2*(N-1)+1):3*(N-1));
%
%core permit./permab.
er(N)= 2.1;
%er(N)= 18.7;
%
mr(N)=1;
%core type indication
ind = 'dil';
flag = objective_function_rcs_total_sph_new(ind,r0,a,er,mr);
t_res = 10.*log10(flag);
z = t_res;
end
