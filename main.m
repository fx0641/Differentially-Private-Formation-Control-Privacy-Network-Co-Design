clc; clear; close all;
%% Inputs
N=10;
eps_min = [0.0151 0.0505 0.1149 0.2254 0.2340 0.1209 0.1628 0.0913 0.0043 0.1313]';
eps_max = [0.6812 0.8625 0.7827 0.5163 0.7130 0.3870 0.6376 0.4721 0.7256 0.5130]';
max_error = 50;
gamma = 1/2/N;
delta = 0.05;
adjparam = 1;
edge_budget = 50;
d=1;

givenL =[4     0     0     0    -1    -1    -1    -1     0     0;
     0     2    -1     0    -1     0     0     0     0     0;
     0    -1     5     0    -1     0    -1    -1    -1     0;
     0     0     0     1     0     0     0    -1     0     0;
    -1    -1    -1     0     4     0     0    -1     0     0;
    -1     0     0     0     0     2     0     0    -1     0;
    -1     0    -1     0     0     0     3     0     0    -1;
    -1     0    -1    -1    -1     0     0     4     0     0;
     0     0    -1     0     0    -1     0     0     2     0;
     0     0     0     0     0     0    -1     0     0     1];


gammai = [0.0625 0.1482 0.1576 0.7350 0.3236 0.4521 0.4035 0.0662 0.3338 0.1290]';
taui = [0.0415 0.0081 0.3451 0.1763 0.0226 0.1746 0.2877 0.2722 0.2138 0.1952]';
rhoi = [14.7942 15.2751 15.0965 15.1682 14.8436 14.7513 14.8464 14.9536 14.7831 15.2329]';

%% Initialize fmincon
lb = [];
ub = [];

A = [];
b = [];
Aeq = [];
beq = [];

initialL = -abs(randn(N*(N-1)/2,1));
initialeps = .1*ones(N,1);
x0 = [initialL' initialeps']';

%% Implement constraints and solve

nonlcon = @(x) constraints(x,N,eps_min, eps_max, max_error, gamma, delta, adjparam, edge_budget, gammai, taui, rhoi, givenL, d);
options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 1000000;
options.MaxIterations =10000;
options.StepTolerance =1e-100;
%options.ConstraintTolerance =1e-10;

x = fmincon(@(x)objectiveFunction(x,N),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
epsilon = x(N*(N-1)/2+1:end)
L = makeLaplacian(x,N)

%% Plot it
%get rid of the edges that shouldn't be there
[m,n] = size(L);
for i = 1:m
   for j = 1:n
       if abs(L(i,j)) <= .0001
           L(i,j) = 0;
        end
   end
end
% arrange agents in a circle
angles = linspace(0,2*pi,N+1);
r = 1;
xs = r*cos(angles);
ys = r*sin(angles);

initialG = graph(abs(givenL),'omitselfloops');
subplot(1,2,1)
p1 = plot(initialG);
p1.XData = xs(1:end-1);
p1.YData = ys(1:end-1);
title('Input Graph')

G = graph(abs(L),'omitselfloops');
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
subplot(1,2,2)
p2 = plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths);
p2.XData = xs(1:end-1);
p2.YData = ys(1:end-1);
title('Optimized Graph')
