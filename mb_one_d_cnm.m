%% Submitted by Debajyoti Basu Sarkar (DBS)
%  For MAT 695 (Fall 2017) offered by Dr. Feng Tian, Hampton University
%-------------------------------------------------------------------------%

%% API provided by Dr. Feng Tian
%  Author: Dr. Feng Tian
%-------------------------------------------------------------------------%

%% MAT 695-01 Exercise 5
% Application of the Crank-Nicolson method (CNM) to simulate the 1D
% material balance equation
% 
%     c(x,u) u_t = (f(x,u) u_x)_x + s(x,t), 
% 
% where
% 
%     c(x,u) = 2+cos(x),
%     f(x,u) = exp(x.*u),
%     s(x,t) = 10.*x.*(1-x)*t, 
%     x in [0,1], 
%     t in [0,2];
%  
% and subject to the following conditions: 
% 
%     u(x,0) = sin(2*pi*x),
%     u(0,t) = 2*sin(pi*t),
%     u(1,t) = 1-cos(pi*t).
%     
%
% The algorithm is implemented in cnm_dir.m. (Do not change this file.)

clear all
close all
clc




%% 1. Illustration
% Figure 1 shows the approximate solution to the initial-boundary value 
% problem.

load data
drawsurf(X,T,soln,'Solution', 'u(x,t)');




%% 2. Initialization

x = 0:0.01:1;
t = 0:0.005:2;

pdefun = {
    @(x,u) 2+cos(x),...
    @(x,u) exp(x.*u),...  
    @(x,t) 10*x.*(1-x)*t...
};

icfun = @(x) sin(2*pi*x);

bcfun = {
    @(t) 2*sin(pi*t),...
    @(t) 1-cos(pi*t)...
};




%% 3. Generate the approximate solution.

[u, xx, tt] = cnm_dir(pdefun,icfun,bcfun,x,t);
drawsurf(xx,tt,u,'CNM Approximation', 'u(x,t)');




%% 4. Approximation Error.

format long
fprintf('\nThe maximum absolute error (MAE) should be 0.005652.\n\n');
fprintf('Your result: %f\n',max(abs(u(:)-soln(:))));

figure
plot(t,max(abs(u-soln)))
xlabel('t')
ylabel('MAE(t)')
title('max |u(t) - soln(t)|')



