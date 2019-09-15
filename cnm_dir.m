%% Submitted by Debajyoti Basu Sarkar (DBS)
%  For MAT 695 (Fall 2017) offered by Dr. Feng Tian, Hampton University
%-------------------------------------------------------------------------%

%% API provided by Dr. Feng Tian
%-------------------------------------------------------------------------%

function  [u xx tt] = cnm_dir(pdefun,icfun,bcfun,x,t)
% The Crank-Nicolson method for the 1-D material balance equation 
% 
%     c(x,u) u_t = (f(x,u) u_x)_x + s(x,t), 
% 
% with initial condition 
% 
%     u(x,0) = icfun(x),
% 
% and boundary condition
% 
%     u(x0,t) = p(t),
%     u(x1,t) = q(t).
% 
% Inputs: 
%     pdefun{1} = c
%     pdefun{2} = f
%     pdefun{3} = s
%     icfun = u(x,0), initial condition
%     bcfun{1} = p, left boundary
%     bcfun{2} = q, right boundary
%     x: space grid (uniform)
%     t: time grid (uniform)
%
% outputs:
%     u: approximate solution; u(:,j) = approx values at time t(j). 
%     xx: space mesh grid
%     tt: time mesh grid


%% Initialization
%
%[xx tt] = ndgrid(x,t);
[xx, tt] = ndgrid(x,t); % MATLAB Suggestion

c = pdefun{1};
f = pdefun{2};
s = pdefun{3};
p = bcfun{1};
q = bcfun{2};

dx = x(2)-x(1);
dt = t(2)-t(1);
b  = dt/dx^2;

u = zeros(size(xx));
u(2:end-1,1) = reshape(icfun(x(2:end-1)),[],1);
u(1,:) = reshape(p(t),1,[]);
u(end,:) = reshape(q(t),1,[]);


%% Implementation (Insert code here.)
%

for j = 1:length(t)-1

%     wplus = weights(f(x(2:end),u(2:end,j)));
%     wminus = weights(f(x(1:end-1),u(1:end-1,j)));
    
%     f1 = f(x(:),u(:,j));
%     w = harmmean([f1(1:end-1), f1(2:end)],2);
    w = weights(f(x(:),u(:,j)));
    wplus = w(2:end-1);
    wminus = w(2:end-1);
    wmid = w(1:end-1)+w(2:end);
    
    L = diag(wminus,1) - diag(wmid) + diag(wplus,-1);
    Cmat = diag(c(x(2:end-1),u(2:end-1,j)));
    Scol = reshape( (s(x(2:end-1),t(j))+s(x(2:end-1),t(j+1)))*(dt/2),[],1);
    Scol(1) = Scol(1) + (b*w(1)/2)*(p(t(j))+p(t(j+1)));
    Scol(end) = Scol(end) + (b*w(end)/2)*(q(t(j))+q(t(j+1)));

    [u(2:end-1,j+1), ~] = bicgstab(Cmat-(b/2)*L,(Cmat+(b/2)*L)*u(2:end-1,j)+Scol,1e-8,1000);

end


end
