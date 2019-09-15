%% Author: Dr. Feng Tian
%-------------------------------------------------------------------------%

function drawsurf(x,t,u,wtitle,ftitle)
figure('Name',wtitle)
surface(x,t,u,'LineStyle','None');
view([-40,25])
zlim([-3,3])
title(ftitle);
xlabel('x')
ylabel('t')
grid on
colorbar
rotate3d
end
