cc = zeros(50,1);
% cc(1) = 20*2*pi;
% cc(26) = 20*2*pi;
cc(37) = 0.5;
cc(38) = -0.5;
% cc(2:25) = randn(1, 24);
[q1, q2]= CFprint(cc, 2, 1, 'off');

t1=0:0.13:2*pi;t2=0:0.13:2*pi;
f1=figure;

subplot(1,2,1);surf(t1,t2,q1','FaceColor','interp');
view([-40 50])
set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
xlabel('\phi_1');ylabel('\phi_2');zlabel('q_1(\phi_1,\phi_2)');axis tight

subplot(1,2,2);surf(t1,t2,q2','FaceColor','interp');
view([-40 50])
set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
xlabel('\phi_1');ylabel('\phi_2');zlabel('q_2(\phi_1,\phi_2)');axis tight
zlim([-1.9 1.9])

colormap(summer)
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
shading interp
title("$-1 \sin(\phi_{1} - \phi_{2})$", 'Interpreter', 'latex', 'Color', [0.0157    0.5098    0.4000])

cc = zeros(50,1);
cc(11) = 0;
cc(37) = 1;
cc(38) = -1;
[q1, q2]= CFprint(cc, 2, 1, 'off');

t1=0:0.13:2*pi;t2=0:0.13:2*pi;
f1=figure;

subplot(1,2,1);surf(t1,t2,q1','FaceColor','interp');
view([-40 50])
set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
xlabel('\phi_1');ylabel('\phi_2');zlabel('q_1(\phi_1,\phi_2)');axis tight
title("$1 \sin(\phi_{1} - \phi_{2})$", 'Interpreter', 'latex', 'Color', [0.0157    0.5098    0.4000])

subplot(1,2,2);surf(t1,t2,q2','FaceColor','interp');
view([-40 50])
set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
xlabel('\phi_1');ylabel('\phi_2');zlabel('q_2(\phi_1,\phi_2)');axis tight
zlim([-2.25 2.25])

colormap(summer)
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
shading interp
title("$\sin(\phi_{1} - \phi_{2}) - \cos(\phi_{1} - \phi_{2})$", 'Interpreter', 'latex', 'Color', [0.0157    0.5098    0.4000])