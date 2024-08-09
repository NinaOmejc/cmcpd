function [q1,q2]= CFprint(cc,bn, make_plots, visible)
%plots the coupling functions from the inferred parameters

%---inputs---
%cc - vector of inferred parameters
%bn - order of Fourier base function

%Note that the input is vector of parameters for one time window
%%
            %---evaluating the coupling functions -----
            t1=0:0.13:2*pi;t2=0:0.13:2*pi; 
            q1(1:length(t1),1:length(t1))=0;q2=q1;
            u=cc; K=length(u)/2;
            for i1=1:length(t1)                
                for j1=1:length(t2)
                    br=2;
                   
                    for ii=1:bn
                        q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t1(i1))+u(br+1)*cos(ii*t1(i1));
                        q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(ii*t2(j1))+u(K+br+1)*cos(ii*t2(j1));
                        br=br+2;  
                    end
                    for ii=1:bn
                        q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t2(j1))+u(br+1)*cos(ii*t2(j1));
                        q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(ii*t1(i1))+u(K+br+1)*cos(ii*t1(i1));
                        br=br+2;
                    end
                    
                    for ii=1:bn
                        for jj=1:bn                            
                                   q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t1(i1)+jj*t2(j1))+u(br+1)*cos(ii*t1(i1)+jj*t2(j1));                                                                
                                   q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(jj*t2(j1)+ii*t1(i1))+u(K+br+1)*cos(jj*t2(j1)+ii*t1(i1));                                   
                                   br=br+2;
                                   
                                   q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t1(i1)-jj*t2(j1))+u(br+1)*cos(ii*t1(i1)-jj*t2(j1));                                                                
                                   q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(jj*t2(j1)-ii*t1(i1))+u(K+br+1)*cos(jj*t2(j1)-ii*t1(i1));     
                                   br=br+2;
                        end
                    end                    
                                               
                end
            end
            if make_plots
                %---plotting -----
                if ~exist("visible", 'var'); visible = 'on'; end
                f1=figure('Visible', visible);

                subplot(1,2,1);surf(t1,t2,q1','FaceColor','interp');                                              
                view([-40 50])
                set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
                xlabel('\phi_1');ylabel('\phi_2');zlabel('q_1(\phi_1,\phi_2)');axis tight
                shading interp
                subplot(1,2,2);surf(t1,t2,q2','FaceColor','interp');                                               
                view([-40 50])
                set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')            
                xlabel('\phi_1');ylabel('\phi_2');zlabel('q_2(\phi_1,\phi_2)');axis tight
                
                colormap("summer") % copper
                set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
                shading interp
            end
 q1=q1';
 q2=q2;
                         %uncomment this lines for saving the figure
                          % saveas(f1,'filename','jpg');
                          % saveas(f1,'filename','fig');

%% 
