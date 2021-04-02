clc;
close all;
clear all;

% Please input the total # of layers
M=input('please input the # of layers:');
% Please choose 1: uniform thickness or 2: geometric sequence of thickness
% flag=input('Please choose 1: uniform thickness or 2: geometric sequence of thickness:');

a=2e-6;
b=14e-6;

    EI=zeros(M,2);
    for I=2:M+1
        radiiVec=linspace(a, b, I)'; % uniform thickness of each layer
        EI(I-1,1)=Jolicoeur(radiiVec(:,1),I); % bending stiffness
    end

    figure
    plot(radiiVec,0,'r+','markers',5);
    ylim([-0.1 0.1]);
    hold on;

    q=1.000001; % layer thickness ratio
    for I=1:M
        radiiVec=zeros(I+1,1);
        l=(b-a)*(1-q)/(1-q^I);
        radiiVec(I+1)=b;
        for j=I:-1:1
            radiiVec(j)=radiiVec(j+1)-l*q^(I-j);
        end
    EI(I,2)=Jolicoeur(radiiVec,I);
    end
    
    plot(radiiVec,0.02,'r+','markers',5); 
    
    figure
    layernum=1:1:M;
    plot(layernum,EI(:,1),'-ro','markers',3,'LineWidth',1.5);
    hold on;
    plot(layernum,EI(:,2),'-b*','markers',3,'LineWidth',1.5);
    hleg1=legend('uniform sequence','geometric sequence');
    xlabel('Number of layers','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
    ylabel('Efftive bending stiffness ($\ensuremath{\, \mathrm{Nm^2}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
    grid on


% EI=zeros(M,1);
% if flag==1
%     for I=2:M+1
%         radiiVec=linspace(a, b, I)'; % uniform thickness of each layer
%         EI(I-1)=Jolicoeur(radiiVec,I); % bending stiffness
%     end
%     
% elseif flag==2
%     q=1.06; % layer thickness ratio
%     for I=1:M
%         radiiVec=zeros(I+1,1);
%         l=(b-a)*(1-q)/(1-q^I);
%         radiiVec(I+1)=b;
%         for j=I:-1:1
%             radiiVec(j)=radiiVec(j+1)-l*q^(I-j);
%         end
%     EI(I)=Jolicoeur(radiiVec,I);
%     end
%     
% else
%     disp('Illegal input!');
%     return;
% end

% layernum=1:1:M;
% 
% figure(1)
% plot(layernum,EI,'-ro','markers',3);
% xlabel('Number of layers','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
% ylabel('Efftive bending stiffness ($\ensuremath{\, \mathrm{Nm^2}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
% grid on
% 
% figure(2)
% plot(radiiVec,0,'r*','markers',3);
