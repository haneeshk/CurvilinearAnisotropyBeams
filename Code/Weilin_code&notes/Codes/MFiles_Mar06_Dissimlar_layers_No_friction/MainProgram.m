clc;
close all;
clear all;

% Please input the total # of layers
M=input('Please input the # of layers: ');
flag=input('Homogeneous/Heterogeneous materials(1/2): ');
% Please choose 1: uniform thickness or 2: geometric sequence of thickness
% flag=input('Please choose 1: uniform thickness or 2: geometric sequence of thickness:');

delete('A.dat','B.dat','X.dat','Residual_norm.dat')

a=2e-0; b=14e-0;% inner and outer radii of the cylinder, unit: mm

    EI=zeros(M,3);
    for I=2:M+1
        radiiVec=linspace(a, b, I)'; % uniform thickness of each layer
        EI(I-1,1)=Jolicoeur(radiiVec(:,1),I,flag)*1e-6; % bending stiffness
        % the unit of compliance matrix cij is mm^2/N, need to time 1e-6
        % to convert the unit of EI to N*m^2
    end
    
%    EIinf=Infinity(a,b);
    
%     figure(1)
%     plot(radiiVec,0,'k+','markers',5);
%     ylim([-0.1 0.1]);
%     hold on;
% 
%     q=1.1; % layer thickness ratio
%     for I=1:M
%         radiiVec=zeros(I+1,1);
%         l=(b-a)*(1-q)/(1-q^I);
%         radiiVec(I+1)=b;
%         for j=I:-1:1
%             radiiVec(j)=radiiVec(j+1)-l*q^(I-j);
%         end
%     EI(I,2)=Jolicoeur(radiiVec,I);
%     end
%     
%     plot(radiiVec,0.02,'r+','markers',5); 
%     hold on;
%     
%         q=0.9; % layer thickness ratio
%     for I=1:M
%         radiiVec=zeros(I+1,1);
%         l=(b-a)*(1-q)/(1-q^I);
%         radiiVec(I+1)=b;
%         for j=I:-1:1
%             radiiVec(j)=radiiVec(j+1)-l*q^(I-j);
%         end
%     EI(I,3)=Jolicoeur(radiiVec,I);
%     end
%     
%     plot(radiiVec,-0.02,'b+','markers',5); 
    

    figure(2)
    layernum=1:1:M;
    plot(layernum,EI(:,1),'-bo','markers',3,'LineWidth',1.5);
    hold on;
%     plot(layernum,EI(:,2),'-r*','markers',3,'LineWidth',1.5);
%     hold on;
%     plot(layernum,EI(:,3),'-b*','markers',3,'LineWidth',1.5);
%     hold on;
%    plot([1,M],[EIinf,EIinf],'--r','LineWidth',1.5);
%    hleg1=legend('uniform thickness-two materials','infinity # of layers- one material');
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
