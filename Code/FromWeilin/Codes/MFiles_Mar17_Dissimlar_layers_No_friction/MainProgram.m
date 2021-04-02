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
% 
if flag==2
    EIinf=zeros(3,1);
    for I=1:2
        EIinf(I)=Infinity(a,b,I,flag)*1e-6;
    end
        EIinf(3)=(EIinf(1)+EIinf(2))/2;
    figure(2)
    layernum=1:1:M;
    plot(layernum,EI(:,1),'-bo','markers',3,'LineWidth',1.5);
    hold on;
    plot([1,M],[EIinf(1),EIinf(1)],'--m','LineWidth',1);
    plot([1,M],[EIinf(2),EIinf(2)],'--k','LineWidth',1);   
    plot([1,M],[EIinf(3),EIinf(3)],'--r','LineWidth',1.5);
    hleg1=legend('uniform thickness - mat. 1&2','infinity # of layers - mat. 1',...
        'infinity # of layers - mat. 2','infinity # of layers - mat. 1&2');
    xlabel('Number of layers','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
    ylabel('Efftive bending stiffness ($\ensuremath{\, \mathrm{Nm^2}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
    grid on
    
elseif flag==1
    EIinf=Infinity(a,b,1,flag)*1e-6;
    figure(2)
    layernum=1:1:M;
    plot(layernum,EI(:,1),'-bo','markers',3,'LineWidth',1.5);
    hold on;
    plot([1,M],[EIinf,EIinf],'--r','LineWidth',1.5);
    hleg1=legend('uniform thickness','infinity # of layers');
    xlabel('Number of layers','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
    ylabel('Efftive bending stiffness ($\ensuremath{\, \mathrm{Nm^2}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
    grid on
    
end

