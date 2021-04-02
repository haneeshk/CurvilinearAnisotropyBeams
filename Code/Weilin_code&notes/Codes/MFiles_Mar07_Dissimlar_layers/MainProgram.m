clc;
clear all;
close all;
% Please input the total # of layers
display('Note: Need to specify layer #, homo/hetrogeneous layers')
display('      && elastic symmetries, interface conditions.');
display('For heterogeneous case, only consider orthotropic matrials.');
display('-----------------------Enter info-----------------------');

M=input('Please input the layer number: ');
% Please indicate material cases
flag_mat=input('Homogeneous/Heterogeneous materials(1/2): ');
% Please specify elastic sysmmetry
flag_sym=input('Orthontropic/TransverseIso/Isotropic(1/2/3): ');
% Please specify the interfacial condition
flag_intf=input('No friction/No slip(1/2): ');

delete('A.dat','B.dat','X.dat','Residual_norm.dat')

a=2e-0; b=14e-0;% inner and outer radii of the cylinder, unit: mm

% The angle between material principle direction and z-axis of the cylinder
% angle = pi/2 - (helix angle), Currently only applies to homogenous cases

 angle=(-90:15:90)/180*pi;
% angle=-15/180*pi;
length=length(angle);

EI=zeros(M,length);
EIinf=zeros(1,length);

for k=1:length
    
    alpha=angle(k);
    
    for I=2:M+1
        radiiVec=linspace(a, b, I)'; % uniform thickness of each layer
        EI(I-1,k)=Jolicoeur(radiiVec(:,1),I,alpha,flag_mat,flag_sym,flag_intf)*1e-6; % bending stiffness
        % the unit of compliance matrix cij is mm^2/N, need to time 1e-6
        % to convert the unit of EI to N*m^2
    end
    
    figure(k);
    layernum=1:1:M;
    plot(layernum,EI(:,k),'-bo','markers',3,'LineWidth',1.5);
    hold on;

    if flag_mat==1
        if flag_intf==1
            EIinf(k)=InfinityNoFriction(a,b,alpha,flag_mat,flag_sym)*1e-6;
        elseif flag_intf==2
            EIinf(k)=InfinityNoSlip(a,b,alpha,flag_mat,flag_sym)*1e-6;
        end
        plot([1,M],[EIinf(k),EIinf(k)],'--r','LineWidth',1.5);
        hleg1=legend('Finite layers','infinite layers');
    end

    xlabel('Number of layers','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
    ylabel('Efftive bending stiffness ($\ensuremath{\, \mathrm{Nm^2}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
    grid on;
    hold off;
    title(sprintf('angle = %5.2f',alpha/pi*180));
    
end

if flag_mat==1
    figure(length+1);
    plot(angle*180/pi,EIinf,'-bo','markers',3,'LineWidth',1.5);
    hold on;
    plot(angle*180/pi,EI(1,:),'-ro','markers',3,'LineWidth',1.5);
    hleg1=legend('infinite layers','One layer');
    xlabel('Angle ($^{\circ}$)','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
    ylabel('Efftive bending stiffness ($\ensuremath{\, \mathrm{Nm^2}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
    grid on;
    hold off;
end

display('All done!');

