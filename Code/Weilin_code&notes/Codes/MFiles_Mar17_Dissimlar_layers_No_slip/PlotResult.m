function PlotResult(rVec,sigmar,sigmatheta,sigmaz,taurz,taurtheta,tauthetaz,ur,utheta,w)

figure(1)
plot(rVec*1e3,sigmar*1e-6,'-k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{m m}})','Interpreter','Latex','FontSize', 14,'rot',00,'FontName','Helvetica');
ylabel('Axial Stress $\sigma_{rr}\ensuremath{\, \mathrm{MPa}}$ ','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(2)
plot(rVec*1e3,sigmatheta*1e-6,'-k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{ mm}})','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Radial Stress $\sigma_{\theta\theta} (\ensuremath{\, \mathrm{MPa}})$','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(3)
plot(rVec*1e3,sigmaz*1e-6,'-k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{ mm}})','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Radial Stress $\sigma_{zz} (\ensuremath{\, \mathrm{MPa}})$','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(4)
plot(rVec*1e3,taurz*1e-6,'-k+','markers',2)
xlabel('Radial Position (mm)','Interpreter','Latex','FontSize', 14,'rot',00,'FontName','Helvetica');
ylabel('Shear Stress $\tau_{rz} \ensuremath{\, \mathrm{MPa}}$','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(5)
plot(rVec*1e3,taurtheta*1e-6,'-k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{m m}})','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Shear Stress $\tau_{r\theta}\ensuremath{\, \mathrm{MPa}}$ ','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(6)
plot(rVec*1e3,tauthetaz*1e-6,'-k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{m m}})','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Shear Stress $\tau_{z\theta}\ensuremath{\, \mathrm{MPa}}$ ','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(7)
plot(rVec*1e3,utheta*1e3,'-k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{mm}})','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Tangential Displacement $u_{\theta}$(mm)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(8)
plot(rVec*1e3,ur*1e3,'-k+','markers',2)
xlabel('Radial Position ($\ensuremath{\, \mathrm{m m}}$)','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Axial Displacement, $u_{r}$(mm)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

figure(9)
plot(rVec*1e3,w*1e3,'-k+','markers',2)
xlabel('Radial Position ($\ensuremath{\, \mathrm{m m}}$)','Interpreter','Latex','FontSize', 14,'rot',0,'FontName','Helvetica');
ylabel('Radial Displacement w($\ensuremath{\, \mathrm{mm}}$)','Interpreter','Latex','FontSize', 14,'rot',90,'FontName','Helvetica');
grid on

end
