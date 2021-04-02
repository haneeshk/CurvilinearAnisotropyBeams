function PlotResult(rVec,sigmaz,sigmar,taurz,taurtheta,ur,utheta,w)

plot(rVec*1e3,sigmaz*1e-6,'--ks','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{m m}})','Interpreter','Latex','FontSize', 18,'rot',00,'FontName','Helvetica');
ylabel('Axial Stress $\sigma_{zz}\ensuremath{\, \mathrm{MPa}}$ ','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
grid on

figure
plot(rVec*1e3,sigmar*1e-6,'-.k+','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{ mm}})','Interpreter','Latex','FontSize', 18,'rot',0,'FontName','Helvetica');
ylabel('Radial Stress $\sigma_{rr} (\ensuremath{\, \mathrm{MPa}})$','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
grid on

figure
plot(rVec*1e3,taurz*1e-6,'--ks','markers',2)
xlabel('Radial Position (mm)','Interpreter','Latex','FontSize', 18,'rot',00,'FontName','Helvetica');
ylabel('Shear Stress $\tau_{rz} \ensuremath{\, \mathrm{MPa}}$','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
grid on

figure
plot(rVec*1e3,taurtheta*1e-6,'--ks','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{m m}})','Interpreter','Latex','FontSize', 18,'rot',0,'FontName','Helvetica');
ylabel('Shear Stress $\tau_{r\theta}\ensuremath{\, \mathrm{MPa}}$ ','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
grid on

figure
plot(rVec*1e3,utheta*1e3,'--ks','markers',2)
xlabel('Radial Position (\ensuremath{\, \mathrm{mm}})','Interpreter','Latex','FontSize', 18,'rot',0,'FontName','Helvetica');
ylabel('Tangential Displacement, $u_{\theta}$(mm)','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
%ylabel(' $\frac{\mathcal{L}^n[\mathbf{r}^n]-\mathcal{L}^0}{\mathcal{L}^0} \times {\scriptstyle{100}}$','Interpreter','Latex','FontSize', 24,'rot',0)
grid on


figure
plot(rVec*1e3,w*1e3,'--ks','markers',2)
xlabel('Radial Position ($\ensuremath{\, \mathrm{m m}}$)','Interpreter','Latex','FontSize', 18,'rot',0,'FontName','Helvetica');
ylabel('Axial Displacement W(\ensuremath{\, \mathrm{mm}})','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
grid on

figure
plot(rVec*1e3,ur*1e3,'--ks','markers',2)
xlabel('Radial Position ($\ensuremath{\, \mathrm{m m}}$)','Interpreter','Latex','FontSize', 18,'rot',0,'FontName','Helvetica');
ylabel('Radial Displacement U($\ensuremath{\, \mathrm{mm}}$)','Interpreter','Latex','FontSize', 18,'rot',90,'FontName','Helvetica');
grid on

end