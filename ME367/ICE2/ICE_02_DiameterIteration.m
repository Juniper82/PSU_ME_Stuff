%% 10 Obtain the shaft diameter based on trial guess at various locations
% safety factor
Nf = 2.5;
fprintf(1,'Factor of Safety = %f\n',Nf);

disp(' ');
disp('The required shaft diameter at point C:')
FUNC = @(d) d-(((32*Nf)/pi)*...
    (((sqrt((Kfb*MCa)^2 + (3/4)*(Kfs*Ta)^2))/Se(d))+...
    ((sqrt((Kfbm*MCm)^2 + (3/4)*(Kfsm*Tm)^2))/Sut)))^(1/3);
dc = fzero(FUNC,0.5);
disp(['Dia @ C = ', num2str(dc),' inches']);
% disp(['Csize(d) =', num2str(Csize(dc))]);    % no units
% disp(['Se(d) =', num2str(Se(dc)), 'psi']);    % psi

disp(' ');
disp('The required shaft diameter at point B:');
Kfb = 1 + q * (Ktk - 1);
Kfs = 1 + qs * (Ktk - 1);

Kfsm = Kfs;
Kfbm = Kfb;
fprintf(1,'Fatigue Stress Conc. @B (keyway) Bending = %g, Torsion = %g\n',Kfb, Kfs);
FUNC = @(d) d-(((32*Nf)/pi)*...
    (((sqrt((Kfb*MBa)^2 + (3/4)*(Kfs*Ta)^2))/Se(d))+...
    ((sqrt((Kfbm*MBm)^2 + (3/4)*(Kfsm*Tm)^2))/Sut)))^(1/3);
db = fzero(FUNC,0.5);
disp(['Dia @ B = ', num2str(db),' inches']); 
% disp(['Csize(d) =', num2str(Csize(db))])    % no units
% disp(['Se(d) =', num2str(Se(db)), 'psi'])    % psi

disp(' ');
disp('The required shaft diameter at point D:');

Kfb = 1 + q * (Ktb - 1);
Kfs = 1 + qs * (Kts - 1);

Kfsm = Kfs;
Kfbm = Kfb;
fprintf(1,'Fatigue Stress Conc. @D Bending = %g, Torsion = %g\n',Kfb, Kfs);
FUNC = @(d) d-(((32*Nf)/pi)*...
    (((sqrt((Kfb*MDa)^2 + (3/4)*(Kfs*Ta)^2))/Se(d))+...
    ((sqrt((Kfbm*MDm)^2 + (3/4)*(Kfsm*Tm)^2))/Sut)))^(1/3);
dd = fzero(FUNC,0.5);
disp(['Dia @ D = ', num2str(dd),' inches']);
% disp(['Csize(d) =', num2str(Csize(dd))])    % no units
% disp(['Se(d) =', num2str(Se(dd)), 'psi'])    % psi

fprintf(1,'================================================================\n');
fprintf(1,'Dia @ C = %g, Dia @ B = %g, Dia @ D = %g\n',dc,db,dd);
fprintf(1, 'Se @C = %g psi, Se @B = %g psi, Se @D = %g psi\n',Se(dc),Se(db),Se(dd));
fprintf(1,'================================================================\n');

