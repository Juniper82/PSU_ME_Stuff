%% ICE 1
% Shaft Design for Steady Torsion and Fully Reversed Bending
% Revision: Spring 2019
% Shaft diameter iteration - can use one of two techniques
%%
% 10 Obtain the shaft diameter based on trial guess
disp(' ');
Nf = 2.5
disp('The shaft diameter at point C:');

%begin approach 1
RHS = @(d) (((32*Nf)/pi)*((Kfb*(MC/Se(d))).^2 +...
    (3/4)*(Kfsm*(T/(Sy))).^2).^(1/2)).^(1/3);

d = 0.5     % in
cnt = 1;
while ((cnt < 10) && (abs(d-RHS(d))> 0.0001))
    dc(cnt) = d;
    rhsc(cnt) = RHS(d);
    d = rhsc(cnt);
    cnt=cnt+1;
end
fprintf(1,'dc = %g in after %d trial\n',d,cnt);
%end approach 1
%====================================================================
% Begin Alternate Approach
% An alternate way for doing the above would be to find zero crossing
% Uncomment this and comment above to try this approach
% RHS = @(d) d - (((32*Nf)/pi)*((Kfb*(MC/Se(d))).^2 +...
%     (3/4)*(Kfsm*(T/(Sy))).^2).^(1/2)).^(1/3);
% d = fzero(RHS,0.5);
% dc = d;
% rhsc = RHS(dc)+dc;
% End Alternate Approach
%====================================================================
fprintf(1,'========================================\n');
fprintf(1,'dc = %g,rhsc = %g \n',dc,rhsc);
fprintf(1,'========================================\n');
disp(['Csize(d) =', num2str(Csize(d))]);
disp(['Se(d) =', num2str(Se(d)),' psi']);
disp(' ')

% d2 = 0.544   % in - This was done after the fact
%
%% 11  Stress concentration for section B on shaft
disp('Fatigue stress concentration for section B :');
Kfb = 1 + q*(Ktk - 1)
Kfs = 1 + qs*(Ktk - 1)
Kfsm = Kfs

%
%% 12 Find the shaft diameter at B
% =======      begin Approach 1
RHS = @(d) (((32*Nf)/pi)*((Kfb*(MB/Se(d))).^2 +...
    (3/4)*(Kfsm*(T/(Sy))).^2).^(1/2)).^(1/3);

d = 0.5     % in
cnt = 1;
while ((cnt < 10) && (abs(d-RHS(d))> 0.0001))
    db(cnt) = d;
    rhsb(cnt) = RHS(d);
    d = rhsb(cnt);
    cnt=cnt+1;
end
% ==========   end approach 1
%====================================================================
% Begin Alternate Approach
% An alternate way for doing the above would be to find zero crossing
% Uncomment this and comment above to try this approach
% RHS = @(d) d - (((32*Nf)/pi)*((Kfb*(MB/Se(d))).^2 +...
%     (3/4)*(Kfsm*(T/(Sy))).^2).^(1/2)).^(1/3);
% d = fzero(RHS,0.5);
% db = d;
% rhsb = RHS(db)+db;
% End Alternate Approach
%====================================================================
fprintf(1,'========================================\n');
fprintf(1,'db = %g,rhsb = %g \n',dc,rhsb);
fprintf(1,'========================================\n');
disp(' ');
disp('The shaft diameter at point B:');
disp(['Csize(d) =', num2str(Csize(d))]);
disp(['Se(d) =', num2str(Se(d)),' psi']);
disp(' ')

% d1 = 0.625  % in -  normally done after iterations

%
%% 13 find shaft diameter at D
disp('Fatigue stress concentration for section D :');
Kfb = 1 + q*(Ktb - 1)
Kfs = 1 + qs*(Kts - 1)
Kfsm = Kfs
% ========== begin approach 1
RHS = @(d) (((32*Nf)/pi)*((Kfb*(MD/Se(d))).^2 +...
    (3/4)*(Kfsm*(T/(Sy))).^2).^(1/2)).^(1/3);

d = 0.5     % in
cnt = 1;
while ((cnt < 10) && (abs(d-RHS(d))> 0.0001))
    dd(cnt) = d;
    rhsd(cnt) = RHS(d);
    d = rhsd(cnt);
    cnt=cnt+1;
end
% ========== end approach 1
%====================================================================
% Begin Alternate Approach
% An alternate way for doing the above would be to find zero crossing
% Uncomment this and comment above to try this approach
% RHS = @(d) d - (((32*Nf)/pi)*((Kfb*(MD/Se(d))).^2 +...
%     (3/4)*(Kfsm*(T/(Sy))).^2).^(1/2)).^(1/3);
% d = fzero(RHS,0.5);
% dd = d;
% rhsd = RHS(dd)+dd;
% End Alternate Approach
%====================================================================

fprintf(1,'========================================\n');
fprintf(1,'dd = %g in,rhsd = %g in\n',dd,rhsd);
fprintf(1,'========================================\n');
disp(' ');
disp('The shaft diameter at point D:');
% disp(['Csize(d) =', num2str(Csize(d)), 'inch']);
disp(['Se(d) =', num2str(Se(d)),' psi']);
disp(' ')

% d3 = 0.500  % in normally done after the fact
fprintf(1,'================================================================\n');
fprintf(1,'Dia @ C = %g in, Dia @ B = %g in, Dia @ D = %g in\n',dc(end),db(end),dd(end));
fprintf(1, 'Se @C = %g psi, Se @B = %g psi, Se @D = %g psi\n',Se(dc(end)),Se(db(end)),Se(dd(end)));
fprintf(1,'================================================================\n');


