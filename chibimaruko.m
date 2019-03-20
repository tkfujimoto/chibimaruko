function chibimaruko
%   Chibimaruko  Run a demo of a chibimaruko pendulum
%   with collision at an amplitude limitation.
%   This is an example of repeated event location, where the initial
%   conditions are changed after each terminal event.  This demo computes 
%   the vibration with ten collisions using calls to ODE23.  
%   The speed of the pendulum is attenuated by 0.9 after each collision.
%   The trajectory is plotted using the output function ODEPLOT.
%
%   See also ODE23, ODE45, ODESET, ODEPLOT, FUNCTION_HANDLE.

%   T. Fujimoto made this program based on the program "ballode" made by Mark W. Reichelt
%   and Lawrence F. Shampine, 1/3/95 (Copyright 1984-2014 The MathWorks, Inc.

clear;close all
freq_pendulum=2.03; %振子の固有振動数[Hz]
om=2*pi*freq_pendulum;%振子の固有角振動数[rad/s]
ze=0.0001;%減衰比
magForce=250; %電磁力パルスの強度[rad/s^2]
freq_pulse=0.6;%電磁パルスの振動数[Hz]
Tpulse=1/freq_pulse%電磁パルスの周期[s]
dutyratio=0.01;%電磁パルスのデューティ比
Tforce=Tpulse*dutyratio;%電磁パルスの作用時間[s]
ek=0.75;%反発係数
max_angle=15;
y0 = [0;0];%初期値[rad;rad/s]
refine = 4;%ode23の精度指定
options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,...
   'Refine',refine);

% Preparation of Accumulate output
tout = [0];
yout = y0.';
fout=[0.1];
teout = [];
yeout = [];
ieout = [];

nstep=1;
nmax=5

fig = figure;
ax = axes;
ax.XLim = [0 nmax*Tpulse];  %表示時間幅
ax.YLim = [-30/180*pi 30/180*pi];  %表示振幅幅
box on
hold on;

while(nstep<=nmax)
tstart=(nstep-1)*Tpulse;
tfinal=(nstep-1)*Tpulse+Tforce;%nstep番目の強制力が作用する時間スパン  

while(tstart<tfinal)
 % Solve until the first terminal event.
   [t,y,te,ye,ie] = ode23(@(t,y) ODE_chibimaruko_pendulum(t,y,om,ze,magForce),[tstart tfinal],y0,options);
   if ~ishold
      hold on
   end
   % Accumulation output:  This could be passed out as output arguments.
   nt = length(t);
   tout = [tout; t(2:nt)];
   yout = [yout; y(2:nt,:)];
   fout=[fout;10*ones(nt-1,1)];
   if (ie==1)
       teout = [teout; te];          % Events at tstart are never reported.
       yeout = [yeout; ye];
       ieout = [ieout; ie];
       % Set the new initial conditions, with .9 attenuation.
       y0(1) =y(nt,1);
       y0(2) = -ek*y(nt,2);
   else
       y0(1) =y(nt,1);
       y0(2) =y(nt,2);
   end
   ud = fig.UserData;
   if ud.stop
      break;
   end
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
      'MaxStep',t(nt)-t(1));
    tstart = t(nt);
end

tfinal=nstep*Tpulse;%nstep番目の自由振動をする時間スパン
while(tstart<tfinal)
 % Solve until the first terminal event.
   [t,y,te,ye,ie] = ode23(@(t,y) ODE_chibimaruko_pendulum(t,y,om,ze,0),[tstart tfinal],y0,options);
   if ~ishold
      hold on
   end
   % Accumulation output:  This could be passed out as output arguments.
   nt = length(t);
   tout = [tout; t(2:nt)];
   yout = [yout; y(2:nt,:)];
   fout=[fout;zeros(nt-1,1)];
   if(ie==1)
       teout = [teout; te];% Events at tstart are never reported.
       yeout = [yeout; ye];
       ieout = [ieout; ie];
       % Set the new initial conditions, with .9 attenuation.
       y0(1) =y(nt,1);
       y0(2) = -ek*y(nt,2);
   else
       y0(1) =y(nt,1);
       y0(2) =y(nt,2);
   end
   ud = fig.UserData;
   if ud.stop
      break;
   end
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
      'MaxStep',t(nt)-t(1));
    tstart = t(nt);
end
nstep=nstep+1;
end
hold off

figure(2)
ax = axes;
ax.XLim = [0 nmax*Tpulse];  %表示時間幅
%ax.YLim = [-25/180*pi 25/180*pi];  %表示振幅幅[rad]
ax.YLim = [-max_angle*2 max_angle*2];  %表示振幅幅[deg]

box on

hold on
plot(tout,yout(:,1)*180/pi)
if (~isempty(ieout))
plot(teout,yeout(:,1)*180/pi,'bo')
end
plot(tout,fout,'r')
xlabel('time [s]');
ylabel('Angle [deg]');
title('chibimaruko trajectory and the events');
hold off
%odeplot([],[],'done');

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = abs(y(1))-15/180*pi;  %detect the limit angle (角度振幅制限)
isterminal = 1;   % stop the integration
direction = 1;   % positive direction
