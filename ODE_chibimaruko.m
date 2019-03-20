function dydt = ODE_chibimaruko( t,y,om,ze,magForce)
%　振り子の振動関数　ODE_chibimaruko_pendulum.m
%       y(1)=θ
%       y(2)=dθ/dt
%　入力引数
%       om : 振子の固有角振動数
%       ze : 減衰比
%       magForce: 電磁力パルス強度
dydt=[y(2);
    -om^2*y(1)-2*ze*om*y(2)+magForce];