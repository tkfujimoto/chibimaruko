function dydt = ODE_chibimaruko( t,y,om,ze,magForce)
%�@�U��q�̐U���֐��@ODE_chibimaruko_pendulum.m
%       y(1)=��
%       y(2)=d��/dt
%�@���͈���
%       om : �U�q�̌ŗL�p�U����
%       ze : ������
%       magForce: �d���̓p���X���x
dydt=[y(2);
    -om^2*y(1)-2*ze*om*y(2)+magForce];