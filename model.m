l1 = 0.187;
l2 = 0.34;
m1 = 0.25;
m2 = 0.15;
b1 = 0.038;
r2 = 0.0065;
g  = 9.8;

kb = 0.00767;
km = kb;
ra = 2.6;
la = 0.18e-3;
N  = 70;
je = 2.8e-6 + 2.27e-5 + 5e-7;
jm = 3.86e-7;

f = 0.005;

x0s = [[0;0.01;0;0] [0;pi/2-0.01;0;0] [0;pi/2-0.1;0;0]];
mos = [0;0;1];

x0 = x0s(:,2);
mo = 1;
va = 0;
[T,X] =ode45(@(t,x) deriv(t,mo,x,g,jm,kb,km,ra,la,je,l1,l2,b1,r2,m1,m2,...
    f,va,N),[0:0.02:20], x0);
figure(1);
plot(T,X(:,1)*180/pi,data02(:,1),data02(:,2),T,X(:,2)*180/pi,data02(:,1),...
    data02(:,3),"LineWidth",1);
xlabel("tempo (s)");
ylabel("graus")
legend(["theta1 sim","theta1 exp","theta2 sim","theta2 exp"],"fontSize",...
    10)
grid();

function dx= deriv(t,mo,x,g,jm,kb,km,ra,la,je,l1,l2,b1,r2,m1,m2,f,va,N)

the1 = x(1);
the2 = x(2);
dthe1= x(3);
dthe2= x(4);

ml11 = (m1/12)*((l1^2) + b1^2) + m2*(l1^2) + (jm*(N^2) + je)...
    + (m2/2)*((3*(r2^2) + l2^2)/6 - ((l2^2)/2)*((cos(the2)^2) -1))...
    + (m1*l1^2)/4;
ml12 = -(m2*l1*l2/2)*cos(the2);
ml21 = ml12;
ml22 = (m2/12)*6*(r2^2) + (m2*(l2^2))/4;

cl11 = mo*(N^2)*km*kb/ra;
cl12 = ((m2*l2^2)/2)*sin(the2)*cos(the2)*dthe1 + m2*l1*l2*sin(the2)...
    *dthe2/2;
cl21 = -(m2*l2*sin(the2)/4)*(l2*cos(the2)*dthe1 + l1*dthe2);
cl22 = (m2*l1*l2/4)*sin(the2)*dthe1;

cl=[cl11 cl12; cl21 cl22];

gv = [0;-m2*l2*g*sin(the2)/2];

d = [mo*N*km/ra;0];

alfa = ml11*ml22-ml21*ml12;

mlinv = (1/alfa)*[ml22 -ml12; -ml21 ml11];

x1 = [the1;the2];
x2 = [dthe1;dthe2];

dx = [zeros(2,2) eye(2); zeros(2,2) -mlinv*(cl+f*eye(2))]*[x1;x2] + ...
    [zeros(2,2) zeros(2,2); zeros(2,2) -mlinv*f*eye(2)]*[x1;x2] + ...
    [0;0;mlinv*(d*va - gv)]; 

end