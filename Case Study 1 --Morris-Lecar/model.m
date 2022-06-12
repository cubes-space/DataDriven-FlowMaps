function dxdt = model(t,x,theta) 

V = x(1) ;
N = x(2) ;

Cm = theta(1);
gl = theta(2);
Vl = theta(3); 
gca = theta(4);
Vca = theta(5);
gk = theta(6) ;
Vk = theta(7); 
V1 = theta(8);
V2 = theta(9);
V3 = theta(10);
V4 = theta(11);
phi = theta(12);
Iapp = theta(13);

Minf = 0.5*(1+tanh((V-V1)/V2));
Ninf = 0.5*(1+tanh((V-V3)/V4));
lambdaN = phi*cosh((V-V3)/(2*V4));


%Original dynamics
dxdt(1) =  (-gl*(V-Vl) - gca*(V - Vca) * Minf - gk*(V - Vk)*N + Iapp)/Cm; 
dxdt(2) = lambdaN*(Ninf - N) ;

dxdt = dxdt';


end
