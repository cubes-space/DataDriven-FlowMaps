function dxdt = model(t,x,theta) 

x = x(1) ;
y = x(2) ;
z = x(3) ; 

sigma = theta(1);
rho   = theta(2);
beta  = theta(3); 


%Original dynamics
dxdt(1) = sigma*(y-x);
dxdt(2) = x*(rho-z) - y;
dxdt(3) = x*y - beta*z;

dxdt = dxdt';


end
