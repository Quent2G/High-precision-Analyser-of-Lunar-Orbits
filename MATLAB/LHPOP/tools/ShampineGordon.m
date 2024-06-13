%----------------------------------------------------------------------------
%
% Purpose:
%   Numerical integration methods for ordinaray differential equations
%
%   This module provides implemenation of the variable order variable 
%   stepsize multistep method of Shampine & Gordon.
% 
% Last modified:   2015/08/25   M. Mahooti
% 
% Reference:
%   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
%   Freeman and Comp., San Francisco (1975).
%
%----------------------------------------------------------------------------
function y = ShampineGordon(func,t,tout,relerr,abserr,n_eqn,y,model)

twou  = 2*eps;
fouru = 4*eps;

PermitTOUT = 1;

% Powers of two (two(n)=2^n)
two  = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,...
        256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0];

gstr = [1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,...
        0.0143, 0.0114, 0.00936, 0.00789, 0.00679,...
        0.00592, 0.00524, 0.00468];

% Allocate vectors with proper dimension
wt    = zeros(n_eqn,1);
p     = zeros(n_eqn,1);
yp    = zeros(n_eqn,1);
phi   = zeros(n_eqn,17);
g     = zeros(14,1);
sig   = zeros(14,1);
rho   = zeros(14,1);
w     = zeros(13,1);
alpha = zeros(13,1);
beta  = zeros(13,1);
v     = zeros(13,1);
psi_  = zeros(13,1);

% Return, if output time equals input time
if (t==tout)    % No integration
    return;
end

epsilon = max(relerr,abserr);
del    = tout - t;
absdel = abs(del);
tend   = t + 100.0*del;

if (~PermitTOUT)
    tend = tout;
end

kle4   = 0;
releps = epsilon;
abseps = epsilon;
start  = 1;
x      = t;
yy     = y;
h      = sign_( max(fouru*abs(x), abs(tout-x)), tout-x );

while (1)   % Start step loop

  % If already past output point, interpolate solution and return
  if (abs(x-t) >= absdel)
      yout  = zeros(n_eqn,1);
      ypout = zeros(n_eqn,1);
      g(2)   = 1.0;
      rho(2) = 1.0;
      hi = tout - x;
      ki = kold + 1;
      
      % Initialize w(*) for computing g(*)
      for i=1:ki
          temp1 = i;
          w(i+1) = 1.0/temp1;
      end
      % Compute g(*)
      term = 0.0;
      for j=2:ki
          psijm1 = psi_(j);
          gamma = (hi + term)/psijm1;
          eta = hi/psijm1;
          for i=1:ki+1-j
              w(i+1) = gamma*w(i+1) - eta*w(i+2);
          end
          g(j+1) = w(2);
          rho(j+1) = gamma*rho(j);
          term = psijm1;
      end
      
      % Interpolate for the solution yout and for
      % the derivative of the solution ypout      
      for j=1:ki
          i = ki+1-j;
          yout  = yout  + g(i+1)*phi(:,i+1);
          ypout = ypout + rho(i+1)*phi(:,i+1);
      end
      yout = y + hi*yout;
      y    = yout;
      t         = tout;             % Set independent variable
      told      = t;                % Store independent variable
      OldPermit = PermitTOUT;
      return                        % Normal exit
  end                         
  
  % If cannot go past output point and sufficiently close,
  % extrapolate and return
  if ( ~PermitTOUT && ( abs(tout-x) < fouru*abs(x) ) )
      h = tout - x;
      yp = func(x,yy,model);          % Compute derivative yp(x)
      y = yy + h*yp;            % Extrapolate vector from x to tout
      t         = tout;         % Set independent variable
      told      = t;            % Store independent variable
      OldPermit = PermitTOUT;
      return                    % Normal exit
  end
  
  % Limit step size, set weight vector and take a step
  h  = sign_(min(abs(h), abs(tend-x)), h);
  for l=1:n_eqn
      wt(l) = releps*abs(yy(l)) + abseps;
  end
  
%   Step
%                                                                   
% Begin block 0                                                     
%                                                                   
% Check if step size or error tolerance is too small for machine    
% precision.  If first step, initialize phi array and estimate a    
% starting step size. If step size is too small, determine an       
% acceptable one.                                                   
%                                                                   

if (abs(h) < fouru*abs(x))
    h = sign_(fouru*abs(x),h);
    crash = 1;
    return           % Exit 
end

p5eps  = 0.5*epsilon;
crash  = 0;
g(2)   = 1.0;
g(3)   = 0.5;
sig(2) = 1.0;

ifail = 0;

% If error tolerance is too small, increase it to an 
% acceptable value.                                  
round = 0.0;
for l=1:n_eqn
    round = round + (y(l)*y(l))/(wt(l)*wt(l));
end
round = twou*sqrt(round);
if (p5eps<round)
    epsilon = 2.0*round*(1.0+fouru);
    crash = 1;
    return
end

if (start)
  % Initialize. Compute appropriate step size for first step. 
  yp = func(x,y,model);
  sum = 0.0;
  for l=1:n_eqn
      phi(l,2) = yp(l);
      phi(l,3) = 0.0;
      sum = sum + (yp(l)*yp(l))/(wt(l)*wt(l));
  end
  sum  = sqrt(sum);
  absh = abs(h);
  if (epsilon<16.0*sum*h*h)
      absh=0.25*sqrt(epsilon/sum);
  end
  h    = sign_(max(absh, fouru*abs(x)), h);
  hold = 0.0;
  hnew = 0.0;
  k    = 1;
  kold = 0;
  start  = 0;
  phase1 = 1;
  nornd  = 1;
  if (p5eps<=100.0*round)
      nornd = 0;
      for l=1:n_eqn
          phi(l,16)=0.0;
      end
  end
end
%                                                                   
% End block 0                                                       
%                                                                   

%                                                                   
% Repeat blocks 1, 2 (and 3) until step is successful               
%                                                                   
while(1)
  
  %                                                                 
  % Begin block 1                                                   
  %                                                                 
  % Compute coefficients of formulas for this step. Avoid computing 
  % those quantities not changed when step size is not changed.     
  %                                                                 
  kp1 = k+1;
  kp2 = k+2;
  km1 = k-1;
  km2 = k-2;
  
  % ns is the number of steps taken with size h, including the 
  % current one. When k<ns, no coefficients change.           
  
  if (h ~=hold)
      ns=0;
  end
  if (ns<=kold)
      ns=ns+1;
  end
  nsp1 = ns+1;
  
  if (k>=ns)
      % Compute those components of alpha(*),beta(*),psi(*),sig(*) 
      % which are changed                                          
      beta(ns+1) = 1.0;
      realns = ns;
      alpha(ns+1) = 1.0/realns;
      temp1 = h*realns;
      sig(nsp1+1) = 1.0;
      if (k>=nsp1)
          for i=nsp1:k
              im1   = i-1;
              temp2 = psi_(im1+1);
              psi_(im1+1) = temp1;
              beta(i+1)  = beta(im1+1)*psi_(im1+1)/temp2;
              temp1    = temp2 + h;
              alpha(i+1) = h/temp1;
              reali = i;
              sig(i+2) = reali*alpha(i+1)*sig(i+1);
          end
      end
      psi_(k+1) = temp1;
      
      % Compute coefficients g(*); initialize v(*) and set w(*).
      if (ns>1)
          % If order was raised, update diagonal part of v(*)
          if (k>kold)
              temp4 = k*kp1;
              v(k+1) = 1.0/temp4;
              nsm2 = ns-2;
              for j=1:nsm2
                  i = k-j;
                  v(i+1) = v(i+1) - alpha(j+2)*v(i+2);
              end
          end
          
          % Update V(*) and set W(*)
          limit1 = kp1 - ns;
          temp5  = alpha(ns+1);
          for iq=1:limit1
              v(iq+1) = v(iq+1) - temp5*v(iq+2);
              w(iq+1) = v(iq+1);
          end
          g(nsp1+1) = w(2);
      else
          for iq=1:k
              temp3 = iq*(iq+1);
              v(iq+1) = 1.0/temp3;
              w(iq+1) = v(iq+1);
          end
      end
      
      % Compute the g(*) in the work vector w(*)
      nsp2 = ns + 2;
      if (kp1>=nsp2)
          for i=nsp2:kp1
              limit2 = kp2 - i;
              temp6  = alpha(i);
              for iq=1:limit2
                  w(iq+1) = w(iq+1) - temp6*w(iq+2);
              end
              g(i+1) = w(2);
          end
      end
  end % if K>=NS
  %
  % End block 1
  %
  
  %
  % Begin block 2
  %
  % Predict a solution p(*), evaluate derivatives using predicted
  % solution, estimate local error at order k and errors at orders
  % k, k-1, k-2 as if constant step size were used.
  %   
  
  % Change phi to phi star
  if (k>=nsp1)
      for i=nsp1:k
          temp1 = beta(i+1);
          for l=1:n_eqn
              phi(l,i+1) = temp1 * phi(l,i+1);
          end
      end
  end
  
  % Predict solution and differences 
  for l=1:n_eqn
      phi(l,kp2+1) = phi(l,kp1+1);
      phi(l,kp1+1) = 0.0;
      p(l)       = 0.0;
  end
  for j=1:k
      i     = kp1 - j;
      ip1   = i+1;
      temp2 = g(i+1);
      for l=1:n_eqn
          p(l)     = p(l) + temp2*phi(l,i+1);
          phi(l,i+1) = phi(l,i+1) + phi(l,ip1+1);
      end
  end
  if (nornd)
      p = y + h*p;
  else
      for l=1:n_eqn
          tau = h*p(l) - phi(l,16);
          p(l) = y(l) + tau;
          phi(l,17) = (p(l) - y(l)) - tau;
      end
  end
  xold = x;
  x = x + h;
  absh = abs(h);
  yp = func(x,p,model);
  
  % Estimate errors at orders k, k-1, k-2 
  erkm2 = 0.0;
  erkm1 = 0.0;
  erk = 0.0;
  
  for l=1:n_eqn
      temp3 = 1.0/wt(l);
      temp4 = yp(l) - phi(l,1+1);
      if (km2> 0)
          erkm2 = erkm2 + ((phi(l,km1+1)+temp4)*temp3)...
                         *((phi(l,km1+1)+temp4)*temp3);
      end
      if (km2>=0)
          erkm1 = erkm1 + ((phi(l,k+1)+temp4)*temp3)...
                         *((phi(l,k+1)+temp4)*temp3);
      end
      erk = erk + (temp4*temp3)*(temp4*temp3);
  end
  
  if (km2> 0)
      erkm2 = absh*sig(km1+1)*gstr(km2+1)*sqrt(erkm2);
  end
  if (km2>=0)
      erkm1 = absh*sig(k+1)*gstr(km1+1)*sqrt(erkm1);
  end
  
  temp5 = absh*sqrt(erk);
  err = temp5*(g(k+1)-g(kp1+1));
  erk = temp5*sig(kp1+1)*gstr(k+1);
  knew = k;
  
  % Test if order should be lowered 
  if (km2 >0)
      if (max(erkm1,erkm2)<=erk)
          knew=km1;
      end
  end
  if (km2==0)
      if (erkm1<=0.5*erk)
          knew=km1;
      end
  end
  %
  % End block 2
  %
  
  %
  % If step is successful continue with block 4, otherwise repeat
  % blocks 1 and 2 after executing block 3
  %
  success = (err<=epsilon);
  
  if (~success)
  
    %
    % Begin block 3
    %
    % The step is unsuccessful. Restore x, phi(*,*), psi(*). If
    % 3rd consecutive failure, set order to 1. If step fails more
    % than 3 times, consider an optimal step size. Double error
    % tolerance and return if estimated step size is too small
    % for machine precision.
    %
    
    % Restore x, phi(*,*) and psi(*)
    phase1 = 0; 
    x = xold;
    for i=1:k
        temp1 = 1.0/beta(i+1);
        ip1 = i+1;
        for l=1:n_eqn
            phi(l,i+1)=temp1*(phi(l,i+1)-phi(l,ip1+1));
        end
    end
    
    if (k>=2)
        for i=2:k
            psi_(i) = psi_(i+1) - h;
        end
    end
    
    % On third failure, set order to one. 
    % Thereafter, use optimal step size   
    ifail = ifail+1;
    temp2 = 0.5;
    if (ifail>3) 
      if (p5eps < 0.25*erk)
          temp2 = sqrt(p5eps/erk);
      end
    end
    if (ifail>=3)
        knew = 1;
    end
    h = temp2*h;
    k = knew;
    if (abs(h)<fouru*abs(x))
        crash = 1;
        h = sign_(fouru*abs(x), h);
        epsilon = epsilon*2.0;
        return                 % Exit 
    end    
    %
    % End block 3, return to start of block 1
    %
    
  end  % end if(success)
  
  if (success)
      break
  end
  
end

%
% Begin block 4
%
% The step is successful. Correct the predicted solution, evaluate
% the derivatives using the corrected solution and update the
% differences. Determine best order and step size for next step.
%

kold = k;
hold = h;

% Correct and evaluate
temp1 = h*g(kp1+1);
if (nornd)
    for l=1:n_eqn
        y(l) = p(l) + temp1*(yp(l) - phi(l,2));
    end
else
    for l=1:n_eqn
        rho = temp1*(yp(l) - phi(l,2)) - phi(l,17);
        y(l) = p(l) + rho;
        phi(l,16) = (y(l) - p(l)) - rho;
    end
end
yp = func(x,y,model);

% Update differences for next step 
for l=1:n_eqn
    phi(l,kp1+1) = yp(l) - phi(l,2);
    phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
end
for i=1:k
    for l=1:n_eqn
        phi(l,i+1) = phi(l,i+1) + phi(l,kp1+1);
    end
end

% Estimate error at order k+1 unless               
% - in first phase when always raise order,        
% - already decided to lower order,                
% - step size not constant so estimate unreliable  
erkp1 = 0.0;
if ( (knew==km1) || (k==12) )
    phase1 = 0;
end

if (phase1)
    k = kp1;
    erk = erkp1;
else
    if (knew==km1)
        % lower order 
        k = km1;
        erk = erkm1;
    else
        if (kp1<=ns)
            for l=1:n_eqn
                erkp1 = erkp1 + (phi(l,kp2+1)/wt(l))*(phi(l,kp2+1)/wt(l));
            end
            erkp1 = absh*gstr(kp1+1)*sqrt(erkp1);
            % Using estimated error at order k+1, determine 
            % appropriate order for next step               
            if (k>1)
                if ( erkm1<=min(erk,erkp1))
                    % lower order
                    k=km1; erk=erkm1;
                else
                    if ( (erkp1<erk) && (k~=12) )
                        % raise order 
                        k=kp1;
                        erk=erkp1;
                    end
                end
            elseif (erkp1<0.5*erk)
                % raise order
                % Here erkp1 < erk < max(erkm1,ermk2) else    
                % order would have been lowered in block 2.   
                % Thus order is to be raised                  
                k = kp1;
                erk = erkp1;
            end
        end % end if kp1<=ns
    end % end if knew!=km1
end % end if !phase1

% With new order determine appropriate step size for next step
if ( phase1 || (p5eps>=erk*two(k+2)) )
    hnew = 2.0*h;
else
    if (p5eps<erk)
        temp2 = k+1;
        r = p5eps/erk^(1.0/temp2);
        hnew = absh*max(0.5, min(0.9,r));
        hnew = sign_(max(hnew, fouru*abs(x)), h);
    else
        hnew = h;
    end
end
h = hnew;
%
% End block 4
%

  % Count number of consecutive steps taken with the order of
  % the method being less or equal to four and test for stiffness
  kle4 = kle4+1;
  if (kold>  4)
      kle4 = 0;
  end
  if (kle4>=50)
      stiff = 1;
  end
  
end % End step loop

