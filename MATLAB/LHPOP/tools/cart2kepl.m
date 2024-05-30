function [elts] = cart2kepl(X_INER,t,varargin)
    if nargin>2
        GM = varargin{1};
    else
        GM = 4.9027926024e+03;
    end
    SMAtol = 1e-6; 
    ECCtol = 1e-6;
    ANGtol = 1e-6; 

    if size(X_INER,1)~=6
        X_INER = X_INER';
        if size(X_INER,1)~=6
            error('size not correct for input X_INER');
        end
    end
    
    elts = cspice_oscelt(X_INER,t,GM);
    ECC = elts(2,:);
    ECC(ECC<ECCtol) = 0;
    SMA = elts(1,:)./(1-ECC);
    INC = elts(3,:).*(180/pi);
    LAN = elts(4,:).*(180/pi);
    AOP = elts(5,:).*(180/pi);
    MA  = elts(6,:).*(180/pi);
    AOP(AOP<0) = AOP(AOP<0)+360;
    MA(MA<0)  = MA(MA<0)+360;
    SMA = round(SMA./SMAtol).*SMAtol;
    ECC = round(ECC./ECCtol).*ECCtol;
    INC = round(INC./ANGtol).*ANGtol;
    LAN = round(LAN./ANGtol).*ANGtol;
    AOP = round(AOP./ANGtol).*ANGtol;
    MA  = round(MA./ANGtol).*ANGtol;
    elts = [SMA;ECC;INC;LAN;AOP;MA];
      
end
