clear all
clc

syms fa1 fa2 fs M N T fc

for i = -22:22
    c = coeff(i);
end

function c = coeff(m)
    fa1 = 1500;
    fa2 = 2500;
    fs = 8000;
    M = 22;
    N = (M * 2) + 1;
    T = 1/fs;
    fc = (fa2 - fa1)/2;
    
    if (m ~= 0)
        c = (1/(m*pi))*(sin(2*pi*m*T*fa2)-sin(2*pi*m*T*fa1));
    else
        c = 2*T*(fa2 - fa1);
    end
    
    fprintf('At m = %d,', m);
    fprintf(' c = %d', c);
    disp(' ');
end