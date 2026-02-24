function Z = methane_compression_factor(T, P)

%Z = ones(size(T));
%greturn

T = max(T, 150);

gammag = 0.5654;
Pc = 756.8 - 131 .* gammag - 3.6 * gammag.^2;
Tc = 169.2 + 349.5 .* gammag - 74 * gammag.^2;

a1=0.317842;
a2=0.382216;
a3=-7.76835;
a4=14.2905;
a5=2e-6;
a6=-0.004693;
a7=0.096254;
a8=0.16672;
a9=0.96691;
a10=0.063069;
a11=-1.96685;
a12=21.0581;
a13=-27.0246;
a14=16.23;
a15=207.783;
a16=-488.161;
a17=176.29;
a18=1.88453;
a19=3.05921;

Tpr = (T * 1.8) ./ Tc;
Ppr = (P / 100000 * 14.5) ./ Pc;

t = 1 ./ Tpr;

A = a1.*exp(a2.*(1-t).^2).*Ppr.*t;
B = a3.*t+a4.*t.^2+a5.*Ppr.^6.*t.^6;
Cc = a9+a8.*Ppr.*t+a7.*Ppr.^2.*t.^2+a6.*Ppr.^3.*t.^3;
Dc = a10.*exp(a11.*(1-t).^2).*t;
Ec = a12.*t+a13.*t.^2+a14.*t.^3;
F = a15.*t+a16.*t.^2+a17.*t.^3;
G = a18+a19.*t;

y = (Dc.*Ppr)./(-((A.^2.*B)./Cc.^3)+(1+A.^2)./Cc);

Z = (Dc.*Ppr.*(1+y+y.^2-y.^3))./((1-y).^3.*(Dc.*Ppr+Ec.*y.^2-F.*y.^G));

if (any(imag(Z) ~= 0))
    u = 1;
end

end
