clear;

tmin=0;
tmax=30;

syms t
xt=sin(t);
yt=t^2/10;
zt=t;
fplot3(xt,yt,zt, [tmin tmax/2])

hold on
vxt=diff(xt)
vyt=diff(yt)
vzt=diff(zt)
fplot3(vxt,vyt,vzt, [tmin tmax])

axt=diff(vxt)
ayt=diff(vyt)
azt=diff(vzt)
fplot3(axt,ayt,azt, [tmin tmax])
hold off
