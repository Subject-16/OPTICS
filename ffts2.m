%creating function for 2d scaled fft
%f is the argument ,of which we want to calculate scaled fft
%nx and ny  is the dimension of f 
function f_scaled=ffts2(f,nx,ny)
%defining dimension of f
mx=-nx/2:1:nx/2;
my=-ny/2:1:ny/2;
[Mx My]=meshgrid(mx,my); 
%scaling factor
s=1;
xs=1i*pi*s/nx;
ks=exp(-xs*(Mx.^2+My.^2));
fs=fft2(f).*ifft2(f);
ns=-nx:1:nx;
[Mxs Mys]=meshgrid(ns,ns);
gs=ks.*(fs);
%padding symmetrically with zero
bs=zeros(2*ny+1,2*nx+1);
bs(ny/2:3*ny/2,nx/2:3*ny/2)=gs;
%creating function h
hs=exp(xs*(Mxs.^2+Mys.^2));
%multiplying element by element and doing ifft
Fs1=fft2(ifft2(bs).*ifft2(hs));
Fs=fftshift(Fs1);
%multiplying by factor
ks1=exp(-xs*(Mx.^2+My.^2));
bs1=zeros(2*ny+1,2*nx+1);
bs1(ny/2:3*ny/2,nx/2:3*nx/2)=ks1;
Fss=Fs.*bs1;
f_scaled=Fss(ny/2:3*ny/2,nx/2:3*ny/2);
end
 

