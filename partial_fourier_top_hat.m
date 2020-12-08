n=512; 
s=8; 
pf=7/8;

x=zeros(n,1); 
x(2*(n/8):(6*n/8)) = 1;

subplot (3,1,1);
plot (x); 
xlim ([0 n]); ylim ([-0.3 1.3])

y=fft(fftshift(x)); 
s = s*2;
y((n/s):(s-1)*n/s) = 0;
subplot (3,1,2);
x1 = ifftshift (ifft (y));
plot ([ real(x1) imag(x1) abs(x1) ]); 
xlim ([0 n]); ylim ([-0.3 1.3])
legend ({ 'real', 'imag', 'abs' })

2*(pf-0.5)*n/s
y((2*(pf-0.5)*n/s):(n/2)) = 0;
subplot (3,1,3);
x2 = ifftshift (ifft (y));
%plot ([ real(y) imag(y) abs(y) ]); 
plot ([ real(x2) imag(x2) abs(x2) ]); 
xlim ([0 n]); ylim ([-0.3 1.3])
legend ({ 'real', 'imag', 'abs' })

%
  %pf=0.5+0.1*(6-k); 
  %x=zeros(n,1); 
  %x((1+(n/2)-n*(pf-0.5)/s):((n/2)+(n/2)/s+1))=1; 
  %y=ifftshift(fft(fftshift(x))); 
%
  %subplot(3,2,k); 
  %plot ([ real(y) imag(y) abs(y)]); 
  %ylim ([-15 35]); 
  %grid on; 
  %title([ 'partial fourier: ' num2str(pf) ]);
%end
%
%legend ({ 'real', 'imag', 'abs' })


