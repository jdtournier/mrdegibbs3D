n=512; 
s=16; 

for k=1:6
  pf=0.5+0.1*(6-k); 
  x=zeros(n,1); 
  x((1+(n/2)-n*(pf-0.5)/s):((n/2)+(n/2)/s+1))=1; 
  y=ifftshift(fft(fftshift(x))); 
  
  subplot(3,2,k); 
  plot ([ real(y) imag(y) abs(y)]); 
  ylim ([-15 35]); 
  grid on; 
  title([ 'partial fourier: ' num2str(pf) ]);
end

legend ({ 'real', 'imag', 'abs' })

