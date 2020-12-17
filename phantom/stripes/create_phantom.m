addpath ~/mrtrix3/matlab/

image = read_mrtrix ('phantom.mif')

x = 1;
w = 1;
xmax = size(image.data,1)

high = 100;
low = 50;

while x <= xmax
  x2 = x+w;
  if x2 > xmax
    x2 = xmax;
  end
  image.data(x:x2,:,:) = high;
  x = x+w;
  if x < xmax
    x2 = x+w;
    if x2 > xmax
      x2 = xmax;
    end
    image.data(x:(x+w),:,:) = low;
    x = x+w;
  end
  w =w+1;
end


write_mrtrix (image, 'phantom_stripes.mif')

