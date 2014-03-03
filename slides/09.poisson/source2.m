function e = source2(x)
  e = zeros(length(x)*length(x),1);
  for j=1:length(x)
    for i=1:length(x)
      e(i+(j-1)*length(x)) = 8*pi*pi*sin(2*pi*x(i))*sin(2*pi*x(j));
    end
  end
