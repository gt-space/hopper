close all; clc;
clear CL; clear X; clear lyapEigs; clear poles;
clear CL_noquat;clear poles_noquat;

Q = ones(13,13);
Q_noquat = ones(9,9);

for i = 1:length(t)-1

  CL(:,:,i) = A(:,:,i) - B(:,:,i)*K1(:,:,i);


  X(:,:,i) = lyap(CL(:,:,i)', Q);
  lyapEigs(:,i) = eig(X(:,:,i));
  poles(:,i)    = eig(CL(:,:,i));
  [V,D]    = eig(CL(:,:,i));

  if any(real(poles(:,i)) >= 0)
      check(i) = 0;
  else
      check(i) = 1;
  end

end

NOTSTABLE = find(check == 0);

%%
Aeigs = eig(A(:,:,1));
%Beigs = eig(B(:,:,1));
C = zeros(13,13); D = zeros(13,4)
sys = ss(A(:,:,1),B(:,:,1),C,D)



%%

for i = 1:13
    disp(['Pole ', num2str(i), ' = ', num2str(poles(i))])
    contribution(:,i) = (abs(V(:,i))/max(abs(V(:,i))));   % normalized participation by state
    disp(contribution(:,i))
end

%%




time = 1:size(poles,2);

figure
plot3(real(poles(:)), imag(poles(:)), repelem(time, size(poles,1))', 'rx')
xlabel('Real')
ylabel('Imag')
zlabel('Time index')
grid on

figure
plot(real(poles(:,1)), imag(poles(:,1)), 'rx', 'MarkerSize', 10)
hold on
xline(0)
yline(0)
grid on
xlabel('Real')
ylabel('Imag')

