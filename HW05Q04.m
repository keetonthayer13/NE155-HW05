%HW04 Q06

%here we will take a perscribed matrix which is in fact very similar to the
%matrix found in Q04 and write a program to implement the Jacobi,
%Gauss-Seidel, and Successive Over-Relaxation (SOR) methods....

%General iterative knowledge: We can say (A+S)x = Sx + b,(equiv. to Ax = b)
%and then say Q = A+S and x = Q^-1*Sx + Q^-1*b which leads to x = Px + b~
%and then we can use the iterative method but know that our error must
%converge

%NOTE: THIS FUNCTION WILL OUTPUT ALL THREE METHODS
function [x_J, iterations_J, x_GS, iterations_GS, x_SOR,...
                            iterations_SOR]= HW04Q06(n,x_0,w)
                        
%inputs -- n being the size of our matrix A, x_0 (NOTE: MUST BE A VECTOR!!)
%being our initial guess,and w (omega) which will be used in our final 
%iterative solution method as the prescaling factor

%outputs -- we will output the Jacobi, GS, and SOR results by x 
%approximation and number of iterations necessary in that order...
%i.e. first output is x_Jacobi, Jacobi iterations, x_GS, and so on

%first define our matrix A... similar to HW04Q04.. tri diagonal system with
%4 on the diagonal and -1 hugging it on both sides
format long
A = zeros(n);

A(1,1:2) = [4 -1]; A(end,end-1:end) = [-1 4];
b = 100*ones(n,1);

for j = 1:n-2
    A(j+1,j:j+2) = [-1 4 -1];
end

%Jacobi Iteration Solver: P = D^-1(D - A) = I - D^-1*A; D = diag(A)
    %we will simply define our needed matrices, with the appropriate
    %relationships provided above and then execute the iteration methods
D = diag(diag(A)); D_inv = inv(D); P = D_inv*(D-A); b_tilda = D_inv*b;

x_k = x_0;
x_knext = (P*x_0)+(b_tilda);
error = norm(abs(x_knext - x_k))/norm(abs(x_knext));;
count = 0;
while error > 10^-8
    x_knext = (P*x_k)+(b_tilda);
    error = norm(abs(x_knext - x_k))/norm(abs(x_knext));
    x_k = x_knext;
    count = count+1;
end
x_J = x_knext; iterations_J = count;
% if count < 28
% x_J = x_knext; iterations_J = count;
% else
%     x_J = 0; iterations_J = 0;
% end

%Gauss-Seidel Iteration Solver: A = L + U + D; (D+L)x_k = -Ux_k-1 + b;
%P = -(D+L)^-1*U
L = tril(A) - D; U = triu(A) - D;
P = -inv(D+L)*U; b_tilda = inv(D+L)*b;

x_k = x_0;
x_knext = (P*x_0)+(b_tilda);
error = norm(abs(x_knext - x_k))/norm(abs(x_knext));
count = 0;
while error > 10^-8
    x_knext = (P*x_k)+(b_tilda);
    error = norm(abs(x_knext - x_k))/norm(abs(x_knext));
    x_k = x_knext;
    count = count+1;
end
x_GS = x_knext; iterations_GS = count;
% if count < 16
% x_GS = x_knext; iterations_GS = count;
% else
%     x_GS = 0; iterations_GS = 0;
% end

%Successive Over Relaxation Iteration Solver: let's multiply through by
%omega (0<w<2) to try and speed up our method a bit...
%(D+wL)x = ((1-w)D-wU)x + wb  -> P = (D+wL)^-1*((1-w)D-wU)
P = inv(D+w*L)*((1-w)*D-(w*U)); b_tilda = inv(D+w*L)*w*b;

x_k = x_0;
x_knext = (P*x_0)+(b_tilda);
error = norm(abs(x_knext - x_k))/norm(abs(x_knext));
count = 0;
while error > 10^-8 %& count < 15
    x_knext = (P*x_k)+(b_tilda);
    error = norm(abs(x_knext - x_k))/norm(abs(x_knext));
    x_k = x_knext;
    count = count+1;
end

x_SOR = x_knext; iterations_SOR = count;
% w = w
% if count < 14
% x_SOR = x_knext; iterations_SOR = count;
% else
%     x_SOR = 0; iterations_SOR = 0;
% end
end