function [ X_orig ] = generate_signal_n1_K3( Trepeat )

x_one_period = [1*ones(1,100) 2*ones(1,80) 1*ones(1,50) 3*ones(1,70) 2*ones(1,90) 1*ones(1,50) 3*ones(1,60)];
X_orig = kron(ones(1,Trepeat),x_one_period);

end

