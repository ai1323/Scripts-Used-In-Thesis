%%Scriot for Random shape generation (without input)


num_shapes = 20;%Number of Atoms 

% Range of Paremetrs 
a_range = [1, 1];
b_range = [1, 1];

m_range = [4,10];%Try larger values

n_range = [0.1, 100];%No more than 100
% Generate shapes
for i = 1:num_shapes
% Draw a random number for each parameters 
 a1 = rand_range(a_range); %Keep constant 
 b1 = rand_range(b_range); %Keep constant 
 
 m1 = rand_range(m_range);

 n11 = rand_range(n_range);
 n21 = rand_range(n_range);

 n31 = rand_range(n_range);
 a2 = rand_range(a_range);%Keep Constant 
 b2 = rand_range(b_range);%Keep Constant 
 
 m2 = rand_range(m_range);

 n12 = rand_range(n_range);
 n22 = rand_range(n_range);

 n32 = rand_range(n_range);
% Recall Function (if error massage here restart matlab)
 SuperFormula3D(a1, b1, m1, n11, n21, n31, a2, b2, m2, n12, n22, n32);
 close all;
 fprintf('Generated shape %d of %d\n', i, num_shapes);
end
%Get the random value  within the specified range 
function val = rand_range(range)
 val = range(1) + (range(2) - range(1)) * rand();
end
%