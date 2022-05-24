function [E,A,B,C,D] = ctdsx(nr,parin)
%CTDSX
%
% Usage:  [E,A,B,C,D] = ctdsx(nr,parin)
%         [E,A,B,C,D] = ctdsx(nr)
%
% Main routine of the benchmark library CTDSX (Version 1.0) described 
% in [1]. It generates benchmark examples for time-invariant, 
% continuous-time, dynamical systems
%         .
%       E x(t)  =  A x(t) + B u(t)
%                                                                  (1)
%         y(t)  =  C x(t) + D u(t)
%
% E, A are real n-by-n matrices, B is n-by-m, C is p-by-n, and 
% D is p-by-m. 
%
% Input:
%  - nr    : index of the desired example according to [1]; 
%            nr is a 1-by-2 matrix;
%            nr(1) defines the group:
%             = 1 : parameter-free problems of fixed size
%             = 2 : parameter-dependent problems of fixed size
%             = 3 : parameter-free problems of scalable size
%             = 4 : parameter-dependent problems of scalable size
%            nr(2) defines the number of the benchmark example within
%            a certain group.
%  - parin : parameters of the chosen example; 
%            referring to [1], the entries in parin have the following 
%            meaning:
%            Ex. 2.1 : parin(1)   = epsilon
%            Ex. 2.2 : parin(1)   = epsilon
%            Ex. 2.3 : parin(1)   = s
%            Ex. 2.4 : parin(1:7) = [b mu r r_c k_l sigma a]
%            Ex. 2.5 : parin(1)   = s
%            Ex. 2.6 : parin(1)   = s
%            Ex. 2.7 : parin(1:2) = [mu nu]
%            Ex. 3.1 : parin(1)   = q
%            Ex. 3.2 : parin(1)   = n
%            Ex. 3.3 : parin(1)   = n
%            Ex. 3.4 : parin(1)   = l
%            Ex. 4.1 : parin(1)   = n
%	               parin(2:4) = [a, b, c]
%                      parin(5:6) = [beta_1,beta_2]
%	               parin(7:8) = [gamma_1,gamma_2]
%            Ex. 4.2 : parin(1)   = l
%                      parin(2:4) = [mu delta kappa].
%            parin is optional; default values as defined in [1] are
%            used as example parameters if parin is omitted. Note that
%            parin is not referenced if nr(1) = 1.
%
% Output:
%  - E, A, B, C, D :  matrices of the dynamical system (1).
%
% References: 
%
% [1] D. Kressner, V. Mehrmann, and T. Penzl.
%     CTDSX - a Collection of Benchmark Examples for State-Space 
%     Realizations of Continuous-Time Dynamical Systems.
%     SLICOT working note 1998-9. 1998.
%
%     For questions concerning the collection or for the submission of
%     test examples, please contact Volker Mehrmann
%     (Email: volker.mehrmann@mathematik.tu-chemnitz.de).

%  D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz).
%  Dec 1, 1998. 

if nargin < 1,
  error('Not enough input arguments.');
end;

if length(nr) < 2,
  error('Please use the nr = [group, example] notation.');
end;


if nr(1) == 1,
  if nr(2) == 1,
    %  Example 1.1:  Laub 1979, Ex.1
    E = eye(2);
    A = [0 1; 0 0];
    B = [0; 1];
    C = eye(2);
    D = zeros(2,1);

  elseif nr(2) == 2,
    %  Example 1.2:  Laub 1979, Ex.2: uncontrollable-unobservable data
    E = eye(2);
    A = [4 3; -4.5 -3.5];
    B = [1; -1];
    C = eye(2);
    D = zeros(2,1);

  elseif nr(2) == 3,
    %  Example 1.3:  Beale/Shafai 1989: model of L-1011 aircraft
    E = eye(4);
    A = [0      1       0     0;
         0     -1.89    0.39 -5.53;
         0     -0.034  -2.98  2.43;
         0.034 -0.0011 -0.99 -0.21];
    B = [0 0; 0.36 -1.6; -0.95 -0.032; 0.03 0];
    C = eye(4);
    D = zeros(4,2);

  elseif nr(2) == 4,
    %  Example 1.4:  Bhattacharyya et al. 1983: binary distillation column
    E = eye(8);
    A = [-0.991   0.529   0       0       0       0       0       0     ;
          0.522  -1.051   0.596   0       0       0       0       0     ;
          0       0.522  -1.118   0.596   0       0       0       0     ;
          0       0       0.522  -1.548   0.718   0       0       0     ;
          0       0       0       0.922  -1.640   0.799   0       0     ;
          0       0       0       0       0.922  -1.721   0.901   0     ;
          0       0       0       0       0       0.922  -1.823   1.021 ;
          0       0       0       0       0       0       0.922  -1.943];
    B = 0.001*[ 3.84  4.00 37.60  3.08  2.36  2.88  3.08  3.00;
               -2.88 -3.04 -2.80 -2.32 -3.32 -3.82 -4.12 -3.96]';
    C = eye(8);
    D = zeros(8,2);

  elseif nr(2) == 5,
    %  Example 1.5:  Patnaik et al. 1980: tubular ammonia reactor
    E = eye(9);
    A = [-4.019  5.12   0      0     -2.082   0      0     0    0.87;
         -0.346  0.986  0      0     -2.34    0      0     0    0.97;
         -7.909 15.407 -4.069  0     -6.45    0      0     0    2.68;
        -21.816 35.606 -0.339 -3.87 -17.8     0      0     0    7.39;
        -60.196 98.188 -7.907  0.34 -53.008   0      0     0    20.4;
          0      0      0      0     94.0  -147.2    0    53.2   0;
          0      0      0      0      0      94.0 -147.2   0     0;
          0      0      0      0      0      12.8    0   -31.6   0;
          0      0      0      0     12.8     0      0    18.8 -31.6];
    B = [ 0.010  0.003  0.009  0.024  0.068   0 0 0 0;
         -0.011 -0.021 -0.059 -0.162 -0.445   0 0 0 0;
         -0.151  0      0      0      0       0 0 0 0]';
    C = eye(9);
    D = zeros(9,3);

  elseif nr(2) == 6,
    %  Example 1.6:  Davison/Gesing 1978: J-100 jet engine
    A1 = 0; A121 = 0; A122 = 0; A123 = 0; A21 = 0; A22 = 0; A23 = 0; A3 = 0;
    A31 = 0; B21 = 0; B22 = 0; B23 = 0; C = 0;
    load ctds106
    E = eye(30);
    A = [A1, A121, zeros(16,1), A122, zeros(16,2), A123, zeros(16,8);
         zeros(2,16), A21, zeros(2,12);
         zeros(3,18), A22, zeros(3,9);
         zeros(3,21), A23, zeros(3,6);
         A31, zeros(6,8), A3];
    B = [zeros(16,3);  B21, zeros(2,2);  zeros(3,1), B22, zeros(3,1);
         zeros(3,2), B23;  zeros(6,3)];
    D = zeros(5,3);

  elseif nr(2) == 7,
    %  Example 1.7:  Davison 1967: binary distillation column
    E = eye(11);
    A = [-1.40e-02 4.30e-03 0 0 0 0 0 0 0 0 0;
          9.50e-03 -1.38e-02 4.60e-03 0 0 0 0 0 0 0 5.00e-04;
          0 9.50e-03 -1.41e-02 6.30e-03 0 0 0 0 0 0 2.00e-04;         
          0 0 9.50e-03 -1.58e-02 1.10e-02 0 0 0 0 0 0;           
          0 0 0 9.50e-03 -3.12e-02 1.50e-02 2.20e-02 0 0 0 0;
          0 0 0 0 2.02e-02 -3.52e-02 2.20e-02 0 0 0 0;            
          0 0 0 0 0 2.02e-02 -4.22e-02 2.80e-02 0 0 0;
          0 0 0 0 0 0 2.02e-02 -4.82e-02 3.70e-02 0 2.00e-04;
          0 0 0 0 0 0 0 2.02e-02 -5.72e-02 4.20e-02 5.00e-04;
          0 0 0 0 0 0 0 0 2.02e-02 -4.83e-02 5.00e-04;
          2.55e-02 0 0 0 0 0 0 0 0 2.55e-02 -1.85e-02];
    B = [ 0         0        0;
          5.00e-06 -4.00e-05 2.50e-03;
          2.00e-06 -2.00e-05 5.00e-03;
          1.00e-06 -1.00e-05 5.00e-03;
          0         0        5.00e-03;
          0         0        5.00e-03;
         -5.00e-06  1.00e-05 5.00e-03;
         -1.00e-05  3.00e-05 5.00e-03;
         -4.00e-05  5.00e-06 2.50e-03;
         -2.00e-05  2.00e-06 2.50e-03;
          4.60e-04  4.60e-04 0];
    C = [ 0 0 0 0 0 0 0 0 0 1 0;
          1 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 1];
    D = zeros(3);

  elseif nr(2) == 8,
    %  Example 1.8:  Chien/Ergin/Ling/Lee 1958: drum boiler
    E = eye(9);
    A = [-3.93e+00 -3.15e-03  0         0         0         4.03e-05  0         0         0;
          3.68e+02 -3.05e+00  3.03e+00  0         0        -3.77e-03  0         0         0;
          2.74e+01  7.87e-02 -5.96e-02  0         0        -2.81e-04  0         0         0;
         -6.47e-02 -5.20e-05  0        -2.55e-01 -3.35e-06  3.60e-07  6.33e-05  1.94e-04  0;
          3.85e+03  1.73e+01 -1.28e+01 -1.26e+04 -2.91e+00 -1.05e-01  1.27e+01  4.31e+01  0;
          2.24e+04  1.80e+01  0        -3.56e+01 -1.04e-04 -4.14e-01  9.00e+01  5.69e+01  0;
          0         0         2.34e-03  0         0         2.22e-04 -2.03e-01  0         0;
          0         0         0        -1.27e+00  1.00e-03  7.86e-05  0        -7.17e-02  0;
         -2.20e+00 -1.77e-03  0        -8.44e+00 -1.11e-04  1.38e-05  1.49e-03  6.02e-03 -1.00e-10];
    B = [ 0         0        0;
          0         0        0;
          1.56e+00  0        0;
          0        -5.13e-06 0;
          8.28e+00 -1.50e+00 3.95e-02; 
          0         1.78e+00 0;
          2.33e+00  0        0;
          0        -2.45e-02 2.84e-03; 
          0         2.94e-05 0 ];
    C = [ 0 0 0 0 0 1 0 0 0;
          0 0 0 0 0 0 0 0 1 ];
    D = zeros(2,3);

  elseif nr(2) == 9,
    %  Example 1.9:  Ly, Gangsaas 1981: B-767 airplane
    A = 0; C = 0;
    load ctds109
    E = eye(55);
    B = zeros(55,2);
    B(48,1) = 8.0E+05;
    B(51,2) = 8.0E+05;
    D = zeros(2);

  elseif nr(2) == 10,
    %  Example 1.10:  control surface servo for an underwater vehicle
    E = eye(8);
    A = [ 0   85     0    0    0  0  0     0;
        -85 -120 -4100    0    0  0  0     0;
         33    0   -33    0 -700  0  0     0;
          0    0     0    0 1400  0  0     0;
          0    0  1600 -450 -110  0  0     0;
          0    0     0   81    0 -1  0  -900;
          0    0     0    0    0  0  0   110;
          0    0     0    0    0 12 -1.1 -22 ];
    B = zeros(8,2);
    B(2,:) = [4.6 9.9e+04];
    C = [0 0 0 0 0 0 1 0];
    D = [0 0];

  else
    str = sprintf('Example #%i is not available in Group #%i !',nr(2),nr(1));
    error(str);
  end;


elseif nr(1) == 2,
  if nr(2) == 1,
    %  Example 2.1:  Chow/Kokotovic 1976: magnetic tape control system
    if nargin < 2, d = 10^(-6);  else d = parin(1); end
    if d==0, error('Parameter must not be zero.'); end
    E = eye(4);
    A = [0 0.4 0 0; 0 0 0.345 0; 0 -0.524/d -0.465/d 0.262/d; 0 0 0 -1/d];
    B = [0; 0; 0; 1/d];
    C = [1 0 0 0; 0 0 1 0];
    D = zeros(2,1);

  elseif nr(2) == 2,
    %  Example 2.2:  Arnold/Laub 1984
    if nargin < 2, d = 10^(-6);  d = 10^(-6);  else d = parin(1); end

    E = eye(4);
    A = [-d 1 0 0; -1 -d 0 0; 0 0 d 1; 0 0 -1 d];
    B = ones(4,1);
    C = ones(1,4);
    D = 0;

  elseif nr(2) == 3,
    %  Example 2.3:  Vertical acceleration of a rigid guided missile
    if nargin < 2, s = 1;  else s = parin(1); end
    if s ~= round(s) | s < 1 | s > 10, error('Invalid value of parameter s.'); end
    PM = [-0.327  0.391 -0.688 -0.738 -0.886 -1.364 -0.333 -0.337 -0.369 -0.402;
          -63.94 -130.29 -619.27 -651.57 -1068.85 -92.82 -163.24 -224.03 -253.71 -277.2;
          -1.0 -1.42 -2.27 -2.75 -3.38 -4.68 -0.666 -0.663 -0.80 -0.884;
          -155.96 -186.5 -552.9 -604.18 -1004.39 -128.46 -153.32 -228.72 -249.87 -419.35;
          -0.237 -0.337 -0.429 -0.532 -0.582 -0.087 -0.124 -0.112 -0.135 -0.166;
           0.326 0.35 0.65 0.66 0.79 1.36 0.298 0.319 0.33 0.36;
          -208.5 -272.38 -651.11 -913.64 -1926.45 -184.26 -247.75 -375.75 -500.59 -796.18;
           90.93 75.06 283.44 250.5 402.96 76.43 63.77 117.4 103.76 178.59];

    E = eye(3);
    A = [PM(1,s) PM(2,s) PM(4,s);1 PM(3,s) PM(5,s);0 0 -190];
    B = [0; 0; 190 ];
    C = [PM(6,s) PM(7,s) PM(8,s) ];
    D = 0;

  elseif nr(2) == 4,
    %  Example 2.4:  Senning 1980: hydraulic positioning system
    if nargin < 2,
      parin = [14000 0.1287 0.15 0.01 0.002 0.24 10.75];
    elseif length(parin) < 7,
      error('Not enough parameters.');
    end;
    b     = parin(1);
    mu    = parin(2);
    r     = parin(3);
    rc    = parin(4);
    kl    = parin(5);
    sigma = parin(6);
    a     = parin(7);

    if ~(9000 < b & b < 16000),        error('Invalid value of parameter 1.'); end
    if ~(0.05 < mu & mu < 0.3),        error('Invalid value of parameter 2.'); end
    if ~(0.05 < r & r < 5),            error('Invalid value of parameter 3.'); end
    if ~(0 < rc & rc < 0.05),          error('Invalid value of parameter 4.'); end
    if ~(1.03e-04 < kl & kl < 0.0035), error('Invalid value of parameter 5.'); end
    if ~(0.001 < sigma & sigma < 15),  error('Invalid value of parameter 6.'); end
    if ~(10.5 < a & a < 11.10),        error('Invalid value of parameter 7.'); end
    v = 874;

    E = eye(3);
    A = [0 1 0; 0 (-r-4*rc/pi)/mu a/mu; 0 -4*b*a/v -4*b*(sigma+kl)/v];
    B = [0; 0; -4*b/v];
    C = [1 0 0];
    D = 0;

  elseif nr(2) == 5,
    %  Example 2.5:  Kwakernaak/Westdyk 1985: cascade of inverted pendula
    if nargin < 2, s = 1;  else s = parin(1); end
    if s ~= round(s) | s < 1 | s > 7, error('Invalid value of parameter s.'); end

    if s == 1
      m_s = 1;
      A = [0   1; 
           9.8 0];
      B = [0; 1];
      C = [1  0];

    elseif s == 2
      m_s = 2;
      A = [0   1  0   0;
           9.8 0 -9.8 0;
           0   0  0   1;
          -9.8 0 29.4 0];
      B = [0 0; 1 -2; 0 0; -2 5];
      C = [1 0 0 0; 0 0 1 0];
  
    elseif s == 3  
      m_s = 3;
      A = [zeros(3) eye(3);
           29.4 -19.6 -3.6285e-16 0 0 0;
          -29.4  39.2 -9.8        0 0 0;
            0   -19.6 19.6        0 0 0];
      B = [zeros(3);
           1.6667  -2.6667  1;
          -2.3333   4.8333 -3.5;
           0.66667 -2.6667  4];
      C = [eye(3) zeros(3)];

    elseif s == 4
      m_s = 4;
      A = [zeros(4) eye(4);
           [39.2 -29.4 0 0; -39.2 58.8 -19.6 -3.6285e-16; 
            -4.3521e-15 -29.4  39.2 -9.8; 0 0 -19.6 19.6] zeros(4)];
      B = [zeros(4);
           1.75      -2.75        1       6.4185e-17;
          -2.5        5.1667     -3.6667  1;
           0.75      -3.0833      4.8333 -3.5;
           3.5586e-17 6.6667e-01 -2.6667  4];
      C = [eye(4) zeros(4)];

    elseif s == 5
      m_s = 5;
      A = [zeros(5) eye(5);
           [49 -39.2 -4.5686e-15 1.0878e-15 0;
           -49 78.4 -29.4 0 0;
           -3.6630e-14 -39.2 58.8 -19.6 -3.6285e-16;
            1.2240e-14 -1.3056e-14 -29.4 39.2 -9.8;
            0 0 0 -19.6 19.6] zeros(5)];
      B = [zeros(5);
           1.8 -2.8 1 3.1160e-16 -2.2755e-16;
          -2.6 5.35 -3.75 1 3.0015e-16;
           0.8 -3.3 5.1667 -3.6667 1;
           5.1282e-16 0.75 -3.0833 4.8333 -3.5;
           5.3006e-17 -1.2003e-17 0.66667 -2.6667 4];
      C = [eye(5) zeros(5)];

    elseif s == 6  
      m_s = 6;
      A = [zeros(6) eye(6);
           [58.8 -49 2.9028e-15 -1.6328e-15  5.4427e-16  -1.0086e-32;
            -58.8 98 -39.2 -4.8961e-15 5.4401e-16 2.0127e-32;
            0 -49 78.4 -29.4 0 0;
            -2.0310e-14 2.0491e-14 -39.2 58.8 -19.6 -3.6285e-16;
            1.1968e-14 -1.5504e-14 5.4401e-15 -29.4 39.2 -9.8; 
            2.1755e-15 -1.6323e-15 2.1760e-15 -1.0885e-15 -19.6 19.6] zeros(6)];
      B = [zeros(6);
           1.8333 -2.8333 1 -3.0442e-16 1.2911e-16 -1.0963e-16;
          -2.6667 5.4667 -3.8 1 6.7758e-17 4.2643e-16;
           0.83333 -3.4333 5.35 -3.75 1 -7.5488e-16;
          -5.5618e-16 0.8 -3.3 5.1667 -3.6667 1;
           4.3692e-16 -9.2492e-16 0.75 -3.0833 4.8333 -3.5;
           1.1031e-16 -1.6825e-16  2.2085e-17 0.66667 -2.6667 4];
      C = [eye(6) zeros(6)];

    elseif s == 7  
      m_s = 10;
      A = [zeros(10) eye(10)
           [98 -88.2 -1.7383e-15 -3.7188e-18 -6.5249e-15 5.4388e-15 1.3084e-31 -1.6316e-15 5.4388e-16  -1.0042e-32;
           -98 176.4 -78.4 0 0 0 0 0 0 0;
           -1.0880e-14 -88.2 156.8 -68.6 6.5281e-15 -1.0880e-14 4.3521e-15 1.4495e-30 -4.8318e-31   8.9449e-48;
            1.5855e-14 -3.6370e-14 -78.4 137.2 -58.8 3.4194e-14 -1.9585e-14 3.2649e-15 -1.0883e-15 2.0172e-32;
           -1.1062e-13 1.0699e-13 -2.5751e-14 -68.6 117.6 -49 1.5959e-14 -4.8969e-15 1.6323e-15 4.0225e-32;
            1.0880e-13 -9.7922e-14 1.3056e-14 -7.6177e-15 -58.8 98 -39.2 -1.3056e-14 4.3521e-15 -8.0569e-32;
           -8.1596e-15 6.5289e-15 -8.1601e-15 -5.4454e-16 3.1009e-14 -49 78.4 -29.4 -4.8961e-15 8.0569e-32;
           -4.9868e-14 4.2342e-14 3.9892e-15 -1.4506e-15 -1.8859e-14  2.9013e-14 -39.2 58.8 -19.6 3.6285e-16;
            2.7200e-14 -2.3256e-14 0 0 8.1601e-15 -9.5202e-15 -4.3521e-15 -29.4 39.2 -9.8;
            0 0 0 0 0 0 0 0 -19.6 19.6] zeros(10)];
      B = [zeros(10);
           1.9 -2.9 1 2.1959e-16 -2.8176e-16 4.8626e-16 -1.0559e-16 -2.5389e-16 7.9671e-17 3.8476e-17;
          -2.8 5.6889 -3.8889 1 -8.3224e-16 3.5204e-16 -3.3739e-16 5.2418e-16 8.8172e-17 -2.6452e-16;
           0.9 -3.6778 5.6528 -3.875 1 -2.5226e-15 1.5710e-15 -8.8962e-16 -9.7319e-17 2.9196e-16;
          -5.1230e-16 0.88889 -3.6389 5.6071 -3.8571 1 -3.6967e-15 1.8224e-15 -2.0082e-16 4.7207e-17;
          -2.3590e-15 3.7146e-15 0.875 -3.5893 5.5476 -3.8333 1 -2.3237e-15 4.6098e-16 -5.5013e-16;
           2.6600e-15 -4.0482e-15 2.3538e-15 0.85714 -3.5238 5.4667 -3.8 1 2.8896e-16 1.0290e-15;
          -5.6230e-16 8.7784e-16 -5.4097e-16 -7.1966e-16 0.83333 -3.4333 5.35 -3.75 1 -1.6135e-15;
          -8.0494e-16 1.1525e-15 -2.2138e-16 3.3885e-16 -1.0433e-15 0.8 -3.3 5.1667 -3.6667 1;
           5.2901e-16 -7.6095e-16 1.7267e-16 -1.0799e-16 4.5900e-16 -1.2346e-16 0.75 -3.0833 4.8333 -3.5;
          -3.7739e-17 1.5780e-17 5.5573e-17 -9.0729e-18 -4.6091e-17 7.8485e-17 -2.2627e-16 0.66667 -2.6667 4];
      C = [eye(10) zeros(10)];
    end
    E = eye(2*m_s);
    D = zeros(m_s);

  elseif nr(2) == 6,
    %  Example 2.6:  Kallstrom/Astrom 1981: regulation of a ship heading
    if nargin < 2, s = 1;  else s = parin(1); end
    if s ~= round(s) | s < 1 | s > 5, error('Invalid value of parameter s.'); end
    PM = [-0.895     -0.770    -0.597       -0.298     -0.454;
          -0.286     -0.335    -0.372       -0.279     -0.433;
          -4.367     -3.394    -3.651       -4.370     -4.005;
          -0.918     -1.627    -0.792       -0.773     -0.807;
           0.108      0.170     0.103        0.116      0.097;   
          -0.918     -1.627    -0.792       -0.773     -0.807];

    E = eye(3);
    A = [PM(1,s) PM(2,s)      0;
         PM(3,s) PM(4,s)      0;
                0      1      0];
    B = [PM(5,s); PM(6,s);    0];
    C = [0 0 1];
    D = 0;

  elseif nr(2) == 7,
    %  Example 2.7:  Ackermann 1989: track-guided bus
    if nargin < 2,
      parin = [15 10];
    elseif length(parin) < 2, error('Not enough input parameters.'); end;
    mu = parin(1);
    nu = parin(2);
    if ~(9.95 <= mu & mu <= 16),  error('Invalid value of parameter 1.'); end
    if ~(1    <= nu & nu <= 20),  error('Invalid value of parameter 2.'); end
    l_F     = 3.67;
    l_R     = 1.93;
    delta_F = 99;
    delta_R = 235;
    c = 10.86;

    E = eye(5);
    A = [-2 * (delta_F+delta_R) / (mu*nu), ...
         -1 - 2 * (delta_F*l_F - delta_R*l_R) / (nu^2*mu) 0 0, ...
          2 * delta_F / (mu*nu); ...
         -2 * (delta_F*l_F - delta_R*l_R) / (mu*c), ...
         -2 * (delta_F*l_F^2 + delta_R*l_R^2) / (mu*nu*c) 0 0, ...
          2 * delta_F*l_F / (mu*c); ...
         nu 0 0 nu 0; ...
          0 1 0 0 0; ...
          0 0 0 0 0 ];
    B = [0; 0; 0; 0; 1];
    C = [0 0 1 6.12 0];
    D = 0;

  else
    str = sprintf('Example #%i is not available in Group #%i !',nr(2),nr(1));
    error(str);
  end;


elseif nr(1) == 3,
  if nr(2) ==1,
    %  Example 3.1:  Laub 1979, Ex.4: string of high speed vehicles
    if nargin < 2, l = 20;  else, l = parin(1); end 
    if l ~= round(l) | l < 2, error('Invalid value of parameter N.'); end
    n = 2*l - 1;
    
    E = eye(n);
    A = zeros(n);
    B = zeros(n,l);
    C = zeros(l-1,n);
    D = zeros(l-1,l);
    for i=1:l-1
      j = 2*i;
      k = j-1;  
      A(k,k)   = -1;
      A(j,k)   = 1;
      A(j,j+1) = -1;
      B(k,i)   = 1;
      C(i,j)   = 1;
    end
    A(n,n) = -1;
    B(n,l) = 1;

  elseif nr(2) == 2,
    %  Example 3.2:  heat flow in a thin rod
    if nargin < 2, n = 100;  else, n = parin(1); end 
    if n ~= round(n) | n < 1, error('Invalid value of parameter n.'); end
    h = n+1;

    E = eye(n);
    A = diag([-h -2*h*ones(1,n-1)]) + diag(h*ones(1,n-1),1) + diag(h*ones(1,n-1),-1);
    B = [zeros(n-1,1); h];
    C = eye(n);
    D = zeros(n,1);

  elseif nr(2) == 3,
    %  Example 3.3:  Laub 1979, Ex.6
    if nargin < 2, n = 21;  else, n = parin(1); end 
    if n ~= round(n) | n < 1, error('Invalid value of parameter n.'); end

    E = eye(n);
    if n > 1,  A = diag(ones(n-1,1),1);  else,  A = 0;  end
    B = flipud(eye(n,1));
    C = eye(1,n);
    D = 0;

  elseif nr(2) == 4,
    %  Example 3.4:  Lang/Penzl 1994: rotating axle
    mu = 0; gamma = 0; delta = 0; kappa = 0;
    load ctds304
    if nargin < 2, l = 211;  else l = parin(1); end
    if l ~= round(l) | l < 3 | l > 211, error('Invalid value of parameter l.'); end
    n = 2*l-1;
    E = zeros(n);
    E(1,1) = mu(1);
    for i = 2:l
      E(i,i-1) = mu(i-1);
      E(i,i) = -mu(i);
    end
    E(1:l,l) = -E(1:l,l);
    for i = l-1:-1:2
      E(1:l,i) = -E(1:l,i)+E(1:l,i+1);
    end
    E(1:l,1) = E(1:l,1)-E(1:l,2);
    E(l+1:n,l+1:n) = eye(l-1);
    A          = zeros(n);
    A(1,2)     = -gamma(1);
    A(2:l,2:l) = -(2*diag(gamma(1:l-1)) - diag(gamma(2:l-1),1) - ...
                   diag(gamma(1:l-2),-1));
    A(1:l,1:l) = A(1:l,1:l) - diag(delta(1:l));
    A(2:l,1)   = delta(2:l) - delta(1:l-1);
    for i = 2:l-1,
      A(i+1:l,i) = A(i+1:l,i) + delta(i:l-1) - delta(i+1:l);
    end
    A(1,l+1) = -kappa(1);
    A(2:l,l+1:n) = -(2*diag(kappa(1:l-1)) - diag(kappa(2:l-1),1) - ...
                     diag(kappa(1:l-2),-1));
    A(l+1:n,:) = [zeros(l-1,1), eye(l-1), zeros(l-1)];
    B = [diag([1;-ones(l-1,1)]) + diag(ones(l-1,1),-1); zeros(l-1,l)];
    C = [diag([1; gamma(1:l-1)]), [zeros(1,l-1); diag(kappa(1:l-1))]]; 
    D = zeros(l);

  else
    str = sprintf('Example #%i is not available in Group #%i !',nr(2),nr(1));
    error(str);
  end;

elseif nr(1) == 4,
  if nr(2) == 1,
    %  Example 4.1:  Rosen/Wang 1995: control of 1-dim. heat flow
    if nargin < 2,
      parin = [100 0.01 1 1 0.2 0.3 0.2 0.3];
    elseif length(parin) < 8, error('Not enough input parameters.'); end;
    n     = parin(1);
    alpha = parin(2);
    beta  = parin(3);
    gamma = parin(4);
    b_int = parin(5:6);
    c_int = parin(7:8);
    if n ~= round(n) | n < 2,  error('Invalid value of parameter n.'); end

    evec = ones(n-1,1);
    % Gram matrix
    E = (diag(evec,-1) + diag(4*ones(n,1)) + diag(evec,1)) / (6*(n+1));
    % generate system matrices
    A = (alpha*(n+1)*(diag(evec,-1) - diag(2*ones(n,1)) + diag(evec,1)));
    b = zeros(n,1);
    C = zeros(1,n);
    % evaluate integrals
    for i = 1:n,
      b1 = max((i-1)/(n+1),b_int(1));	b2 = min((i+1)/(n+1),b_int(2));
      c1 = max((i-1)/(n+1),c_int(1));	c2 = min((i+1)/(n+1),c_int(2));
      if b1 >= b2,   
        b(i) = 0;
      else
        b(i) = b2 - b1;
        temp = min(b2,i/(n+1));
        if b1 < temp,
          b(i) = b(i) + (n+1)*(temp^2 - b1^2)/2 + i*(b1 - temp);
        end
        temp = max(b1,i/(n+1));
        if temp < b2,
          b(i) = b(i) - (n+1)*(b2^2 - temp^2)/2 - i*(temp - b2);
        end
        b(i) = beta * b(i);
      end 
      if c1 >= c2, 
        C(i) = 0;
      else
        C(i) = c2 - c1;
        temp = min(c2,i/(n+1));
        if c1 < temp,
          C(i) = C(i) + (n+1)*(temp^2 - c1^2)/2 + i*(c1 - temp);
        end
        temp = max(c1,i/(n+1));
        if temp < c2,
          C(i) = C(i) - (n+1)*(c2^2 - temp^2)/2 - i*(temp - c2);
        end
        C(i) = gamma * C(i);
      end 
    end
    B = b;
    D = 0;

  elseif nr(2) == 2
    %  Example 4.2:  Hench et al. 1995: coupled springs, dashpots, masses
    if nargin < 2,
      parin = [30 4 4 1];
    elseif length(parin) < 4, error('Not enough input parameters.'); end;
    l     = parin(1);
    mu    = parin(2);
    delta = parin(3);
    kappa = parin(4);
    if l ~= round(l) | l < 2,  error('Invalid value of parameter l.'); end
    n = 2*l;
    m = 2;
    p = 2*l;

    E = [ eye(l), zeros(l); zeros(l), mu*eye(l) ];
    A = [zeros(l), eye(l);
        -kappa*(diag([1; 2*ones(l-2,1); 1]) - diag(ones(l-1,1),1) - ...
         diag(ones(l-1,1),-1)), -delta*eye(l) ];
    B = [ zeros(l,m); eye(l,1), flipud(-eye(l,1))];
    C = eye(n);
    D = zeros(p,2);

  else
    str = sprintf('Example #%i is not available in Group #%i !',nr(2),nr(1));
    error(str);
  end;

else
  str = sprintf('Group #%i is not available !',nr(1));
  error(str);
end;
