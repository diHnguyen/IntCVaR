
edge = [1 8; 1 12; 1 13; 2 9; 2 11; 3 5; 3 8; 3 12; 3 14; 3 15; 4 1; 4 3; 4 6; 4 7; 4 14; 5 1; 5 2; 5 3; 5 7; 5 13; 6 8; 6 10; 6 13; 6 14; 6 15; 7 15; 8 10; 8 15; 10 7; 10 15; 11 6; 11 7; 11 12; 12 6; 12 10; 13 3; 13 14; 14 1; 14 2; 14 6; 14 7; 14 11; 15 5; 15 8; 15 9; 15 13; 15 14]
cL_orig = [47.0, 84.0, 13.0, 41.0, 5.0, 9.0, 22.0, 7.0, 72.0, 65.0, 36.0, 49.0, 45.0, 36.0, 100.0, 77.0, 14.0, 81.0, 9.0, 67.0, 70.0, 92.0, 95.0, 91.0, 26.0, 75.0, 51.0, 35.0, 50.0, 82.0, 10.0, 32.0, 41.0, 93.0, 60.0, 40.0, 96.0, 35.0, 44.0, 97.0, 67.0, 7.0, 49.0, 2.0, 71.0, 92.0, 28.0]
cU_orig = [47.0, 84.0, 20.0, 41.0, 5.0, 9.0, 22.0, 7.0, 72.0, 98.0, 57.0, 49.0, 45.0, 45.0, 100.0, 77.0, 14.0, 81.0, 9.0, 86.0, 70.0, 92.0, 95.0, 91.0, 26.0, 75.0, 51.0, 35.0, 58.0, 100.0, 10.0, 32.0, 41.0, 93.0, 60.0, 40.0, 96.0, 46.0, 44.0, 97.0, 67.0, 7.0, 49.0, 2.0, 71.0, 92.0, 31.0]
d = [19.0, 51.0, 17.0, 20.0, 5.0, 5.0, 3.0, 5.0, 56.0, 46.0, 2.0, 49.0, 3.0, 33.0, 36.0, 59.0, 11.0, 81.0, 9.0, 42.0, 14.0, 56.0, 31.0, 67.0, 15.0, 63.0, 20.0, 5.0, 21.0, 96.0, 8.0, 25.0, 39.0, 79.0, 42.0, 21.0, 17.0, 24.0, 43.0, 27.0, 37.0, 3.0, 26.0, 2.0, 63.0, 81.0, 24.0]

Len = length(cL_orig)

c_orig = 0.5*(cL_orig+cU_orig)

#STARTING SOLUTION:
yy = zeros(Len)

yy[15] = 1
yy[39] = 1
yy[4] = 1

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)

p = [1.0]
g = [SP_init]
h = [0.0]

origin = 4
destination = 9

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 4.0
α = 1.0
last_node = maximum(edge)
