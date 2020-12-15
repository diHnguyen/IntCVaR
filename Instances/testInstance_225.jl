edge = [1 7; 1 13; 2 4; 2 9; 3 2; 3 8; 3 10; 3 11; 3 13; 3 14; 4 2; 4 3; 4 7; 4 11; 5 3; 5 6; 6 8; 6 13; 6 15; 8 4; 8 7; 9 4; 9 5; 9 7; 9 8; 9 10; 9 13; 9 14; 10 3; 10 4; 10 5; 10 9; 10 12; 10 14; 10 15; 11 2; 11 4; 11 5; 11 10; 11 13; 11 14; 11 15; 12 2; 12 4; 12 15; 13 14; 14 4; 14 8; 14 13; 14 15]
cL_orig = [108.0, 247.0, 15.0, 336.0, 38.0, 74.0, 2.0, 161.0, 257.0, 485.0, 88.0, 49.0, 3.0, 154.0, 13.0, 10.0, 66.0, 96.0, 345.0, 56.0, 9.0, 231.0, 1.0, 90.0, 14.0, 46.0, 194.0, 43.0, 57.0, 3.0, 202.0, 50.0, 5.0, 5.0, 96.0, 242.0, 333.0, 33.0, 9.0, 18.0, 56.0, 10.0, 498.0, 314.0, 41.0, 45.0, 323.0, 230.0, 3.0, 1.0]
cU_orig = [108.0, 247.0, 147.0, 336.0, 38.0, 74.0, 386.0, 161.0, 257.0, 485.0, 88.0, 49.0, 3.0, 154.0, 13.0, 10.0, 66.0, 96.0, 345.0, 56.0, 9.0, 231.0, 3.0, 90.0, 14.0, 46.0, 194.0, 365.0, 57.0, 21.0, 202.0, 50.0, 89.0, 5.0, 96.0, 242.0, 333.0, 33.0, 81.0, 18.0, 56.0, 10.0, 498.0, 314.0, 53.0, 45.0, 323.0, 230.0, 13.0, 1.0]
d = [15.0, 38.0, 21.0, 17.0, 3.0, 9.0, 13.0, 15.0, 50.0, 6.0, 45.0, 24.0, 50.0, 50.0, 30.0, 47.0, 15.0, 11.0, 31.0, 37.0, 17.0, 3.0, 12.0, 7.0, 18.0, 36.0, 6.0, 36.0, 48.0, 31.0, 12.0, 47.0, 32.0, 14.0, 6.0, 24.0, 27.0, 8.0, 23.0, 43.0, 38.0, 6.0, 37.0, 41.0, 43.0, 48.0, 18.0, 48.0, 12.0, 38.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1]

c_orig = 0.5*(cL_orig+cU_orig)

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)


p = [1.0]

g = [SP_init]

h = [0.0]


origin = 1

destination =15

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 5
last_node = maximum(edge)
