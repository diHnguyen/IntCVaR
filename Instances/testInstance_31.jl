edge = [1 3; 1 7; 1 9; 1 10; 1 11; 1 12; 2 3; 2 4; 2 9; 2 10; 2 11; 3 2; 3 12; 3 14; 4 12; 5 4; 5 6; 5 13; 5 14; 5 15; 6 2; 6 3; 6 7; 6 10; 6 12; 6 13; 7 2; 7 5; 7 9; 8 2; 8 4; 8 5; 8 6; 8 12; 9 5; 9 8; 9 12; 10 3; 10 11; 10 12; 10 13; 10 15; 11 9; 11 14; 11 15; 12 4; 12 8; 12 9; 12 11; 13 2; 13 5; 13 6; 13 8; 13 15; 14 2; 14 3; 14 5; 14 8; 14 9; 14 13]
cL_orig = [8.0, 33.0, 123.0, 152.0, 363.0, 63.0, 18.0, 85.0, 196.0, 15.0, 87.0, 0.0, 314.0, 443.0, 37.0, 18.0, 6.0, 277.0, 236.0, 6.0, 155.0, 136.0, 3.0, 198.0, 289.0, 273.0, 238.0, 3.0, 13.0, 168.0, 72.0, 76.0, 89.0, 14.0, 169.0, 6.0, 141.0, 96.0, 36.0, 39.0, 10.0, 225.0, 3.0, 144.0, 87.0, 78.0, 190.0, 28.0, 39.0, 12.0, 24.0, 97.0, 7.0, 7.0, 324.0, 398.0, 73.0, 55.0, 227.0, 28.0]
cU_orig = [14.0, 105.0, 123.0, 152.0, 363.0, 131.0, 18.0, 85.0, 196.0, 15.0, 425.0, 6.0, 386.0, 443.0, 109.0, 18.0, 6.0, 277.0, 236.0, 6.0, 155.0, 136.0, 3.0, 198.0, 289.0, 273.0, 238.0, 3.0, 41.0, 168.0, 72.0, 76.0, 89.0, 18.0, 169.0, 6.0, 141.0, 96.0, 36.0, 73.0, 10.0, 225.0, 3.0, 144.0, 87.0, 78.0, 190.0, 262.0, 39.0, 12.0, 24.0, 97.0, 7.0, 27.0, 614.0, 398.0, 133.0, 205.0, 227.0, 70.0]
d = [23.0, 46.0, 17.0, 20.0, 21.0, 25.0, 42.0, 1.0, 44.0, 49.0, 21.0, 31.0, 42.0, 7.0, 16.0, 12.0, 5.0, 1.0, 40.0, 14.0, 3.0, 20.0, 7.0, 30.0, 20.0, 6.0, 19.0, 31.0, 21.0, 26.0, 45.0, 18.0, 13.0, 19.0, 38.0, 34.0, 20.0, 22.0, 22.0, 22.0, 49.0, 8.0, 49.0, 45.0, 33.0, 27.0, 15.0, 22.0, 40.0, 29.0, 3.0, 8.0, 7.0, 46.0, 17.0, 28.0, 41.0, 13.0, 43.0, 15.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]

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
