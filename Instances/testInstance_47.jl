edge = [1 2; 1 6; 1 7; 1 10; 1 11; 1 13; 1 15; 2 3; 2 8; 2 12; 3 9; 3 11; 3 13; 3 15; 4 6; 4 7; 4 12; 4 14; 4 15; 5 3; 5 6; 5 8; 5 15; 6 14; 6 15; 7 3; 7 4; 7 8; 7 9; 7 10; 8 9; 8 14; 9 4; 9 5; 9 12; 9 15; 10 2; 10 6; 10 7; 10 15; 11 2; 11 3; 11 4; 11 7; 11 12; 11 13; 11 14; 12 13; 12 14; 13 6; 13 7; 13 8; 13 10; 13 11; 13 15; 14 2; 14 4; 14 8; 14 10; 14 13]
cL_orig = [5.0, 150.0, 38.0, 278.0, 346.0, 314.0, 41.0, 43.0, 141.0, 315.0, 2.0, 217.0, 205.0, 87.0, 21.0, 4.0, 84.0, 74.0, 145.0, 33.0, 23.0, 3.0, 143.0, 56.0, 280.0, 21.0, 111.0, 47.0, 10.0, 18.0, 22.0, 133.0, 242.0, 189.0, 118.0, 271.0, 32.0, 47.0, 76.0, 139.0, 40.0, 184.0, 21.0, 39.0, 26.0, 10.0, 8.0, 0.0, 78.0, 50.0, 187.0, 27.0, 43.0, 65.0, 78.0, 591.0, 148.0, 126.0, 197.0, 23.0]
cU_orig = [19.0, 150.0, 38.0, 278.0, 346.0, 314.0, 41.0, 43.0, 141.0, 315.0, 24.0, 217.0, 205.0, 87.0, 61.0, 4.0, 84.0, 74.0, 785.0, 33.0, 23.0, 127.0, 143.0, 56.0, 280.0, 21.0, 111.0, 47.0, 54.0, 26.0, 22.0, 133.0, 242.0, 189.0, 118.0, 271.0, 32.0, 47.0, 76.0, 139.0, 130.0, 362.0, 53.0, 197.0, 58.0, 18.0, 10.0, 76.0, 78.0, 488.0, 187.0, 27.0, 43.0, 99.0, 78.0, 591.0, 148.0, 216.0, 197.0, 23.0]
d = [44.0, 8.0, 24.0, 22.0, 8.0, 31.0, 13.0, 46.0, 46.0, 4.0, 16.0, 39.0, 32.0, 36.0, 2.0, 33.0, 30.0, 3.0, 14.0, 29.0, 14.0, 21.0, 40.0, 1.0, 29.0, 30.0, 36.0, 44.0, 17.0, 12.0, 40.0, 48.0, 28.0, 34.0, 18.0, 23.0, 24.0, 48.0, 17.0, 10.0, 44.0, 9.0, 13.0, 28.0, 35.0, 5.0, 9.0, 29.0, 37.0, 36.0, 2.0, 15.0, 50.0, 2.0, 44.0, 29.0, 28.0, 19.0, 42.0, 24.0]
Len = length(d)

yy = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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