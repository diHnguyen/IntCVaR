edge = [1 6; 1 8; 1 9; 1 10; 1 11; 2 5; 2 6; 2 9; 2 14; 3 5; 3 7; 3 10; 3 11; 3 12; 4 5; 4 7; 4 10; 4 11; 4 12; 4 13; 5 2; 5 6; 5 8; 5 10; 5 14; 6 5; 6 8; 6 9; 6 10; 6 13; 6 15; 7 4; 7 6; 7 9; 7 10; 7 11; 8 6; 8 11; 8 12; 8 14; 9 4; 9 11; 9 13; 10 8; 10 11; 10 12; 10 14; 11 3; 11 7; 11 10; 11 12; 11 15; 12 6; 12 13; 12 14; 13 5; 13 7; 14 2; 14 6; 14 9; 14 15]
cL_orig = [86.0, 55.0, 270.0, 442.0, 209.0, 135.0, 193.0, 290.0, 44.0, 51.0, 149.0, 126.0, 353.0, 213.0, 39.0, 4.0, 77.0, 320.0, 338.0, 51.0, 10.0, 5.0, 87.0, 21.0, 238.0, 1.0, 70.0, 25.0, 10.0, 98.0, 172.0, 17.0, 46.0, 15.0, 13.0, 25.0, 0.0, 107.0, 112.0, 261.0, 115.0, 95.0, 78.0, 29.0, 20.0, 18.0, 97.0, 175.0, 78.0, 14.0, 44.0, 193.0, 243.0, 38.0, 44.0, 44.0, 221.0, 312.0, 308.0, 172.0, 3.0]
cU_orig = [86.0, 55.0, 270.0, 442.0, 209.0, 135.0, 193.0, 290.0, 44.0, 51.0, 149.0, 126.0, 353.0, 213.0, 39.0, 8.0, 279.0, 320.0, 338.0, 475.0, 148.0, 23.0, 87.0, 81.0, 238.0, 79.0, 70.0, 25.0, 10.0, 110.0, 236.0, 17.0, 46.0, 15.0, 13.0, 261.0, 18.0, 139.0, 188.0, 261.0, 325.0, 95.0, 78.0, 29.0, 20.0, 42.0, 97.0, 175.0, 78.0, 38.0, 44.0, 193.0, 243.0, 38.0, 44.0, 44.0, 221.0, 312.0, 308.0, 326.0, 79.0]
d = [41.0, 13.0, 4.0, 15.0, 50.0, 20.0, 23.0, 48.0, 23.0, 40.0, 21.0, 32.0, 5.0, 3.0, 11.0, 15.0, 3.0, 42.0, 18.0, 46.0, 36.0, 40.0, 2.0, 30.0, 27.0, 12.0, 8.0, 32.0, 9.0, 23.0, 33.0, 18.0, 31.0, 23.0, 10.0, 48.0, 41.0, 40.0, 7.0, 39.0, 6.0, 5.0, 27.0, 43.0, 41.0, 33.0, 32.0, 40.0, 11.0, 23.0, 28.0, 22.0, 38.0, 40.0, 37.0, 37.0, 43.0, 20.0, 45.0, 14.0, 7.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]

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