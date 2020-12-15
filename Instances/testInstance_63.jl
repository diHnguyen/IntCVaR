edge = [1 6; 1 7; 2 4; 2 9; 2 12; 2 13; 3 4; 3 9; 3 10; 3 11; 3 12; 3 13; 3 15; 4 5; 4 7; 4 11; 4 12; 4 13; 5 3; 5 4; 5 8; 5 9; 5 11; 5 12; 5 13; 6 12; 7 8; 7 11; 7 14; 8 3; 8 6; 8 7; 8 10; 8 14; 9 4; 9 8; 9 11; 9 13; 9 14; 9 15; 10 9; 10 11; 10 12; 10 13; 10 15; 11 2; 11 4; 11 10; 11 12; 11 13; 11 15; 12 2; 12 5; 13 5; 13 9; 13 10; 14 9; 14 11; 14 12; 14 13]
cL_orig = [118.0, 25.0, 80.0, 260.0, 113.0, 153.0, 22.0, 42.0, 212.0, 395.0, 61.0, 207.0, 17.0, 30.0, 46.0, 93.0, 183.0, 386.0, 88.0, 31.0, 22.0, 168.0, 167.0, 0.0, 139.0, 266.0, 5.0, 112.0, 287.0, 33.0, 20.0, 33.0, 95.0, 28.0, 206.0, 26.0, 60.0, 11.0, 156.0, 160.0, 1.0, 46.0, 62.0, 36.0, 45.0, 154.0, 44.0, 37.0, 40.0, 69.0, 88.0, 322.0, 77.0, 268.0, 151.0, 15.0, 41.0, 116.0, 60.0, 6.0]
cU_orig = [118.0, 133.0, 80.0, 260.0, 113.0, 153.0, 22.0, 122.0, 212.0, 395.0, 245.0, 207.0, 669.0, 52.0, 80.0, 541.0, 183.0, 386.0, 88.0, 31.0, 22.0, 168.0, 167.0, 4.0, 317.0, 266.0, 5.0, 138.0, 287.0, 33.0, 62.0, 33.0, 95.0, 40.0, 206.0, 26.0, 60.0, 11.0, 156.0, 160.0, 1.0, 46.0, 62.0, 68.0, 175.0, 154.0, 44.0, 37.0, 40.0, 69.0, 88.0, 322.0, 573.0, 268.0, 151.0, 15.0, 69.0, 116.0, 60.0, 32.0]
d = [14.0, 4.0, 40.0, 45.0, 25.0, 2.0, 44.0, 40.0, 49.0, 20.0, 36.0, 3.0, 49.0, 13.0, 28.0, 6.0, 40.0, 19.0, 19.0, 9.0, 43.0, 45.0, 26.0, 34.0, 32.0, 12.0, 5.0, 22.0, 24.0, 47.0, 42.0, 43.0, 17.0, 42.0, 23.0, 43.0, 31.0, 30.0, 44.0, 41.0, 31.0, 6.0, 42.0, 39.0, 17.0, 35.0, 21.0, 3.0, 31.0, 4.0, 32.0, 6.0, 43.0, 18.0, 40.0, 26.0, 17.0, 43.0, 18.0, 40.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1]

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
