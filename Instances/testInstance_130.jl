edge = [1 2; 1 7; 1 8; 1 9; 1 15; 2 3; 2 4; 2 7; 2 9; 2 11; 2 14; 2 15; 3 4; 3 7; 3 9; 3 10; 3 12; 4 2; 4 5; 4 6; 4 7; 4 10; 4 11; 4 13; 5 4; 5 8; 5 11; 5 13; 6 3; 6 5; 6 8; 6 10; 6 11; 6 15; 7 5; 7 8; 7 13; 8 2; 8 3; 8 5; 8 7; 8 14; 8 15; 9 2; 9 3; 9 4; 9 7; 9 11; 9 12; 10 5; 10 7; 10 11; 10 15; 11 2; 11 3; 11 4; 11 5; 11 10; 11 13; 12 10; 12 13; 12 14; 13 4; 13 5; 13 7; 13 12; 13 14; 14 2; 14 7; 14 9; 14 15]
cL_orig = [28.0, 274.0, 35.0, 96.0, 653.0, 48.0, 36.0, 42.0, 42.0, 108.0, 292.0, 280.0, 21.0, 83.0, 39.0, 0.0, 434.0, 32.0, 40.0, 38.0, 68.0, 18.0, 309.0, 35.0, 44.0, 102.0, 7.0, 179.0, 56.0, 25.0, 60.0, 167.0, 239.0, 7.0, 80.0, 43.0, 46.0, 148.0, 101.0, 3.0, 41.0, 63.0, 308.0, 33.0, 11.0, 134.0, 11.0, 51.0, 32.0, 107.0, 150.0, 2.0, 173.0, 168.0, 358.0, 4.0, 50.0, 10.0, 6.0, 79.0, 10.0, 97.0, 100.0, 54.0, 0.0, 22.0, 16.0, 18.0, 39.0, 1.0, 45.0]
cU_orig = [28.0, 274.0, 193.0, 96.0, 653.0, 48.0, 102.0, 42.0, 42.0, 108.0, 292.0, 280.0, 21.0, 83.0, 41.0, 86.0, 434.0, 158.0, 46.0, 38.0, 68.0, 24.0, 309.0, 59.0, 44.0, 102.0, 221.0, 179.0, 160.0, 25.0, 60.0, 167.0, 239.0, 7.0, 80.0, 43.0, 46.0, 148.0, 101.0, 3.0, 41.0, 63.0, 308.0, 653.0, 11.0, 134.0, 11.0, 51.0, 94.0, 107.0, 150.0, 2.0, 173.0, 168.0, 358.0, 14.0, 216.0, 10.0, 6.0, 79.0, 10.0, 97.0, 100.0, 388.0, 2.0, 22.0, 16.0, 164.0, 39.0, 1.0, 45.0]
d = [17.0, 36.0, 14.0, 4.0, 45.0, 3.0, 44.0, 18.0, 28.0, 6.0, 5.0, 15.0, 40.0, 6.0, 12.0, 48.0, 7.0, 23.0, 4.0, 48.0, 36.0, 32.0, 10.0, 42.0, 40.0, 21.0, 4.0, 29.0, 36.0, 46.0, 47.0, 28.0, 41.0, 10.0, 4.0, 16.0, 25.0, 34.0, 10.0, 7.0, 1.0, 22.0, 37.0, 47.0, 38.0, 25.0, 15.0, 1.0, 14.0, 44.0, 44.0, 38.0, 19.0, 50.0, 2.0, 42.0, 29.0, 5.0, 50.0, 36.0, 50.0, 15.0, 24.0, 11.0, 35.0, 42.0, 27.0, 4.0, 26.0, 24.0, 22.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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