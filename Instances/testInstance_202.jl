edge = [1 3; 1 5; 1 9; 2 5; 2 12; 2 14; 2 15; 3 5; 3 6; 3 8; 3 11; 3 14; 3 15; 4 3; 4 5; 4 11; 4 12; 5 2; 5 6; 5 8; 5 11; 5 12; 5 14; 6 2; 6 4; 6 5; 6 7; 6 8; 6 12; 7 2; 7 3; 7 9; 7 13; 8 4; 8 7; 8 9; 8 11; 8 14; 9 2; 9 5; 9 13; 10 2; 10 3; 10 4; 10 7; 10 8; 10 12; 11 4; 11 5; 11 7; 11 13; 11 15; 12 4; 12 13; 12 14; 13 3; 13 6; 13 8; 13 11; 13 15; 14 2; 14 5; 14 7; 14 9; 14 10; 14 15]
cL_orig = [98.0, 53.0, 60.0, 12.0, 261.0, 28.0, 429.0, 5.0, 143.0, 206.0, 180.0, 530.0, 84.0, 7.0, 46.0, 184.0, 77.0, 68.0, 7.0, 48.0, 183.0, 263.0, 23.0, 17.0, 44.0, 38.0, 35.0, 10.0, 113.0, 162.0, 181.0, 82.0, 255.0, 122.0, 28.0, 15.0, 132.0, 185.0, 53.0, 108.0, 14.0, 98.0, 277.0, 251.0, 119.0, 56.0, 38.0, 168.0, 289.0, 7.0, 9.0, 123.0, 140.0, 6.0, 88.0, 228.0, 9.0, 138.0, 61.0, 71.0, 106.0, 44.0, 322.0, 186.0, 139.0, 13.0]
cU_orig = [98.0, 53.0, 132.0, 12.0, 261.0, 28.0, 429.0, 39.0, 143.0, 206.0, 180.0, 530.0, 84.0, 37.0, 46.0, 184.0, 95.0, 68.0, 7.0, 48.0, 183.0, 383.0, 23.0, 33.0, 44.0, 38.0, 35.0, 10.0, 113.0, 162.0, 181.0, 82.0, 255.0, 122.0, 50.0, 15.0, 132.0, 239.0, 55.0, 108.0, 14.0, 98.0, 277.0, 251.0, 119.0, 56.0, 38.0, 168.0, 289.0, 235.0, 191.0, 143.0, 148.0, 6.0, 88.0, 228.0, 415.0, 138.0, 61.0, 71.0, 106.0, 44.0, 322.0, 186.0, 139.0, 13.0]
d = [32.0, 34.0, 45.0, 38.0, 48.0, 31.0, 8.0, 11.0, 44.0, 19.0, 15.0, 35.0, 17.0, 43.0, 20.0, 1.0, 5.0, 1.0, 44.0, 42.0, 35.0, 11.0, 39.0, 20.0, 30.0, 31.0, 38.0, 2.0, 25.0, 6.0, 34.0, 33.0, 46.0, 16.0, 21.0, 43.0, 24.0, 14.0, 13.0, 18.0, 16.0, 43.0, 8.0, 9.0, 31.0, 33.0, 2.0, 40.0, 28.0, 29.0, 31.0, 4.0, 31.0, 40.0, 37.0, 38.0, 13.0, 20.0, 38.0, 31.0, 30.0, 36.0, 2.0, 17.0, 32.0, 19.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
