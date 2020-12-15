edge = [1 6; 1 15; 2 3; 2 7; 2 10; 2 12; 3 4; 3 8; 4 3; 4 5; 4 9; 4 10; 5 2; 5 7; 5 8; 5 9; 5 12; 6 2; 6 3; 6 7; 6 10; 6 11; 6 12; 7 5; 7 9; 7 13; 7 14; 8 11; 8 13; 8 15; 9 2; 9 4; 9 8; 9 13; 10 3; 10 6; 10 8; 10 9; 10 13; 10 14; 10 15; 11 2; 11 3; 11 9; 11 13; 12 2; 12 3; 12 5; 12 6; 12 8; 12 10; 12 13; 12 14; 12 15; 13 3; 13 4; 13 6; 13 9; 13 10; 13 11; 13 14; 14 7; 14 9; 14 13]
cL_orig = [88.0, 8.0, 28.0, 190.0, 162.0, 168.0, 22.0, 111.0, 43.0, 10.0, 209.0, 185.0, 6.0, 90.0, 33.0, 15.0, 76.0, 33.0, 25.0, 17.0, 119.0, 4.0, 272.0, 29.0, 42.0, 114.0, 22.0, 17.0, 91.0, 216.0, 146.0, 201.0, 32.0, 78.0, 322.0, 59.0, 42.0, 20.0, 100.0, 138.0, 52.0, 367.0, 310.0, 0.0, 66.0, 477.0, 115.0, 286.0, 1.0, 148.0, 23.0, 45.0, 32.0, 93.0, 144.0, 188.0, 103.0, 3.0, 3.0, 32.0, 27.0, 32.0, 12.0, 11.0]
cU_orig = [88.0, 8.0, 28.0, 220.0, 162.0, 444.0, 22.0, 111.0, 43.0, 36.0, 209.0, 185.0, 6.0, 90.0, 33.0, 15.0, 222.0, 33.0, 119.0, 17.0, 119.0, 4.0, 272.0, 35.0, 42.0, 144.0, 46.0, 37.0, 139.0, 216.0, 146.0, 201.0, 32.0, 78.0, 322.0, 63.0, 42.0, 20.0, 100.0, 138.0, 204.0, 367.0, 310.0, 6.0, 66.0, 477.0, 679.0, 286.0, 103.0, 148.0, 23.0, 45.0, 32.0, 93.0, 236.0, 188.0, 103.0, 3.0, 41.0, 46.0, 27.0, 126.0, 12.0, 11.0]
d = [47.0, 8.0, 40.0, 9.0, 26.0, 46.0, 6.0, 36.0, 8.0, 43.0, 1.0, 16.0, 34.0, 18.0, 23.0, 2.0, 2.0, 26.0, 33.0, 8.0, 7.0, 47.0, 3.0, 45.0, 23.0, 4.0, 28.0, 48.0, 1.0, 18.0, 4.0, 47.0, 8.0, 13.0, 50.0, 2.0, 13.0, 32.0, 27.0, 43.0, 10.0, 11.0, 24.0, 24.0, 15.0, 38.0, 8.0, 16.0, 23.0, 24.0, 34.0, 11.0, 36.0, 14.0, 13.0, 41.0, 1.0, 20.0, 50.0, 33.0, 49.0, 27.0, 7.0, 19.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
