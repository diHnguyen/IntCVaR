edge = [1 6; 1 10; 1 12; 1 13; 1 15; 2 5; 2 8; 2 9; 2 10; 3 4; 3 7; 3 9; 3 11; 3 12; 4 2; 4 3; 4 7; 4 10; 4 15; 5 3; 5 6; 5 8; 5 9; 5 11; 5 12; 5 14; 6 2; 6 3; 6 4; 6 7; 6 9; 6 11; 7 2; 7 3; 7 5; 7 10; 7 11; 7 12; 7 14; 7 15; 8 2; 8 4; 8 10; 8 11; 8 14; 8 15; 9 2; 9 4; 9 8; 9 11; 9 12; 10 3; 10 6; 10 12; 10 13; 10 14; 10 15; 11 3; 11 7; 11 8; 11 9; 11 12; 11 13; 11 14; 12 4; 12 5; 12 8; 12 9; 13 10; 13 14; 13 15; 14 2; 14 6; 14 10; 14 11; 14 12]
cL_orig = [131.0, 122.0, 65.0, 17.0, 271.0, 3.0, 182.0, 17.0, 335.0, 49.0, 55.0, 165.0, 250.0, 219.0, 85.0, 33.0, 77.0, 10.0, 333.0, 58.0, 41.0, 83.0, 186.0, 104.0, 84.0, 91.0, 102.0, 68.0, 35.0, 17.0, 1.0, 185.0, 171.0, 41.0, 26.0, 90.0, 0.0, 109.0, 336.0, 134.0, 1.0, 165.0, 48.0, 25.0, 163.0, 76.0, 248.0, 73.0, 13.0, 46.0, 57.0, 124.0, 46.0, 89.0, 79.0, 67.0, 149.0, 331.0, 56.0, 41.0, 50.0, 41.0, 27.0, 97.0, 120.0, 329.0, 34.0, 142.0, 79.0, 15.0, 3.0, 164.0, 364.0, 70.0, 70.0, 40.0]
cU_orig = [131.0, 232.0, 65.0, 47.0, 271.0, 3.0, 182.0, 17.0, 335.0, 49.0, 85.0, 165.0, 250.0, 219.0, 85.0, 33.0, 149.0, 188.0, 665.0, 58.0, 41.0, 83.0, 186.0, 104.0, 84.0, 91.0, 280.0, 68.0, 35.0, 17.0, 5.0, 185.0, 171.0, 193.0, 26.0, 90.0, 22.0, 343.0, 336.0, 148.0, 47.0, 165.0, 98.0, 25.0, 213.0, 76.0, 248.0, 73.0, 13.0, 46.0, 57.0, 124.0, 46.0, 89.0, 79.0, 67.0, 149.0, 331.0, 56.0, 41.0, 50.0, 41.0, 27.0, 97.0, 120.0, 329.0, 34.0, 142.0, 79.0, 15.0, 3.0, 164.0, 364.0, 70.0, 70.0, 40.0]
d = [23.0, 20.0, 6.0, 40.0, 49.0, 31.0, 13.0, 49.0, 5.0, 15.0, 48.0, 49.0, 3.0, 45.0, 47.0, 33.0, 19.0, 5.0, 39.0, 10.0, 37.0, 47.0, 48.0, 12.0, 49.0, 32.0, 50.0, 13.0, 32.0, 18.0, 42.0, 7.0, 44.0, 38.0, 24.0, 32.0, 20.0, 15.0, 21.0, 50.0, 43.0, 3.0, 16.0, 46.0, 22.0, 1.0, 3.0, 28.0, 21.0, 22.0, 9.0, 36.0, 2.0, 25.0, 28.0, 46.0, 6.0, 4.0, 23.0, 16.0, 35.0, 25.0, 25.0, 11.0, 40.0, 48.0, 36.0, 33.0, 28.0, 39.0, 29.0, 50.0, 4.0, 39.0, 27.0, 13.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

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