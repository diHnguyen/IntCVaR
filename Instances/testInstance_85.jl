edge = [1 3; 1 4; 1 6; 1 12; 2 3; 2 5; 2 11; 2 13; 3 5; 3 6; 3 7; 3 10; 3 13; 4 5; 4 6; 4 7; 4 9; 4 14; 5 4; 5 6; 5 7; 5 8; 5 13; 5 14; 6 2; 6 8; 6 12; 7 9; 7 11; 7 14; 7 15; 8 13; 8 15; 9 2; 9 5; 9 10; 9 12; 9 14; 10 7; 10 15; 11 2; 11 4; 11 8; 11 13; 12 2; 12 3; 12 6; 12 8; 12 9; 13 3; 13 6; 13 9; 13 10; 13 11; 13 12; 14 2; 14 5; 14 7; 14 9; 14 10; 14 11; 14 13; 14 15]
cL_orig = [88.0, 69.0, 8.0, 255.0, 40.0, 39.0, 42.0, 286.0, 2.0, 24.0, 11.0, 275.0, 476.0, 19.0, 39.0, 21.0, 0.0, 64.0, 49.0, 8.0, 54.0, 81.0, 107.0, 86.0, 152.0, 63.0, 201.0, 40.0, 40.0, 126.0, 84.0, 87.0, 25.0, 96.0, 13.0, 44.0, 113.0, 96.0, 52.0, 10.0, 25.0, 286.0, 10.0, 55.0, 210.0, 288.0, 60.0, 70.0, 50.0, 463.0, 257.0, 196.0, 64.0, 76.0, 29.0, 383.0, 250.0, 281.0, 72.0, 165.0, 135.0, 2.0, 32.0]
cU_orig = [88.0, 73.0, 56.0, 255.0, 44.0, 65.0, 158.0, 286.0, 2.0, 24.0, 97.0, 275.0, 476.0, 19.0, 39.0, 21.0, 2.0, 64.0, 49.0, 20.0, 140.0, 81.0, 351.0, 86.0, 152.0, 63.0, 201.0, 40.0, 40.0, 190.0, 84.0, 149.0, 25.0, 96.0, 13.0, 44.0, 113.0, 96.0, 52.0, 10.0, 25.0, 286.0, 16.0, 55.0, 210.0, 288.0, 60.0, 70.0, 50.0, 463.0, 257.0, 196.0, 64.0, 76.0, 29.0, 383.0, 250.0, 281.0, 278.0, 165.0, 135.0, 22.0, 64.0]
d = [32.0, 13.0, 36.0, 13.0, 43.0, 11.0, 18.0, 18.0, 45.0, 9.0, 29.0, 24.0, 35.0, 22.0, 40.0, 40.0, 41.0, 33.0, 12.0, 31.0, 25.0, 36.0, 31.0, 43.0, 45.0, 4.0, 42.0, 34.0, 41.0, 46.0, 29.0, 47.0, 7.0, 19.0, 39.0, 37.0, 38.0, 20.0, 18.0, 30.0, 28.0, 16.0, 3.0, 16.0, 15.0, 17.0, 25.0, 4.0, 16.0, 27.0, 32.0, 41.0, 37.0, 26.0, 22.0, 48.0, 1.0, 7.0, 3.0, 40.0, 18.0, 33.0, 50.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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