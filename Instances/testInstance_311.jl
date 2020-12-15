edge = [1 3; 1 5; 1 7; 1 10; 1 14; 2 6; 2 15; 3 2; 3 5; 3 7; 3 12; 3 13; 4 2; 4 3; 4 5; 4 7; 4 8; 4 9; 4 13; 4 14; 5 11; 6 10; 6 12; 7 2; 7 3; 7 5; 7 6; 7 8; 7 14; 8 3; 8 5; 8 9; 8 10; 8 11; 8 12; 8 14; 8 15; 9 2; 9 6; 9 13; 9 14; 10 4; 10 5; 10 11; 10 15; 11 4; 11 6; 11 10; 11 12; 11 13; 11 14; 12 2; 12 5; 12 6; 12 8; 12 9; 12 13; 13 2; 13 4; 13 7; 13 8; 13 9; 13 11; 13 14; 13 15; 14 5; 14 7; 14 8; 14 9; 14 10; 14 13]
cL_orig = [53.0, 187.0, 229.0, 174.0, 411.0, 128.0, 356.0, 16.0, 81.0, 5.0, 323.0, 286.0, 48.0, 15.0, 0.0, 6.0, 87.0, 9.0, 280.0, 165.0, 88.0, 94.0, 203.0, 49.0, 16.0, 72.0, 48.0, 42.0, 126.0, 158.0, 120.0, 11.0, 99.0, 75.0, 32.0, 2.0, 245.0, 142.0, 17.0, 146.0, 222.0, 25.0, 249.0, 17.0, 18.0, 199.0, 94.0, 29.0, 13.0, 85.0, 15.0, 244.0, 96.0, 177.0, 57.0, 62.0, 0.0, 495.0, 97.0, 148.0, 103.0, 56.0, 34.0, 1.0, 43.0, 94.0, 251.0, 67.0, 64.0, 25.0, 32.0]
cU_orig = [53.0, 187.0, 229.0, 174.0, 411.0, 128.0, 356.0, 16.0, 81.0, 5.0, 323.0, 286.0, 112.0, 15.0, 16.0, 88.0, 87.0, 9.0, 280.0, 165.0, 410.0, 208.0, 325.0, 149.0, 22.0, 72.0, 48.0, 42.0, 126.0, 158.0, 120.0, 11.0, 99.0, 75.0, 32.0, 4.0, 245.0, 374.0, 17.0, 146.0, 222.0, 249.0, 249.0, 17.0, 118.0, 199.0, 94.0, 29.0, 13.0, 85.0, 27.0, 244.0, 380.0, 177.0, 163.0, 62.0, 4.0, 495.0, 97.0, 148.0, 381.0, 56.0, 34.0, 11.0, 43.0, 94.0, 251.0, 67.0, 64.0, 25.0, 32.0]
d = [49.0, 49.0, 39.0, 23.0, 31.0, 48.0, 11.0, 48.0, 26.0, 10.0, 1.0, 20.0, 40.0, 13.0, 29.0, 42.0, 7.0, 18.0, 18.0, 1.0, 30.0, 20.0, 6.0, 38.0, 19.0, 32.0, 2.0, 39.0, 49.0, 3.0, 7.0, 13.0, 50.0, 36.0, 16.0, 45.0, 47.0, 11.0, 4.0, 10.0, 18.0, 11.0, 49.0, 37.0, 28.0, 34.0, 31.0, 48.0, 45.0, 34.0, 28.0, 2.0, 29.0, 28.0, 25.0, 33.0, 44.0, 30.0, 26.0, 15.0, 29.0, 23.0, 23.0, 21.0, 15.0, 22.0, 29.0, 47.0, 29.0, 41.0, 50.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]

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