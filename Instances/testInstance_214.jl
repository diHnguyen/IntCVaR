edge = [1 2; 1 5; 1 7; 1 8; 1 12; 2 3; 2 5; 2 8; 2 13; 3 4; 3 10; 3 11; 3 12; 3 13; 3 15; 4 3; 4 5; 4 10; 4 13; 6 2; 6 7; 6 11; 6 12; 6 13; 6 14; 7 3; 7 5; 7 8; 7 9; 7 10; 7 11; 7 12; 8 2; 8 4; 8 5; 8 6; 8 12; 8 14; 9 2; 9 5; 9 7; 9 10; 9 12; 9 13; 10 2; 10 7; 10 9; 10 11; 10 14; 10 15; 11 2; 11 4; 11 5; 11 9; 11 10; 11 14; 12 2; 12 4; 12 6; 12 9; 13 2; 13 3; 13 5; 13 6; 13 15; 14 2; 14 4; 14 10]
cL_orig = [16.0, 173.0, 287.0, 163.0, 178.0, 3.0, 11.0, 232.0, 232.0, 8.0, 224.0, 53.0, 82.0, 344.0, 422.0, 7.0, 14.0, 53.0, 139.0, 168.0, 9.0, 220.0, 158.0, 325.0, 9.0, 197.0, 100.0, 32.0, 20.0, 41.0, 42.0, 102.0, 21.0, 61.0, 34.0, 0.0, 88.0, 126.0, 136.0, 35.0, 20.0, 37.0, 15.0, 128.0, 246.0, 6.0, 23.0, 12.0, 3.0, 220.0, 45.0, 257.0, 284.0, 25.0, 1.0, 127.0, 127.0, 212.0, 169.0, 56.0, 393.0, 65.0, 235.0, 104.0, 38.0, 189.0, 386.0, 53.0]
cU_orig = [78.0, 173.0, 287.0, 163.0, 178.0, 3.0, 11.0, 232.0, 232.0, 16.0, 224.0, 53.0, 82.0, 344.0, 770.0, 29.0, 28.0, 53.0, 139.0, 168.0, 21.0, 220.0, 222.0, 325.0, 237.0, 197.0, 100.0, 32.0, 20.0, 41.0, 42.0, 102.0, 21.0, 61.0, 34.0, 4.0, 98.0, 126.0, 136.0, 35.0, 96.0, 47.0, 221.0, 128.0, 246.0, 112.0, 23.0, 46.0, 3.0, 220.0, 45.0, 257.0, 284.0, 25.0, 17.0, 127.0, 143.0, 212.0, 169.0, 56.0, 393.0, 65.0, 235.0, 104.0, 38.0, 189.0, 386.0, 53.0]
d = [10.0, 40.0, 40.0, 50.0, 32.0, 27.0, 25.0, 20.0, 47.0, 49.0, 28.0, 38.0, 43.0, 25.0, 34.0, 8.0, 20.0, 6.0, 33.0, 29.0, 3.0, 15.0, 1.0, 24.0, 34.0, 12.0, 6.0, 29.0, 21.0, 18.0, 44.0, 1.0, 1.0, 47.0, 4.0, 24.0, 18.0, 37.0, 11.0, 15.0, 3.0, 1.0, 5.0, 37.0, 12.0, 44.0, 48.0, 26.0, 40.0, 30.0, 20.0, 37.0, 40.0, 23.0, 20.0, 46.0, 41.0, 22.0, 50.0, 23.0, 3.0, 29.0, 16.0, 21.0, 8.0, 33.0, 26.0, 13.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

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