edge = [1 6; 1 7; 1 8; 1 10; 1 14; 2 3; 2 4; 2 7; 2 11; 2 13; 2 14; 3 2; 3 6; 3 8; 3 9; 3 14; 4 2; 4 7; 4 12; 5 10; 5 11; 5 12; 5 14; 6 3; 6 8; 6 9; 6 12; 7 8; 7 10; 7 15; 8 11; 9 4; 9 5; 9 8; 9 12; 10 3; 10 8; 10 13; 11 4; 11 5; 11 10; 12 4; 12 7; 12 9; 12 11; 12 13; 13 4; 13 6; 13 11; 14 3; 14 4; 14 5; 14 6; 14 7; 14 9; 14 10; 14 12]
cL_orig = [71.0, 1.0, 46.0, 261.0, 523.0, 1.0, 30.0, 50.0, 341.0, 158.0, 120.0, 35.0, 131.0, 72.0, 263.0, 525.0, 38.0, 95.0, 290.0, 116.0, 9.0, 82.0, 399.0, 101.0, 91.0, 115.0, 230.0, 16.0, 70.0, 245.0, 86.0, 233.0, 108.0, 49.0, 79.0, 7.0, 71.0, 85.0, 33.0, 130.0, 45.0, 171.0, 120.0, 48.0, 10.0, 12.0, 281.0, 255.0, 16.0, 336.0, 116.0, 264.0, 371.0, 258.0, 177.0, 144.0, 33.0]
cU_orig = [71.0, 3.0, 46.0, 399.0, 523.0, 1.0, 30.0, 50.0, 341.0, 444.0, 120.0, 35.0, 131.0, 72.0, 263.0, 525.0, 38.0, 95.0, 290.0, 116.0, 9.0, 82.0, 399.0, 101.0, 91.0, 115.0, 230.0, 22.0, 70.0, 245.0, 86.0, 233.0, 254.0, 49.0, 79.0, 11.0, 85.0, 85.0, 615.0, 130.0, 45.0, 171.0, 362.0, 48.0, 10.0, 12.0, 411.0, 255.0, 126.0, 434.0, 116.0, 264.0, 371.0, 258.0, 177.0, 144.0, 33.0]
d = [33.0, 44.0, 33.0, 44.0, 4.0, 4.0, 36.0, 40.0, 17.0, 47.0, 24.0, 36.0, 38.0, 32.0, 5.0, 34.0, 39.0, 34.0, 23.0, 23.0, 29.0, 16.0, 46.0, 13.0, 32.0, 46.0, 16.0, 36.0, 39.0, 27.0, 23.0, 40.0, 44.0, 42.0, 20.0, 8.0, 44.0, 31.0, 8.0, 2.0, 3.0, 34.0, 48.0, 17.0, 17.0, 47.0, 45.0, 9.0, 8.0, 48.0, 25.0, 7.0, 9.0, 22.0, 1.0, 26.0, 38.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
