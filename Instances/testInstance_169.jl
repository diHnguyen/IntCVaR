edge = [1 3; 1 7; 1 9; 1 11; 2 4; 2 6; 2 7; 2 9; 2 11; 3 5; 3 6; 3 10; 3 11; 4 10; 4 11; 5 4; 5 11; 5 12; 6 3; 6 8; 6 9; 6 11; 7 5; 7 6; 7 9; 8 2; 8 4; 8 6; 8 7; 8 12; 9 2; 9 3; 9 11; 9 12; 9 13; 9 14; 10 3; 10 4; 10 9; 10 12; 10 15; 11 3; 11 4; 11 8; 12 10; 12 11; 12 13; 12 15; 13 3; 13 11; 13 12; 13 15; 14 7; 14 10; 14 11; 14 15]
cL_orig = [29.0, 29.0, 29.0, 178.0, 60.0, 78.0, 69.0, 311.0, 195.0, 60.0, 61.0, 319.0, 135.0, 31.0, 152.0, 12.0, 256.0, 67.0, 116.0, 46.0, 19.0, 85.0, 9.0, 27.0, 72.0, 38.0, 126.0, 70.0, 32.0, 33.0, 19.0, 27.0, 11.0, 134.0, 35.0, 221.0, 4.0, 278.0, 3.0, 78.0, 188.0, 328.0, 301.0, 101.0, 96.0, 48.0, 1.0, 78.0, 206.0, 89.0, 10.0, 63.0, 313.0, 105.0, 14.0, 37.0]
cU_orig = [29.0, 29.0, 29.0, 178.0, 60.0, 78.0, 69.0, 311.0, 471.0, 66.0, 61.0, 319.0, 135.0, 275.0, 152.0, 12.0, 256.0, 297.0, 116.0, 46.0, 19.0, 241.0, 9.0, 27.0, 114.0, 38.0, 126.0, 70.0, 32.0, 35.0, 367.0, 111.0, 11.0, 134.0, 35.0, 221.0, 4.0, 278.0, 3.0, 78.0, 188.0, 328.0, 327.0, 101.0, 96.0, 48.0, 1.0, 78.0, 206.0, 89.0, 10.0, 63.0, 313.0, 217.0, 14.0, 37.0]
d = [17.0, 4.0, 33.0, 5.0, 2.0, 46.0, 45.0, 42.0, 4.0, 7.0, 34.0, 5.0, 32.0, 29.0, 40.0, 16.0, 36.0, 4.0, 2.0, 13.0, 10.0, 23.0, 27.0, 43.0, 21.0, 5.0, 15.0, 25.0, 23.0, 49.0, 7.0, 13.0, 18.0, 5.0, 23.0, 38.0, 40.0, 15.0, 32.0, 27.0, 5.0, 41.0, 32.0, 4.0, 3.0, 30.0, 32.0, 11.0, 46.0, 24.0, 41.0, 5.0, 47.0, 46.0, 33.0, 28.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]

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
