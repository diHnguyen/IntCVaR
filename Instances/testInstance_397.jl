edge = [1 5; 1 6; 1 10; 1 12; 1 13; 1 15; 2 5; 2 7; 2 9; 2 10; 2 12; 2 15; 3 4; 3 5; 3 7; 3 8; 3 11; 3 13; 4 6; 4 7; 4 9; 4 10; 4 12; 4 15; 5 2; 5 4; 5 7; 5 8; 5 10; 5 12; 6 7; 6 8; 6 9; 6 10; 6 13; 6 14; 6 15; 7 4; 7 6; 7 9; 7 10; 8 7; 8 15; 9 8; 9 11; 9 12; 9 13; 10 5; 10 7; 10 13; 10 14; 11 10; 12 3; 13 4; 13 6; 13 11; 13 14; 14 2; 14 5; 14 6; 14 10; 14 12]
cL_orig = [123.0, 222.0, 285.0, 497.0, 193.0, 358.0, 144.0, 215.0, 201.0, 273.0, 274.0, 298.0, 46.0, 97.0, 109.0, 128.0, 78.0, 268.0, 2.0, 14.0, 87.0, 230.0, 0.0, 143.0, 83.0, 1.0, 54.0, 34.0, 220.0, 221.0, 22.0, 45.0, 50.0, 69.0, 32.0, 81.0, 202.0, 78.0, 50.0, 59.0, 114.0, 49.0, 38.0, 2.0, 9.0, 140.0, 25.0, 210.0, 95.0, 63.0, 43.0, 7.0, 174.0, 131.0, 74.0, 39.0, 23.0, 44.0, 350.0, 203.0, 69.0, 18.0]
cU_orig = [273.0, 222.0, 285.0, 497.0, 193.0, 358.0, 144.0, 215.0, 201.0, 273.0, 274.0, 450.0, 46.0, 97.0, 109.0, 128.0, 78.0, 268.0, 8.0, 36.0, 277.0, 230.0, 12.0, 143.0, 83.0, 3.0, 138.0, 34.0, 238.0, 221.0, 22.0, 45.0, 50.0, 69.0, 428.0, 81.0, 202.0, 78.0, 50.0, 59.0, 114.0, 49.0, 38.0, 2.0, 9.0, 140.0, 77.0, 210.0, 95.0, 63.0, 43.0, 7.0, 174.0, 131.0, 74.0, 39.0, 23.0, 44.0, 350.0, 293.0, 69.0, 126.0]
d = [49.0, 49.0, 29.0, 39.0, 11.0, 17.0, 21.0, 39.0, 13.0, 2.0, 31.0, 49.0, 24.0, 32.0, 33.0, 21.0, 50.0, 46.0, 33.0, 36.0, 8.0, 24.0, 10.0, 34.0, 17.0, 3.0, 4.0, 45.0, 35.0, 20.0, 45.0, 34.0, 40.0, 47.0, 17.0, 49.0, 8.0, 31.0, 2.0, 27.0, 44.0, 13.0, 6.0, 11.0, 20.0, 35.0, 50.0, 36.0, 4.0, 4.0, 11.0, 38.0, 29.0, 19.0, 22.0, 36.0, 46.0, 40.0, 9.0, 16.0, 38.0, 19.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
