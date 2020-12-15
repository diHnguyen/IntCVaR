edge = [1 8; 1 9; 1 10; 1 11; 1 13; 1 14; 2 5; 2 6; 2 8; 3 4; 3 6; 3 7; 3 8; 3 11; 3 12; 4 3; 4 5; 4 6; 4 8; 4 11; 5 7; 5 8; 5 10; 5 11; 5 14; 6 2; 6 7; 6 8; 6 9; 6 10; 6 11; 6 12; 6 13; 7 2; 7 3; 7 10; 7 11; 7 14; 8 2; 8 5; 8 10; 8 11; 8 15; 9 3; 9 4; 9 5; 9 13; 9 15; 10 5; 10 6; 10 9; 10 12; 10 13; 10 14; 11 13; 12 2; 12 11; 12 15; 13 6; 13 10; 13 12; 14 2; 14 4; 14 5; 14 7; 14 12; 14 15]
cL_orig = [125.0, 292.0, 202.0, 13.0, 271.0, 22.0, 86.0, 142.0, 151.0, 34.0, 9.0, 159.0, 15.0, 268.0, 188.0, 21.0, 24.0, 22.0, 189.0, 138.0, 92.0, 91.0, 118.0, 205.0, 61.0, 176.0, 17.0, 92.0, 31.0, 81.0, 1.0, 204.0, 273.0, 87.0, 197.0, 71.0, 170.0, 152.0, 171.0, 111.0, 13.0, 105.0, 51.0, 28.0, 12.0, 54.0, 77.0, 236.0, 125.0, 179.0, 24.0, 96.0, 101.0, 106.0, 42.0, 240.0, 43.0, 79.0, 8.0, 67.0, 1.0, 0.0, 265.0, 138.0, 336.0, 11.0, 47.0]
cU_orig = [125.0, 292.0, 202.0, 427.0, 395.0, 350.0, 86.0, 142.0, 151.0, 34.0, 9.0, 159.0, 171.0, 268.0, 188.0, 21.0, 44.0, 22.0, 189.0, 138.0, 92.0, 111.0, 118.0, 205.0, 61.0, 176.0, 17.0, 92.0, 31.0, 233.0, 53.0, 204.0, 273.0, 87.0, 197.0, 71.0, 170.0, 152.0, 253.0, 111.0, 13.0, 105.0, 51.0, 28.0, 80.0, 252.0, 77.0, 306.0, 125.0, 179.0, 48.0, 96.0, 101.0, 106.0, 42.0, 240.0, 53.0, 221.0, 8.0, 67.0, 1.0, 4.0, 265.0, 224.0, 336.0, 37.0, 47.0]
d = [26.0, 4.0, 44.0, 34.0, 2.0, 26.0, 19.0, 45.0, 21.0, 23.0, 48.0, 11.0, 29.0, 8.0, 18.0, 2.0, 28.0, 30.0, 31.0, 4.0, 4.0, 25.0, 50.0, 30.0, 17.0, 44.0, 22.0, 46.0, 21.0, 23.0, 9.0, 48.0, 41.0, 47.0, 46.0, 9.0, 4.0, 5.0, 27.0, 28.0, 44.0, 17.0, 19.0, 4.0, 30.0, 11.0, 45.0, 32.0, 43.0, 16.0, 9.0, 24.0, 3.0, 32.0, 49.0, 6.0, 22.0, 37.0, 17.0, 13.0, 42.0, 23.0, 47.0, 12.0, 13.0, 33.0, 32.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
