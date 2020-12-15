edge = [1 2; 1 4; 1 5; 1 9; 1 11; 2 3; 2 4; 2 6; 2 7; 2 8; 2 9; 2 13; 3 6; 3 8; 3 10; 3 13; 3 14; 4 2; 4 7; 4 8; 4 12; 4 14; 5 9; 5 11; 5 13; 5 14; 6 7; 6 11; 6 13; 6 14; 6 15; 7 3; 7 8; 7 11; 8 2; 8 3; 8 7; 8 9; 8 12; 8 14; 8 15; 9 3; 9 5; 9 8; 9 12; 9 15; 10 3; 10 4; 10 5; 10 6; 10 7; 10 9; 10 12; 10 14; 11 5; 11 6; 11 14; 11 15; 12 5; 12 8; 12 9; 12 11; 12 13; 12 14; 13 2; 13 3; 13 5; 13 11; 13 14; 14 3; 14 5; 14 6; 14 8; 14 9]
cL_orig = [1.0, 109.0, 51.0, 300.0, 24.0, 23.0, 78.0, 91.0, 26.0, 35.0, 13.0, 27.0, 52.0, 99.0, 262.0, 123.0, 78.0, 61.0, 36.0, 25.0, 22.0, 22.0, 170.0, 199.0, 327.0, 424.0, 34.0, 146.0, 52.0, 2.0, 155.0, 121.0, 6.0, 96.0, 115.0, 230.0, 13.0, 7.0, 154.0, 263.0, 75.0, 202.0, 96.0, 5.0, 25.0, 8.0, 194.0, 251.0, 85.0, 180.0, 123.0, 45.0, 92.0, 38.0, 293.0, 118.0, 116.0, 198.0, 26.0, 80.0, 36.0, 9.0, 43.0, 75.0, 471.0, 272.0, 215.0, 64.0, 26.0, 451.0, 169.0, 38.0, 42.0, 123.0]
cU_orig = [19.0, 109.0, 51.0, 300.0, 958.0, 23.0, 78.0, 91.0, 26.0, 107.0, 507.0, 27.0, 106.0, 99.0, 364.0, 123.0, 78.0, 61.0, 36.0, 35.0, 22.0, 384.0, 170.0, 199.0, 327.0, 424.0, 34.0, 146.0, 52.0, 6.0, 155.0, 121.0, 80.0, 300.0, 115.0, 230.0, 13.0, 7.0, 154.0, 263.0, 75.0, 202.0, 206.0, 5.0, 163.0, 8.0, 194.0, 251.0, 85.0, 180.0, 123.0, 45.0, 92.0, 118.0, 293.0, 118.0, 144.0, 198.0, 26.0, 92.0, 72.0, 13.0, 43.0, 75.0, 471.0, 272.0, 215.0, 64.0, 26.0, 451.0, 727.0, 38.0, 320.0, 123.0]
d = [12.0, 8.0, 26.0, 23.0, 45.0, 31.0, 2.0, 14.0, 11.0, 25.0, 4.0, 43.0, 48.0, 6.0, 29.0, 5.0, 43.0, 48.0, 42.0, 40.0, 33.0, 18.0, 16.0, 49.0, 13.0, 8.0, 43.0, 31.0, 19.0, 46.0, 28.0, 39.0, 16.0, 41.0, 42.0, 1.0, 36.0, 37.0, 4.0, 11.0, 15.0, 9.0, 28.0, 41.0, 48.0, 14.0, 2.0, 24.0, 29.0, 40.0, 10.0, 4.0, 32.0, 25.0, 31.0, 23.0, 40.0, 36.0, 4.0, 14.0, 40.0, 25.0, 50.0, 38.0, 28.0, 41.0, 28.0, 28.0, 16.0, 11.0, 18.0, 36.0, 16.0, 14.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
