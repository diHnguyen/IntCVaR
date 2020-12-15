edge = [1 3; 1 5; 1 6; 1 7; 1 11; 2 7; 2 13; 3 2; 3 6; 3 12; 3 13; 4 7; 4 8; 5 4; 5 8; 5 12; 6 13; 6 15; 7 2; 7 10; 7 11; 7 14; 8 2; 8 3; 8 7; 8 15; 9 3; 9 4; 9 5; 9 8; 9 11; 9 12; 9 15; 10 2; 10 3; 10 6; 10 7; 10 8; 10 9; 10 13; 10 15; 11 2; 11 6; 11 9; 12 3; 12 4; 12 9; 12 11; 12 13; 13 3; 13 4; 14 4; 14 8; 14 9; 14 11; 14 12; 14 13]
cL_orig = [52.0, 7.0, 245.0, 214.0, 149.0, 240.0, 15.0, 12.0, 74.0, 359.0, 243.0, 93.0, 19.0, 33.0, 115.0, 121.0, 217.0, 421.0, 54.0, 128.0, 157.0, 69.0, 168.0, 181.0, 17.0, 279.0, 158.0, 101.0, 117.0, 33.0, 59.0, 7.0, 102.0, 61.0, 68.0, 75.0, 80.0, 1.0, 17.0, 147.0, 99.0, 339.0, 157.0, 33.0, 112.0, 53.0, 31.0, 26.0, 12.0, 69.0, 29.0, 72.0, 271.0, 16.0, 13.0, 99.0, 10.0]
cU_orig = [52.0, 53.0, 245.0, 214.0, 149.0, 240.0, 453.0, 12.0, 74.0, 359.0, 243.0, 147.0, 19.0, 33.0, 115.0, 167.0, 217.0, 421.0, 130.0, 128.0, 157.0, 69.0, 286.0, 181.0, 17.0, 279.0, 158.0, 133.0, 117.0, 33.0, 67.0, 67.0, 150.0, 61.0, 96.0, 75.0, 80.0, 9.0, 39.0, 147.0, 215.0, 339.0, 157.0, 33.0, 112.0, 53.0, 31.0, 26.0, 12.0, 287.0, 29.0, 72.0, 271.0, 80.0, 13.0, 99.0, 10.0]
d = [10.0, 41.0, 9.0, 36.0, 16.0, 48.0, 28.0, 41.0, 4.0, 50.0, 33.0, 19.0, 2.0, 24.0, 46.0, 9.0, 41.0, 45.0, 11.0, 14.0, 37.0, 14.0, 16.0, 5.0, 6.0, 16.0, 33.0, 7.0, 2.0, 50.0, 39.0, 43.0, 46.0, 45.0, 18.0, 11.0, 37.0, 29.0, 2.0, 30.0, 37.0, 42.0, 12.0, 13.0, 35.0, 11.0, 12.0, 46.0, 23.0, 42.0, 21.0, 23.0, 28.0, 26.0, 21.0, 26.0, 32.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
