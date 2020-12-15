edge = [1 2; 1 3; 1 4; 1 6; 1 8; 1 9; 1 10; 2 4; 2 6; 2 10; 2 13; 3 2; 3 6; 3 8; 3 15; 4 2; 4 5; 4 7; 4 11; 4 13; 4 15; 5 7; 5 11; 5 14; 5 15; 6 3; 6 5; 6 8; 7 2; 7 4; 7 6; 7 8; 7 10; 7 12; 7 15; 8 3; 8 9; 8 11; 9 2; 9 7; 9 8; 9 14; 10 6; 10 8; 10 9; 10 12; 10 15; 11 3; 11 7; 12 3; 12 9; 12 10; 12 13; 12 15; 13 12; 14 5; 14 6; 14 7; 14 15]
cL_orig = [4.0, 12.0, 141.0, 24.0, 173.0, 105.0, 279.0, 88.0, 78.0, 298.0, 100.0, 11.0, 108.0, 97.0, 371.0, 33.0, 13.0, 8.0, 153.0, 8.0, 143.0, 76.0, 143.0, 87.0, 375.0, 111.0, 9.0, 31.0, 46.0, 47.0, 26.0, 17.0, 60.0, 233.0, 267.0, 51.0, 2.0, 123.0, 127.0, 61.0, 28.0, 90.0, 160.0, 31.0, 0.0, 96.0, 201.0, 395.0, 177.0, 141.0, 41.0, 37.0, 1.0, 103.0, 37.0, 17.0, 6.0, 174.0, 28.0]
cU_orig = [8.0, 12.0, 141.0, 138.0, 173.0, 121.0, 279.0, 88.0, 152.0, 298.0, 124.0, 11.0, 108.0, 97.0, 371.0, 53.0, 13.0, 66.0, 153.0, 450.0, 851.0, 76.0, 377.0, 645.0, 375.0, 111.0, 9.0, 31.0, 46.0, 47.0, 26.0, 31.0, 60.0, 233.0, 267.0, 167.0, 2.0, 123.0, 127.0, 61.0, 28.0, 90.0, 160.0, 45.0, 88.0, 96.0, 201.0, 395.0, 177.0, 141.0, 41.0, 37.0, 15.0, 103.0, 37.0, 287.0, 22.0, 174.0, 34.0]
d = [43.0, 33.0, 4.0, 41.0, 50.0, 46.0, 28.0, 12.0, 2.0, 47.0, 30.0, 36.0, 8.0, 43.0, 27.0, 12.0, 26.0, 48.0, 2.0, 23.0, 22.0, 19.0, 27.0, 44.0, 15.0, 4.0, 20.0, 2.0, 4.0, 39.0, 48.0, 7.0, 22.0, 44.0, 36.0, 2.0, 5.0, 29.0, 8.0, 45.0, 34.0, 44.0, 1.0, 7.0, 28.0, 18.0, 14.0, 47.0, 18.0, 18.0, 1.0, 33.0, 20.0, 34.0, 47.0, 8.0, 34.0, 33.0, 40.0]
Len = length(d)

yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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