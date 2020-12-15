edge = [1 2; 1 8; 1 13; 1 15; 2 4; 2 10; 2 13; 3 6; 3 7; 3 14; 3 15; 4 2; 4 6; 4 8; 4 12; 4 15; 5 2; 5 3; 5 4; 5 15; 6 3; 6 7; 6 8; 6 11; 6 14; 6 15; 7 3; 7 5; 7 6; 7 13; 7 14; 8 3; 8 5; 8 9; 8 12; 8 14; 9 8; 9 13; 10 2; 10 3; 10 4; 10 7; 10 8; 10 12; 11 2; 11 6; 11 8; 11 10; 11 14; 12 2; 12 7; 13 5; 13 7; 13 8; 14 2; 14 6; 14 7]
cL_orig = [10.0, 47.0, 7.0, 274.0, 2.0, 377.0, 338.0, 78.0, 134.0, 531.0, 160.0, 97.0, 31.0, 43.0, 90.0, 145.0, 130.0, 39.0, 7.0, 5.0, 77.0, 6.0, 83.0, 181.0, 280.0, 252.0, 122.0, 83.0, 20.0, 97.0, 58.0, 200.0, 41.0, 41.0, 133.0, 92.0, 31.0, 176.0, 105.0, 350.0, 197.0, 55.0, 8.0, 97.0, 166.0, 151.0, 56.0, 11.0, 116.0, 300.0, 88.0, 2.0, 172.0, 106.0, 19.0, 308.0, 155.0]
cU_orig = [58.0, 179.0, 7.0, 274.0, 106.0, 377.0, 338.0, 78.0, 134.0, 531.0, 160.0, 97.0, 31.0, 43.0, 90.0, 145.0, 130.0, 69.0, 7.0, 21.0, 77.0, 6.0, 103.0, 181.0, 280.0, 252.0, 122.0, 83.0, 20.0, 97.0, 442.0, 200.0, 41.0, 41.0, 143.0, 92.0, 31.0, 216.0, 105.0, 350.0, 197.0, 55.0, 8.0, 97.0, 166.0, 151.0, 56.0, 63.0, 116.0, 300.0, 368.0, 66.0, 172.0, 106.0, 145.0, 308.0, 193.0]
d = [25.0, 11.0, 13.0, 44.0, 7.0, 35.0, 11.0, 35.0, 13.0, 37.0, 15.0, 18.0, 26.0, 49.0, 4.0, 42.0, 24.0, 50.0, 29.0, 17.0, 38.0, 47.0, 42.0, 36.0, 16.0, 24.0, 8.0, 27.0, 25.0, 36.0, 36.0, 7.0, 21.0, 1.0, 23.0, 11.0, 23.0, 48.0, 13.0, 20.0, 25.0, 24.0, 3.0, 29.0, 8.0, 43.0, 28.0, 37.0, 50.0, 5.0, 31.0, 21.0, 4.0, 20.0, 24.0, 24.0, 30.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

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
