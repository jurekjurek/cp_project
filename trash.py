
# grav layer sorting by radii

# def grav_layer(r_star, m_star, split = 5):
#     '''
#     Calculate the different radii for the different layers
#     assuming constant density
#     '''
#     m = m_star 
#     r1 = r_star / 5
#     r2 = 2*r1
#     r3 = 3*r1
#     r4 = 4*r1
#     r5 = r_star

#     m1 = m * ((1/5)**3)
#     m2 = m * ((2/5)**3-(1/5)**3)
#     m3 = m * ((3/5)**3-(2/5)**3)
#     m4 = m * ((4/5)**3-(3/5)**3)
#     m5 = m * ((5/5)**3-(4/5)**3)

#     # gravitational pull on i-th layer by star
#     # F5 = G*1/5*m*(m-1/5*m)/((r4)**2)
#     # F4 = G*1/5*m*(m-2/5*m)/((r3)**2)
#     # F3 = G*1/5*m*(m-3/5*m)/((r2)**2)
#     # F2 = G*1/5*m*(m-4/5*m)/((r1)**2)
#     # F_whole_star = G*m**2/(r_star/2)**2

#     F5 = G*m5*(m-m5)/((r4)**2)
#     F4 = G*m4*(m-m4)/((r3)**2)
#     F3 = G*m3*(m-m3)/((r2)**2)
#     F2 = G*m2*(m-m2)/((r1)**2)

#     A2 = G*M*m
#     A3 = G*M*m
#     A4 = G*M*m
#     A5 = G*M*m
#     B2 = 3*(G*M/c)**2*m
#     B3 = 3*(G*M/c)**2*m
#     B4 = 3*(G*M/c)**2*m
#     B5 = 3*(G*M/c)**2*m

#     print(F2, F3, F4, F5) 

#     print(A, B, F3)
#     r2_thr = r_help(A2,B2,F2).real
#     r3_thr = r_help(A3,B3,F3).real
#     r4_thr = r_help(A4,B4,F4).real
#     r5_thr = r_help(A5,B5,F5).real

#     # r_whole_star = r_help(A,B).real
#     # eq to solve: -A/r^2 - B/r^3 = -F_i for r

#     return r2_thr-r2, r3_thr-r3, r4_thr-r4, r5_thr-r5





# def test():
#     '''
#     function that calculates, when the F_BH on the i-th layer of the star is bigger than the gravitational pull of the star itself
#     '''
#     r = np.sqrt(x**2+y**2)
#     F_bh = -A / r**2 - 3*(G*M/c)**2 / r**3
#     F2, F3, F4, F5 = grav_layer(r_star, m)
#     for i in range(len(t)):
#         if F_bh <= F2:
#             print('Distance for first layer to drop: ', r)
#         if F_bh <= F3:
#             print('Distance for second layer to drop: ', r)
#         if F_bh <= F4:
#             print('Distance for third layer to drop: ', r)
#         if F_bh <= F5:
#             print('Distance for fourth layer to drop: ', r)

# test()