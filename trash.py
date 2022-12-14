
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






# first layer, first particle
    x1_1, y1_1, vx1_1, vy1_1 = my_rk4(x_is[0,0], y_is[0,0],vx_star, vy_star,star=False)
    # first layer, second particle
    x1_2, y1_2, vx1_2, vy1_2 = my_rk4(x_is[0,1], y_is[0,1],vx_star, vy_star,star=False)
    # and so on...
    x1_3, y1_3, vx1_3, vy1_3 = my_rk4(x_is[0,2], y_is[0,2],vx_star, vy_star,star=False)
    x1_4, y1_4, vx1_4, vy1_4 = my_rk4(x_is[0,3], y_is[0,3],vx_star, vy_star,star=False)
    x1_5, y1_5, vx1_5, vy1_5 = my_rk4(x_is[0,4], y_is[0,4],vx_star, vy_star,star=False)
    x1_6, y1_6, vx1_6, vy1_6 = my_rk4(x_is[0,5], y_is[0,5],vx_star, vy_star,star=False)
    x1_7, y1_7, vx1_7, vy1_7 = my_rk4(x_is[0,6], y_is[0,6],vx_star, vy_star,star=False)
    x1_8, y1_8, vx1_8, vy1_8 = my_rk4(x_is[0,7], y_is[0,7],vx_star, vy_star,star=False)

    # second layer, first particle and so on....
    x2_1, y2_1, vx2_1, vy2_1 = my_rk4(x_is[1,0], y_is[1,0],vx_star, vy_star,star=False)
    x2_2, y2_2, vx2_2, vy2_2 = my_rk4(x_is[1,1], y_is[1,1],vx_star, vy_star,star=False)
    x2_3, y2_3, vx2_3, vy2_3 = my_rk4(x_is[1,2], y_is[1,2],vx_star, vy_star,star=False)
    x2_4, y2_4, vx2_4, vy2_4 = my_rk4(x_is[1,3], y_is[1,3],vx_star, vy_star,star=False)
    x2_5, y2_5, vx2_5, vy2_5 = my_rk4(x_is[1,4], y_is[1,4],vx_star, vy_star,star=False)
    x2_6, y2_6, vx2_6, vy2_6 = my_rk4(x_is[1,5], y_is[1,5],vx_star, vy_star,star=False)
    x2_7, y2_7, vx2_7, vy2_7 = my_rk4(x_is[1,6], y_is[1,6],vx_star, vy_star,star=False)
    x2_8, y2_8, vx2_8, vy2_8 = my_rk4(x_is[1,7], y_is[1,7],vx_star, vy_star,star=False)


    # third layer, all particles
    x3_1, y3_1, vx3_1, vy3_1 = my_rk4(x_is[2,0], y_is[2,0],vx_star, vy_star,star=False)
    x3_2, y3_2, vx3_2, vy3_2 = my_rk4(x_is[2,1], y_is[2,1],vx_star, vy_star,star=False)
    x3_3, y3_3, vx3_3, vy3_3 = my_rk4(x_is[2,2], y_is[2,2],vx_star, vy_star,star=False)
    x3_4, y3_4, vx3_4, vy3_4 = my_rk4(x_is[2,3], y_is[2,3],vx_star, vy_star,star=False)
    x3_5, y3_5, vx3_5, vy3_5 = my_rk4(x_is[2,4], y_is[2,4],vx_star, vy_star,star=False)
    x3_6, y3_6, vx3_6, vy3_6 = my_rk4(x_is[2,5], y_is[2,5],vx_star, vy_star,star=False)
    x3_7, y3_7, vx3_7, vy3_7 = my_rk4(x_is[2,6], y_is[2,6],vx_star, vy_star,star=False)
    x3_8, y3_8, vx3_8, vy3_8 = my_rk4(x_is[2,7], y_is[2,7],vx_star, vy_star,star=False)


    # fourth layer, all particles
    x4_1, y4_1, vx4_1, vy4_1 = my_rk4(x_is[3,0], y_is[3,0],vx_star, vy_star,star=False)
    x4_2, y4_2, vx4_2, vy4_2 = my_rk4(x_is[3,1], y_is[3,1],vx_star, vy_star,star=False)
    x4_3, y4_3, vx4_3, vy4_3 = my_rk4(x_is[3,2], y_is[3,2],vx_star, vy_star,star=False)
    x4_4, y4_4, vx4_4, vy4_4 = my_rk4(x_is[3,3], y_is[3,3],vx_star, vy_star,star=False)
    x4_5, y4_5, vx4_5, vy4_5 = my_rk4(x_is[3,4], y_is[3,4],vx_star, vy_star,star=False)
    x4_6, y4_6, vx4_6, vy4_6 = my_rk4(x_is[3,5], y_is[3,5],vx_star, vy_star,star=False)
    x4_7, y4_7, vx4_7, vy4_7 = my_rk4(x_is[3,6], y_is[3,6],vx_star, vy_star,star=False)
    x4_8, y4_8, vx4_8, vy4_8 = my_rk4(x_is[3,7], y_is[3,7],vx_star, vy_star,star=False)


    # fifth layer, all particles
    x5_1, y5_1, vx5_1, vy5_1 = my_rk4(x_is[4,0], y_is[4,0],vx_star, vy_star,star=False)
    x5_2, y5_2, vx5_2, vy5_2 = my_rk4(x_is[4,1], y_is[4,1],vx_star, vy_star,star=False)
    x5_3, y5_3, vx5_3, vy5_3 = my_rk4(x_is[4,2], y_is[4,2],vx_star, vy_star,star=False)
    x5_4, y5_4, vx5_4, vy5_4 = my_rk4(x_is[4,3], y_is[4,3],vx_star, vy_star,star=False)
    x5_5, y5_5, vx5_5, vy5_5 = my_rk4(x_is[4,4], y_is[4,4],vx_star, vy_star,star=False)
    x5_6, y5_6, vx5_6, vy5_6 = my_rk4(x_is[4,5], y_is[4,5],vx_star, vy_star,star=False)
    x5_7, y5_7, vx5_7, vy5_7 = my_rk4(x_is[4,6], y_is[4,6],vx_star, vy_star,star=False)
    x5_8, y5_8, vx5_8, vy5_8 = my_rk4(x_is[4,7], y_is[4,7],vx_star, vy_star,star=False)


