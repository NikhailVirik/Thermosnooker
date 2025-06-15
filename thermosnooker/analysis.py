"""Analysis Module."""
# pylint: disable=consider-using-enumerate
import matplotlib.pyplot as plt
from scipy import stats as stats
from scipy.optimize import curve_fit
from thermosnooker.balls import Container
from thermosnooker.balls import Ball
from thermosnooker.simulations import *
from thermosnooker.physics import maxwell
def task9():
    """
    Task 9.

    In this function, you should test your animation. To do this, create a container
    and ball as directed in the project brief. Create a SingleBallSimulation object from these
    and try running your animation. Ensure that this function returns the balls final position and
    velocity.

    Returns:
        tuple[NDArray[np.float64], NDArray[np.float64]]: The balls final position and velocity
    """

    c = Container(radius=10.)
    b = Ball(pos=[-5, 0], vel=[1, 0], radius=1., mass=1.)
    sbs = SingleBallSimulation(container=c, ball=b)
    sbs.run(20, animate = True, pause_time = 0.5)
    return (b.pos(),b.vel())

def task10():
    """
    Task 10.
    Returns:
        [Figure]: animation of simulation
    """
    mbs = MultiBallSimulation(c_radius=10, b_radius=1, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=2)
    mbs.run(10000, animate = True, pause_time=0.03)


def task11():
    """
    Task 11.

    In this function we shall be quantitatively checking that the balls aren't escaping or sticking.
    To do this, create the two histograms as directed in the project script. Ensure that these two
    histogram figures are returned.

    Returns:
        tuple[Figure, Firgure]: The histograms (distance from centre, inter-ball spacing).
    """
    mbs=MultiBallSimulation(c_radius=10,b_radius=1,b_speed=10,b_mass=1,rmax=8,nrings=3,multi=6)
    mbs.run(500, animate = False, pause_time=0)
    c_dist=[]
    for i in range(0,len(mbs.balls())):
        tmpt= mbs.balls()[i].pos()
        mag= np.sqrt((tmpt[0]*tmpt[0])+(tmpt[1]*tmpt[1]))
        c_dist.append(mag)
    b_dist=[]
    for n in range(0,len(mbs.balls())):
        for m in range(0,len(mbs.balls())):
            if n>m:
                rel_d=mbs.balls()[n].pos()-mbs.balls()[m].pos()
                mag=np.sqrt((rel_d[0]*rel_d[0])+(rel_d[1]*rel_d[1]))
                b_dist.append(mag)
            else:
                continue
    #x1=np.linspace(min(c_dist), max(c_dist),100)
    #x2=np.linspace(min(b_dist), max(b_dist),100)
    a=plt.hist(c_dist)
    plt.title('Avg distance from centre')
    plt.xlabel('Distance /m')
    plt.ylabel('No. of Balls')

    plt.show()
    b=plt.hist(b_dist)
    plt.title('Avg inter-ball dist')
    plt.xlabel('Distance /m')
    plt.ylabel('No. of Balls')
    plt.show()
    return [a,b]


def task12():
    """
    Task 12.
    Tests for floating point KE Mom-x mom-y adnd total pressure evolution.]

    Returns:
        tuple[Figure, Figure, Figure, Figure]: KE, momentum_x, momentum_y ratios, pressure.
    """
    mbs = MultiBallSimulation(c_radius=10,b_radius=1,b_speed=10,b_mass=1,rmax=8,nrings=3,multi=6)
    times=[]
    eks=[]
    p=[]
    pxs=[]
    pys=[]
    press=[]
    for n in range(0,1000):
        print(n, end='\r')
        times.append(mbs.time())
        eks.append(mbs.kinetic_energy())
        p.append(mbs.momentum())
        press.append(mbs.pressure())
        mbs.next_collision()
    rpx=[]
    rpy=[]

    for i in range(0,len(p)):
        pxs.append(p[i][0])
        pys.append(p[i][1])
        rpx.append(pxs[i]/pxs[0])
        rpy.append(pys[i]/pys[0])
    rel_e=[]
    for i in range(0,len(eks)):
        rel_e.append(eks[i]/eks[0])
    print('Ks', rel_e)

    a=plt.plot(times,rel_e)
    plt.ylim(1-1E-14,1+1E-14)
    plt.title('Rel KE vs t')
    plt.xlabel('time /s')
    plt.ylabel('KE/KEi')
    plt.show()
    b=plt.plot(times,rpx)
    plt.ylim(1-1E-14,1+1E-14)
    plt.title('Rel Mom-x vs time')
    plt.xlabel('time /s')
    plt.ylabel(r'$p_{x}/p_{xi}')
    plt.show()
    c=plt.plot(times,rpy)
    plt.ylim(1-1E-14,1+1E-14)
    plt.title('Rel Mom-y vs time')
    plt.xlabel('time /s')
    plt.ylabel(r'$p_{y}/p_{yi}')
    plt.show()
    d=plt.plot(times,press)
    plt.title('Pressure vs time')
    plt.xlabel('time /s')
    plt.ylabel('pressure /Pa')
    plt.show()

    return [a,b,c,d]

def task13():
    """
    Task 13.

    In this function we investigate how well our simulation reproduces the distributions of the IGL.
    Create the 3 figures directed by the project script, namely:
    1) PT plot
    2) PV plot
    3) PN plot
    Ensure that this function returns the three matplotlib figures.

    Returns:
        tuple[Figure, Figure, Figure]: The 3 requested figures: (PT, PV, PN)
    """
    vs=np.linspace(0.1,300, 10)
    k=1.380649E-23
    press__0_1=[]
    test1=[]
    test2=[]
    test3=[]

    temps_0_1=[]
    for n in range(0,len(vs)):
        mbs = MultiBallSimulation(c_radius=10, b_radius=0.1, b_speed=vs[n], b_mass=1,rmax=8, nrings=3, multi=6)
        mbs.run(500, animate=False, pause_time=0)
        press__0_1.append(mbs.pressure())
        vol_dt=float(mbs.container().volume())
        n_dt= int(len(mbs.balls()))
        temps_0_1.append(mbs.t_equipartition())
        Container.dp_tot1 = 0
    press__0_5=[]
    pp_dt_0_5=[]
    test1=[]
    test2=[]
    test3=[]

    temps_0_5=[]
    for n in range(0,len(vs)):
        mbs = MultiBallSimulation(c_radius=10, b_radius=0.5,
                                  b_speed=vs[n], b_mass=1,rmax=8, nrings=3, multi=6)
        mbs.run(500, animate=False, pause_time=0)
        press__0_5.append(mbs.pressure())
        temps_0_5.append(mbs.t_equipartition())
        test1.append(mbs.container().dp_tot())
        test2.append(mbs.time())
        Container.dp_tot1 = 0
    def p_vdw(temp, b):
        k=1.380649E-23
        v=vol_dt
        n=n_dt
        temp=np.array(temp)
        return ((n*k*temp)/(v-n*b))
    err=[]
    p_0_5b_dt, p_0_5b_dt_cov = curve_fit(p_vdw, temps_0_5, press__0_5, p0=[5])
    for i in range(0,len(temps_0_5)):
        predict=(n_dt*k*temps_0_5[i])/vol_dt
        pp_dt_0_5.append(predict)
    err.append(p_0_5b_dt_cov)


    press_1=[]

    press1_1=[]

    for n in range(0,len(vs)):
        mbs = MultiBallSimulation(c_radius=10, b_radius=1, b_speed=vs[n], b_mass=1,rmax=8, nrings=3, multi=6)
        for _ in range(0,200):
            press1_1.append(mbs.pressure())
            mbs.next_collision()
        mbs.run(500, animate=False, pause_time=0)
        test3.append(mbs.container().dp_tot())
        press_1.append(mbs.pressure())
        Container.dp_tot1 = 0
    print(test1,test2,test3)

    rads=np.linspace(10,20,10)
    vols=[]
    press_dv_0_1=[]
    for n in range(0,len(rads)):
        mbs = MultiBallSimulation(c_radius=rads[n], b_radius=0.1, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=6)
        mbs.run(500,animate=False, pause_time=0)
        press_dv_0_1.append(mbs.pressure())
        vols.append(mbs.container().volume())
        temp=mbs.t_equipartition()
        Container.dp_tot1=0
        n_dt= int(len(mbs.balls()))
    press_dv_0_5=[]
    press_pred_dv=[]
    test11=[]
    test21=[]
    test31=[]
    test41=[]
    for n in range(0,len(rads)):
        mbs = MultiBallSimulation(c_radius=rads[n], b_radius=0.5, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=6)
        mbs.run(500,animate=False, pause_time=0)
        press_dv_0_5.append(mbs.pressure())
        test11.append(mbs.container().dp_tot())
        test21.append(mbs.time())
        Container.dp_tot1=0
    for i in range(0,len(vols)):
        predict=(n_dt/vols[i])*k*temp
        press_pred_dv.append(predict)
    press_dv_1=[]
    for n in range(0,len(rads)):
        mbs = MultiBallSimulation(c_radius=rads[n], b_radius=1, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=6)
        mbs.run(500,animate=False, pause_time=0)
        press_dv_1.append(mbs.pressure())
        Container.dp_tot1=0
        test31.append(mbs.container().dp_tot())
        test41.append(mbs.time())
    m2=np.linspace(1,7,7)
    multis=[]
    for l in range(0,len(m2)):
        multis.append(int(m2[l]))

    ns=[]
    press_dn_0_1=[]
    for n in range(0,len(multis)):
        mbs = MultiBallSimulation(c_radius=10, b_radius=0.1, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=multis[n])
        ns.append(len(mbs.balls()))
        mbs.run(500, animate=False, pause_time=0)
        press_dn_0_1.append(mbs.pressure())
        voln=mbs.container().volume()
        #E=0.5*mbs.balls()[0].mass()*10*10
        temp=mbs.t_equipartition()
        Container.dp_tot1 = 0
    press_dn_0_5=[]
    press_pred_dn=[]
    for n in range(0,len(multis)):
        mbs = MultiBallSimulation(c_radius=10, b_radius=0.5, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=multis[n])
        mbs.run(500, animate=False, pause_time=0)
        press_dn_0_5.append(mbs.pressure())
        #E=0.5*mbs.balls()[0].mass()*10*10
        Container.dp_tot1 = 0

    for i in range(0,len(ns)):
        predict=(ns[i]/voln)*k*temp
        press_pred_dn.append(predict)
    press_dn_1=[]
    for n in range(0,len(multis)):
        mbs = MultiBallSimulation(c_radius=10, b_radius=1, b_speed=10, b_mass=1,rmax=8, nrings=3, multi=multis[n])
        mbs.run(500, animate=False, pause_time=0)
        press_dn_1.append(mbs.pressure())
        test1.append(mbs.container().dp_tot())
        test2.append(mbs.time())
        Container.dp_tot1 = 0
    print('Len', len(ns), len(press_dn_0_5), len(press_dn_1))

    err_dt01=[]
    fit_dt_0_1,cov_dt_0_1=np.polyfit(temps_0_1,press__0_1,1,cov='unscaled')
    err_dt01.append(cov_dt_0_1)
    func_dt_0_1=np.poly1d(fit_dt_0_1)
    a=plt.scatter(temps_0_1,press__0_1, label=f'r=0.1m. Nk/V={fit_dt_0_1[0]}', color='orange')
    plt.plot(temps_0_1, func_dt_0_1(temps_0_1), color='purple')
    err_dt05=[]
    fit_dt_0_5,cov_dt_0_5=np.polyfit(temps_0_5,press__0_5,1,cov='unscaled')
    err_dt05.append(cov_dt_0_5)
    func_dt_0_5=np.poly1d(fit_dt_0_5)
    plt.scatter(temps_0_5,press__0_5, label=f'r=0.5m. Nk/V={fit_dt_0_5[0]}', color='blue')
    plt.plot(temps_0_5, func_dt_0_5(temps_0_5), color='purple')
    data_fit_0_5b_dt = p_vdw(temps_0_5, *p_0_5b_dt)
    plt.plot(temps_0_5, data_fit_0_5b_dt, color='pink')
    print('0.5 m dT fit params', p_0_5b_dt )
    err_dt_i=[]
    fit_dt_i,cov_dt_i=np.polyfit(temps_0_5, pp_dt_0_5,1,cov='unscaled')
    err_dt_i.append(cov_dt_i)
    func_dt_i=np.poly1d(fit_dt_i)
    plt.scatter(temps_0_5, pp_dt_0_5, label=f'Ideal. Nk/V={fit_dt_i[0]}', color= 'green')
    plt.plot(temps_0_5, func_dt_i(temps_0_5), color='purple')

    err_dt_1=[]
    fit_dt_1,cov_dt_1=np.polyfit(temps_0_5, press_1,1,cov='unscaled')
    err_dt_1.append(cov_dt_1)
    func_dt_1=np.poly1d(fit_dt_1)
    plt.scatter(temps_0_5, press_1, label=f'r=1m. Nk/V={fit_dt_1[0]}', color='red')
    plt.plot(temps_0_5,func_dt_1(temps_0_5), color='purple')
    plt.title('Pressure vs Temp')
    plt.xlabel('Temp/K')
    plt.ylabel('Pressure/ Pa')
    plt.legend()
    plt.show()

###################################################################################################
    def func_in(x,k):
        return k/x
    x_fit=np.linspace(min(vols),max(vols),100)

    fit_dv_0_1,cov_dv_0_1=curve_fit(func_in,vols, press_dv_0_1)
    k_opt01=fit_dv_0_1[0]
    b=plt.scatter(vols,press_dv_0_1, label=f'r=0.1m. NkT={fit_dv_0_1[0]}', color= 'orange')
    y01=func_in(x_fit, k_opt01)
    plt.plot(x_fit, y01, color="purple")

    fit_dv_0_5,cov_dv_0_5=curve_fit(func_in,vols, press_dv_0_5)
    k_opt05=fit_dv_0_5[0]
    plt.scatter(vols,press_dv_0_5, label=f'r=0.5m. NkT={fit_dv_0_5[0]}', color= 'blue')
    y05=func_in(x_fit, k_opt05)
    plt.plot(x_fit, y05, color="purple")

    fit_dv_i,cov_dv_i=curve_fit(func_in,vols,press_pred_dv)
    k_opti=fit_dv_i[0]
    plt.scatter(vols,press_pred_dv, label= f'Ideal. NkT={fit_dv_i[0]}', color= 'green')
    yi=func_in(x_fit, k_opti)
    plt.plot(x_fit,yi, color='purple')

    fit_dv_1,cov_dv_1=curve_fit(func_in,vols,press_dv_1)
    k_opt1=fit_dv_1[0]
    plt.scatter(vols,press_dv_1, label=f'r=1m. NkT={fit_dv_1[0]}', color='red')
    y1=func_in(x_fit, k_opt1)
    plt.plot(x_fit, y1, color='purple')
    plt.title('Pressure vs Vol')
    plt.xlabel('Vol/ m^3')
    plt.ylabel('Pressure/ Pa')
    plt.legend()
    plt.show()

#############################################################################################
    fit_dn_0_1, cov_dn_0_1=np.polyfit(ns,press_dn_0_1,1,cov='unscaled')
    func_dn_0_1 = np.poly1d(fit_dn_0_1)
    c=plt.scatter(ns,press_dn_0_1, label=f'r=0.1m. kT/V={fit_dn_0_1[0]}', color='orange')
    plt.plot(ns, func_dn_0_1(ns), 'purple')

    fit_dn_0_5, cov_dn_0_5=np.polyfit(ns,press_dn_0_5,1,cov='unscaled')
    func_dn_0_5 = np.poly1d(fit_dn_0_5)
    plt.scatter(ns,press_dn_0_5, label=f'r=0.5m. kT/V={fit_dn_0_5[0]}', color='blue')
    plt.plot(ns, func_dn_0_5(ns), 'purple')

    fit_dn_i, cov_dn_i=np.polyfit(ns,press_pred_dn,1,cov='unscaled')
    func_dn_i=np.poly1d(fit_dn_i)
    plt.scatter(ns,press_pred_dn, label=f'Ideal. kT/V={fit_dn_i[0]}', color='green')
    plt.plot(ns,func_dn_i(ns), color='purple')

    fit_dn_1, cov_dn_1=np.polyfit(ns, press_dn_1,1,cov='unscaled')
    fit_dn_2, cov_dn_2=np.polyfit(ns,press_dn_1,2,cov='unscaled')
    func_dn_1=np.poly1d(fit_dn_1)
    func_dn_2=np.poly1d(fit_dn_2)

    plt.scatter(ns, press_dn_1, label='r=1m. 0.0126x^2 - 0.0196x + 1.059', color= 'red')
    plt.plot(ns,func_dn_1(ns), color="purple")
    plt.plot(ns, func_dn_2(ns), color= "blue")
    plt.title('Pressure vs No. of Particles')
    plt.xlabel('No. of part.')
    plt.ylabel('Pressure/ Pa')
    plt.legend()
    plt.show()

    chi2_05_dt=0
    for i in range(0,len(press__0_5)):
        chi2_05_dt += ((press__0_5[i]-pp_dt_0_5[i])*(press__0_5[i]-pp_dt_0_5[i]))/pp_dt_0_5[i]
    df_dt=len(pp_dt_0_5)-1
    pv_05_dt=1-stats.chi2.cdf(chi2_05_dt, df_dt)

    chi2_1_dt=0
    for i in range(0,len(press_1)):
        chi2_1_dt += ((press_1[i]-pp_dt_0_5[i])*(press_1[i]-pp_dt_0_5[i]))/pp_dt_0_5[i]
    pv_1_dt=1-stats.chi2.cdf(chi2_1_dt, df_dt)

    chi2_05_dv=0
    for i in range(0,len(press_dv_0_5)):
        chi2_05_dv += ((press_dv_0_5[i]-press_pred_dv[i])*(press_dv_0_5[i]-press_pred_dv[i]))/press_pred_dv[i]
    df_dv=len(press_pred_dv)-1
    pv_05_dv=1-stats.chi2.cdf(chi2_05_dv,df_dv)

    chi2_1_dv=0
    for i in range(0,len(press_dv_1)):
        chi2_1_dv += ((press_dv_1[i]-press_pred_dv[i])*(press_dv_1[i]-press_pred_dv[i]))/press_pred_dv[i]
    pv_1_dv=1-stats.chi2.cdf(chi2_1_dv, df_dv)

    chi2_05_dn=0
    for i in range(0,len(press_dn_0_5)):
        chi2_05_dn += ((press_dn_0_5[i]-press_pred_dn[i])*(press_dn_0_5[i]-press_pred_dn[i]))/press_pred_dn[i]
    df_dn=len(press_pred_dn)-1
    pv_05_dn=1-stats.chi2.cdf(chi2_1_dv, df_dn)

    chi2_1_dn=0
    for i in range(0,len(press_dn_1)):
        chi2_1_dn += ((press_dn_1[i]-press_pred_dn[i])*(press_dn_1[i]-press_pred_dn[i]))/press_pred_dn[i]
    pv_1_dn=1-stats.chi2.cdf(chi2_1_dn, df_dn)

    print('0.5m dT pval', pv_05_dt)
    print('0.5m dT Chi2',  chi2_05_dt)
    print('1m dT pval', pv_1_dt)
    print('1m dT Chi2',  chi2_1_dt)
    print('0.5m dV pval', pv_05_dv)
    print('0.5m dV Chi2', chi2_05_dv)
    print('1m dV pval', pv_1_dv)
    print('1m dV Chi2', chi2_1_dv)
    print('0.5m dN pval', pv_05_dn)
    print('0.5m dN Chi2', chi2_05_dn)
    print('1m dN pval', pv_1_dn)
    print('1m dN Chi2', chi2_1_dn)

    return [a,b,c]


def task14():
    """
    Task 14.

    Returns:
        Figure: The temperature ratio figure.
    """

    radii=np.linspace(0.01,1, 100)
    t_i=[]
    t_rat=[]

    for n in range(0,len(radii)):
        mbs=MultiBallSimulation(
                                c_radius=10,
                                b_radius=radii[n],b_speed=10,
                                b_mass=1,rmax=8,nrings=3,multi=6)
        mbs.run(500, animate=False, pause_time=0)
        t_i.append(mbs.t_ideal())
        t_e=mbs.t_equipartition()
        Container.dp_tot1 = 0
    for m in range(0,len(t_i)):
        rat=t_e/(t_i[m])
        t_rat.append(rat)
    a=plt.plot(radii, t_rat)
    plt.xlabel('radius/ m')
    plt.ylabel('Te/Ti')
    plt.title('Te/Ti vs rad')
    plt.show()
    return a


def task15():
    """
    Task 15.
    
    Returns:
        Figure: The speed histogram.
    """
    speed=[]
    k=1.380649E-23
    mbs = MultiBallSimulation(c_radius=20, b_radius=0.1, b_speed=10, b_mass=1,rmax=18, nrings=6, multi=6)
    mbs.run(2000, animate=False, pause_time=0)
    speed=mbs.speeds()
    speed.sort()
    kbt=k*mbs.t_ideal()
    mass=mbs.balls()[0].mass()
    mass_adj=mass*len(mbs.balls())
    p_10=[]
    for i in range(0,len(speed)):
        prob=len(mbs.balls())*maxwell(speed=speed[i], kbt=kbt, mass=mass)
        p_10.append(prob)
    plt.plot(speed, p_10, label= '10 m/s')
    counts, bin_edges=np.histogram(speed, bins=int(max(speed)))
    #bin_centers= (bin_edges[:-1] + bin_edges[1:])/2
    # print("Speeds", speed)
    #print("Boltz", boltz)
    #print('Bin centers', bin_centers)
    g=plt.hist(speed, bins=int(max(speed)), label='10 m/s')
    mbs_20 = MultiBallSimulation(c_radius=20,b_radius=0.1,b_speed=20,b_mass=1,rmax=18,nrings=6,multi=6)
    mbs_20.run(2000, animate=False, pause_time=0)
    mass=mbs_20.balls()[0].mass()
    mass_adj=len(mbs_20.balls())*mass
    speed_20=mbs_20.speeds()
    speed_20.sort()
    print(speed_20)
    kbt_20=k*mbs_20.t_ideal()
    p_20=[]
    for i in range(0,len(speed_20)):
        prob=len(mbs_20.balls())*maxwell(speed=speed_20[i], kbt=kbt_20, mass=mass)
        p_20.append(prob)
    plt.plot(speed_20, p_20, label='20 m/s')
    h=plt.hist(speed_20, bins=int(max(speed_20)), label='20 m/s')



    mbs_30=MultiBallSimulation(c_radius=20,b_radius=0.1,b_speed=30,b_mass=1,rmax=18,nrings=6,multi=6)
    mbs_30.run(2000, animate=False, pause_time=0)
    speed_30=mbs_30.speeds()
    speed_30.sort()
    mass_adj=mass*len(mbs_30.balls())
    kbt_30=k*mbs_30.t_ideal()
    mass=mbs_30.balls()[0].mass()
    p_30=[]
    for i in range(0,len(speed_30)):
        prob=len(mbs_30.balls())*maxwell(speed=speed_30[i], kbt=kbt_30, mass=mass)
        p_30.append(prob)
    plt.plot(speed_30, p_30, label='30 m/s')
    j=plt.hist(speed_30, bins=int(max(speed_30)), label='30 m/s')
    plt.xlabel('speed /ms^-1')
    plt.ylabel('No. of Balls')
    plt.legend()
    plt.show()
    return g,h,j



if __name__ == "__main__":

    # Run task 9 function
    #BALL_POS, BALL_VEL = task9()

    # Run task 10 function
    # task10()

    # Run task 11 function
    # FIG11_BALLCENTRE, FIG11_INTERBALL = task11()

    # Run task 12 function
    # FIG12_KE, FIG12_MOMX, FIG12_MOMY, FIG12_PT = task12()

    # Run task 13 function
    # FIG13_PT, FIG13_PV, FIG13_PN = task13()

    # Run task 14 function
    # FIG14 = task14()

    # Run task 15 function
    FIG15 = task15()

    plt.show()