"""
Runs simulation for single and multiballs
"""
import numpy as np
import matplotlib.pyplot as plt
from thermosnooker.balls import Container
from thermosnooker.balls import Ball
class Simulation:
    """
    Backbone of .run
    """
    def next_collision(self):
        """
        next_collision() as defined in SingleBalls.

        Raises:
            NotImplementedError: error raised if next_collision() does not exist."""
        raise NotImplementedError('next_collision() needs to be implemented in derived classes')
    def setup_figure(self):
        """
        setup_figure() as defined in SingleBalls.

        Raises:
            NotImplementedError: error raised if setup_figure() does not exist."""
        raise NotImplementedError('setup_figure() needs to be implemented in derived classes')
    def run(self, num_collisions, animate=False, pause_time=0.001):
        """
        performs simulations
        
        Args:
            num_collisions: total number of next_collision()s called.
            animate: whether to produce figure or not.
            pause_time: speed of figure animation if produced."""
        if animate:
            self.setup_figure()
        for _ in range(num_collisions):
            self.next_collision()
            if animate:
                plt.pause(pause_time)
        if animate:
            plt.show()
class SingleBallSimulation(Simulation):
    """
    Class for set up of singelball simulation
    """
    t=0
    def __init__(self, container, ball):
        """
        initializes singelballs

        Args:
            container: contianer class object.
            ball: ball class object."""
        self._c=container
        self._b=ball
    def container(self):
        """
        Args:
            No parameters required

        Returns:
            Container: container obj as inputted"""
        return self._c
    def ball(self):
        """
        Args:
            No parameters required

        Returns:
            Ball: ball obj as inputed"""
        return self._b
    def setup_figure(self):
        """
        Args:
            No parameters required."""
        rad = self.container().radius()
        ax = plt.axes(xlim=(-rad, rad), ylim=(-rad, rad))
        ax.add_artist(self.container().patch())
        ax.add_patch(self.ball().patch())
    def next_collision(self):
        """
        Simulate the next collision
        """
        dt=self._c.time_to_collision(self._b)
        SingleBallSimulation.t += dt
        self._b.move(dt)
        self._c.move(dt)
        self._c.collide(self._b)

class MultiBallSimulation(Simulation):
    """
    MultiBalls
    """
    def __init__(self, c_radius = 10, b_radius = 1., b_speed = 10, b_mass = 1, rmax = 8, nrings = 3, multi = 6):
        """
        initializes singelballs

        Args:
            c_radius: radius of contianer class object.
            b_radius: radius of ball class object.
            b_speed: speed of ball class object.
            b_mass: mass of ball class object.
            rmax: maximum distance from centre to centre of balls.
            nrings: no. of rings of balls present.
            multi: multiplier of no of balls per ring"""
        self._crad=float(c_radius)
        self._brad=float(b_radius)
        self._bv=float(b_speed)
        self._bm=float(b_mass)
        self._rmax=rmax
        self._nrings=nrings
        self._multi=multi
        self._lballs=self.make_balls()
        self._t=0
        self._c=Container(radius=self._crad, mass=10000000)
    def container(self):
        """
        Args:
            No paramters required.

        Returns:
            Container: container obj as inputted"""
        return self._c
    def make_balls(self):
        """
        make the balls
    
        Returns:
            List: ball objs numerous as parameters nrings and multi"""
        lballs=[]
        #test=[]
        def rtrings(rmax, nrings, multi):

            for n in range(0,nrings+1):

                npoints=n*multi
                if npoints !=0:
                    theta=(2*np.pi)/npoints
                    rad=(rmax/nrings)*n
                    for i in range(0,npoints):
                        rand=np.random.rand(2)
                        m_rand=np.sqrt((rand[0]*rand[0])+(rand[1]*rand[1]))
                        n_rand=rand/m_rand
                        lballs.append( Ball(pos=[(rad*np.cos(i*theta)),(rad*np.sin(i*theta))], vel=self._bv*(n_rand), radius=self._brad, mass=self._bm))

                else:
                    rad=(rmax/nrings)*n
                    theta=0
                    rand=np.random.rand(2)
                    m_rand=np.sqrt((rand[0]*rand[0])+(rand[1]*rand[1]))
                    n_rand=rand/m_rand
                    lballs.append( Ball(pos=[rad*np.cos(theta),rad*np.sin(theta)], vel=self._bv*(n_rand) ,radius=self._brad, mass=self._bm))

        rtrings(rmax=self._rmax, nrings=self._nrings, multi=self._multi)

        return lballs
    def balls(self):
        """
        Args:
            No parameters required

        Returns:
            List: ball objs numerous as parameters nrings and multi"""
        return self._lballs
    def setup_figure(self):
        """
        Args:
            No parameters required."""
        rad = self.container().radius()
        # fig = plt.figure()
        ax = plt.axes(xlim=(-rad, rad), ylim=(-rad, rad))
        ax.add_artist(self.container().patch())
        for n in self.balls():
            ax.add_patch(n.patch())
    def next_collision(self):
        """
        Args:
            No parameters required."""

        min_dt=np.inf
        con_col=False
        colliders=[None,None]
        colliders=np.array(colliders)
        for n in range(0,len(self.balls())):
            for m in range(0,len(self.balls())):
                if n!=m:
                    dt=self.balls()[n].time_to_collision(self.balls()[m])
                    if dt is None:
                        continue
                    elif dt < min_dt:
                        min_dt=dt
                        colliders[0]=self.balls()[n]
                        colliders[1]=self.balls()[m]
                        con_col = False

            dt_c=self.container().time_to_collision(self.balls()[n])
            if dt_c is not None:
                dt_c=dt_c
            else:
                dt_c=0
            if dt_c < min_dt:

                min_dt = dt_c
                colliders[0]=self.balls()[n]
                colliders[1]=self.container()
                con_col = True

                
        for i in range(0,len(self.balls())):
            self.balls()[i].move(min_dt)
            self.container().move(min_dt)
        self._t += min_dt
        if con_col is True:
            self.container().collide(colliders[0])
        else:
            colliders[0].collide(colliders[1])

    def kinetic_energy(self):
        """
        Args:
            No parameters required.

        Returns:
            Float: Total kinetic energy of the entire system."""
        Ek=0
        #end
        for i in range(0,len(self.balls())):
            v=self.balls()[i].vel()
            mag2=(v[0]*v[0])+(v[1]*v[1])
            Eki=0.5*self.balls()[0].mass()*mag2
            Ek += Eki
        vc=self.container().vel()
        magc=(vc[0]*vc[0])+(vc[1]*vc[1])
        Ek_c=0.5*self.container().mass()*magc
        Ek += Ek_c
        return Ek
    def momentum(self):
        """
        Args:
            No parameters required.

        Returns:
            Float: Total magnitude of momentum of the entire system."""
        p=0
        for i in range(0,len(self.balls())):
            v=self.balls()[i].vel()
            pi=self.balls()[i].mass()*v
            p += pi

        p_c=self.container().vel()*self.container().mass()
        p += p_c
        return p
    def time(self):
        """
        Args:
            No parameters required

        Returns:
            Float: Running total time of motion of particles."""
        return self._t
    def pressure(self):
        """
        Args:
            No parameters required

        Returns:
            Float: Total pressure inside container after certain time of motion of particles."""
        dpc=self.container().dp_tot()
        a=self.container().surface_area()
        if self.time() == 0:
            f=0
        else:
            f=dpc/self.time()
        pres=f/a
        return pres
    def t_equipartition(self):
        """
        Args:
            No parameters required

        Returns:
            Float: Tempertaure of system via equipartion."""
        k=1.380649E-23
        E_avg= (self.kinetic_energy())/len(self.balls())
        temp=E_avg/k
        return temp
    def t_ideal(self):
        """
        Args:
            No parameters required

        Returns:
            Float: Ideal gas temperature of system."""
        k=1.380649E-23
        p=self.pressure()
        v=self.container().volume()
        n=len(self.balls())
        temp_i=(p*v)/(n*k)
        return temp_i
    def speeds(self):
        """
        Args:
            No parameters required

        Returns:
            List: Magnitude of speed of all particles in the system."""
        sped=[]
        for i in range(0,len(self.balls())):
            v=self.balls()[i].vel()
            mag_v=np.sqrt((v[0]*v[0])+(v[1]*v[1]))
            sped.append(mag_v)
        return sped