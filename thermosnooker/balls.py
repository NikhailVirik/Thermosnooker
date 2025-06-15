"""
Start of document
"""
import numpy as np
from matplotlib.patches import Circle

class Ball:
    """
    Class for creating balls
    """
    def __init__(self, pos=[0., 0.], vel=[1., 0.], radius=1., mass=1.):
        if len(pos)!=2:
            raise ValueError('only accpet 2-component array')
        else:
            self._pos=np.array(pos)
        if len(vel)!=2:
            raise ValueError('only accpet 2-component array for velocity')
        else:
            self._vel=np.array(vel)
        self._rad=float(radius)
        self._mass=float(mass)
        self._patch = Circle((self._pos[0],self._pos[1]), self._rad, fill=False, ls='solid')
    def pos(self):
        """
        Args:
            No parameters required.

        Returns:
            array: current position of obj in 2D.
        """
        return np.array(self._pos)
    def radius(self):
        """
        Args:
            No parameters required.

        Returns:
            float: Radius of obj."""
        return self._rad
    def mass(self):
        """
        Args:
            No parameters required.

        Returns:
            float: mass of obj."""
        return self._mass
    def set_vel(self, new_vel):
        """
        Args:
            new_vel (array): Updated obj velocity.
        
        Raises:
            ValueError: error raised when length."""
        if len(new_vel) != 2:
            raise ValueError('only array of tuple allowed')
        else:
            self._vel=np.array(new_vel)
    def vel(self):
        """
        Args:
            No parameters required.

        Returns:
            array: current velocity of obj in 2D."""
        return self._vel
    def patch(self):
        """
        Args:
            No parameters required.

        Returns:
            Figure: patch Figure obj from matplotlib patches."""
        return self._patch
    def move(self,dt):
        """
        Args:
            dt (float): time difference to next position of obj."""
        self._pos=np.array(self._pos + (self._vel*dt))
        self._patch.center = self._pos   
    def time_to_collision(self, other):
        """
        Args:
            other (Ball/Container): obj to collide with.

        Returns:
            float: time to next collision."""
        r = self.pos() - other.pos()
        v= self.vel() - other.vel()
        r1 = self.radius() + other.radius()
        r2 = other.radius() - self.radius()
        m_v=np.dot(v,v)
        rv=np.dot(r,v)
        m_r=np.dot(r,r)
        if m_v==0:
            return None
        elif isinstance(self, Container) or isinstance(other, Container):
            dt=((-2*rv)+(2*np.sqrt((rv*rv)-((m_v)*(m_r-(r2*r2))))))/(2*m_v)
        elif rv < 0:
            dt=((-2*rv)-(2*np.sqrt((rv*rv)-((m_v)*(m_r-(r1*r1))))))/(2*m_v)
        elif rv > 0:
            return None
        if np.isclose(dt,0):
            dt=0
        if dt<=0 or isinstance(dt, complex):
            return None
        else:
            return dt
    def collide(self,other):
        """
        Args:
            other: other object to collide with.

        Returns:
            self._vel: velocity of original obj in 2D after collision.
            other._vel: velocity of colliding obj in 2D after collision."""
        a=self.vel()
        b=other.vel()
        m1=self.mass()
        m2=other.mass()
        s=m1+m2
        d=m1-m2
        u1_cm=self.vel()
        u2_cm=other.vel()
        r=self.pos()-other.pos()
        mag_r=np.sqrt((r[0]*r[0])+(r[1]*r[1]))
        n=r/mag_r
        v1n=np.dot(u1_cm,n)*n
        v2n=np.dot(u2_cm,n)*n
        v1t=u1_cm-v1n
        v2t=u2_cm-v2n
        v1_trans=((d/s)*v1n) + ((2*m2/s)*v2n)
        v1_cm=v1t + v1_trans
        v2_trans=((self.mass()*v1n)+(other.mass()*v2n)-(self.mass()*v1_trans))/other.mass()
        v2_cm=v2t + v2_trans
        self.set_vel(v1_cm)
        other.set_vel(v2_cm)
        if isinstance(self, Container):
            dpc=self.vel()-a
            mag_dpc=self.mass()*np.sqrt((dpc[0]*dpc[0])+(dpc[1]*dpc[1]))
            Container.dp_tot1+=mag_dpc
        elif isinstance(other,Container):
            dpc=other.vel()-b
            mag_dpc=other.mass()*np.sqrt((dpc[0]*dpc[0])+(dpc[1]*dpc[1]))
            Container.dp_tot1+=mag_dpc
        return self.vel(), other.vel()
class Container(Ball):
    """
    Class for creating contaienr object
    """
    dp_tot1=0
    def __init__(self, radius= 10., mass= 10000000):
        Ball.__init__(self, pos = [0, 0], vel= [0, 0], radius = radius, mass=mass)
    def volume(self):
        """
        Args:
            No parameters required.

        Returns:
            float: volume of container in 2D."""
        vol= np.pi*self._rad*self._rad
        return vol
    def surface_area(self):
        """
        Args:
            No parameters required

        Returns:
            float: surface area of container in 2D"""
        sa=2*np.pi*self._rad
        return sa
    def dp_tot(self):
        """
        Args:
            No parameters required

        Returns:
            float: Running total magnitude change in momentum of container"""
        return Container.dp_tot1
