/*-----------------------------------------------------------------
 RK4:
 
 - A numerical integration technique, called Runge-Kutta 4, that
 updates the velocity and position of the particle. This is
 better than Euler Cromer and the previous RK's.
 
 Precondition : None
 Postcondition: Updates velocity and position of the particle
 -----------------------------------------------------------------*/
- (void) RK4
{
    // a() -> acceleration
    // v   -> velocity
    // x   -> position
    
    double v = velocity;
    double x = position;
    
    double  dv1 = a(t,x,v) * dt;
    double  dx1 = v*dt;
    // --------------
    
    double  dv2 = a( t+dt/2.0, x + dx1/2.0, v +dv1/2.0) * dt;
    double  dx2 = (v+dv1/2.0) * dt;
    // ---------------
    
    double  dv3 = a( t+dt/2.0, x + dx2/2.0, v +dv2/2.0) * dt;
    double  dx3 = (v+dv2/2.0) * dt;
    // ---------------
    
    double  dv4 = a( t+dt, x + dx3, v +dv3) * dt;
    double  dx4 = (v+dv3) * dt;
    // --------------
    
    double  dv = (dv1 + 2*dv2 + 2*dv3 + dv4) / 6.0;
    double  dx = (dx1 + 2*dx2 + 2*dx3 + dx4) / 6.0;
    // ---------------
    
    velocity += dv;
    position += dx;
    t +=dt;
}

