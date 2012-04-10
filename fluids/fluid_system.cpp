/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#include <conio.h>

#ifdef _MSC_VER
	#include <gl/glut.h>
#else
	#include <GL/glut.h>
#endif

#include "common_defs.h"
#include "mtime.h"
#include "fluid_system.h"

#define EPSILON			0.00001f			//for collision detection

FluidSystem::FluidSystem ()
{
}

void FluidSystem::Initialize ( int mode, int total )
{
	if ( mode != BFLUID ) {
		printf ( "ERROR: FluidSystem not initialized as BFLUID.\n");
	}
	PointSet::Initialize ( mode, total );
	
	FreeBuffers ();
	AddBuffer ( BFLUID, sizeof ( Fluid ), total );
	AddAttribute ( 0, "pos", sizeof ( Vector3DF ), false );	
	AddAttribute ( 0, "color", sizeof ( DWORD ), false );
	AddAttribute ( 0, "vel", sizeof ( Vector3DF ), false );
	AddAttribute ( 0, "ndx", sizeof ( unsigned short ), false );
	AddAttribute ( 0, "age", sizeof ( unsigned short ), false );

	AddAttribute ( 0, "pressure", sizeof ( double ), false );
	AddAttribute ( 0, "density", sizeof ( double ), false );
	AddAttribute ( 0, "sph_force", sizeof ( Vector3DF ), false );
	AddAttribute ( 0, "next", sizeof ( Fluid* ), false );
	AddAttribute ( 0, "tag", sizeof ( bool ), false );

	AddAttribute ( 0, "temp", sizeof ( float ), false );
    AddAttribute ( 0, "state", sizeof ( enum Status ), false );
	AddAttribute ( 0, "mass", sizeof ( float ), false );
    AddAttribute ( 0, "adjacents", sizeof ( int ), false );

	SPH_Setup ();
	Reset ( total );
}

void FluidSystem::Reset ( int nmax )
{
	ResetBuffer ( 0, nmax );

	m_DT = 0.003; //  0.001;			// .001 = for point grav

	// Reset parameters
	m_Param [ MAX_FRAC ] = 1.0;
	m_Param [ POINT_GRAV ] = 0.0;
	m_Param [ PLANE_GRAV ] = 1.0;

	m_Param [ BOUND_ZMIN_SLOPE ] = 0.0;
	m_Param [ FORCE_XMAX_SIN ] = 0.0;
	m_Param [ FORCE_XMIN_SIN ] = 0.0;	
	m_Toggle [ WRAP_X ] = false;
	m_Toggle [ WALL_BARRIER ] = false;
	m_Toggle [ LEVY_BARRIER ] = false;
	m_Toggle [ DRAIN_BARRIER ] = false;
	m_Param [ SPH_INTSTIFF ] = 1.00;
	m_Param [ SPH_VISC ] = 0.2;
	m_Param [ SPH_INTSTIFF ] = 0.50;
	m_Param [ SPH_EXTSTIFF ] = 20000;
	m_Param [ SPH_SMOOTHRADIUS ] = 0.01;
	
	m_Vec [ POINT_GRAV_POS ].Set ( 0, 0, 0 );
	m_Vec [ PLANE_GRAV_DIR ].Set ( 0, 0, -9.8 );
	m_Vec [ EMIT_POS ].Set ( 0, 0, 0 );
	m_Vec [ EMIT_RATE ].Set ( 0, 0, 0 );
	m_Vec [ EMIT_ANG ].Set ( 0, 90, 1.0 );
	m_Vec [ EMIT_DANG ].Set ( 0, 0, 0 );
}

int FluidSystem::AddPoint ()
{
	xref ndx;	
	Fluid* f = (Fluid*) AddElem ( 0, ndx );	
	f->sph_force.Set(0,0,0);
	f->vel.Set(0,0,0);
	f->vel_eval.Set(0,0,0);
	f->next = 0x0;
	f->pressure = 0;
	f->density = 0;
    f->temp = ICE_T;
    f->state = LIQUID; //SOLID;
    f->mass = 0; // mucho problem?
	return ndx;
}

int FluidSystem::AddPointReuse ()
{
	xref ndx;
	Fluid* f;
    if ( NumPoints() <= mBuf[0].max-2 ) {
		f = (Fluid*) AddElem ( 0, ndx );
    } else {
		f = (Fluid*) RandomElem ( 0, ndx );
    }

	f->sph_force.Set(0,0,0);
	f->vel.Set(0,0,0);
	f->vel_eval.Set(0,0,0);
	f->next = 0x0;
	f->pressure = 0;
	f->density = 0; 
	f->temp = ICE_T;
    f->state = SOLID;
	f->mass = 1; // mucho problem?
	return ndx;
}

void FluidSystem::Run ()
{
	bool bTiming = true;

	mint::Time start, stop;
	
	//float ss = vgrid->voxelSize[0]*2;// m_Param [ SPH_PDIST ] / m_Param[ SPH_SIMSCALE ];		// simulation scale (not Schutzstaffel)
	
	#ifdef NOGRID
		// Slow method - O(n^2)
		SPH_ComputePressureSlow ();
		SPH_ComputeForceSlow ();
	#else
    // -- CPU only --


    start.SetSystemTime ( ACC_NSEC );
    Grid_InsertParticles ();
    //if ( bTiming) { stop.SetSystemTime ( ACC_NSEC ); stop = stop - start; printf ( "INSERT: %s\n", stop.GetReadableTime().c_str() ); }
		
    start.SetSystemTime ( ACC_NSEC );
    SPH_ComputePressureGrid ();
    //if ( bTiming) { stop.SetSystemTime ( ACC_NSEC ); stop = stop - start; printf ( "PRESS: %s\n", stop.GetReadableTime().c_str() ); }

    start.SetSystemTime ( ACC_NSEC );
    SPH_ComputeForceGridNC ();
    //if ( bTiming) { stop.SetSystemTime ( ACC_NSEC ); stop = stop - start; printf ( "FORCE: %s\n", stop.GetReadableTime().c_str() ); }

    SPH_DrawDomain();

    start.SetSystemTime ( ACC_NSEC );
    Advance();
    //if ( bTiming) { stop.SetSystemTime ( ACC_NSEC ); stop = stop - start; printf ( "ADV: %s\n", stop.GetReadableTime().c_str() ); }
		
	#endif
}



void FluidSystem::SPH_DrawDomain ()
{
	Vector3DF min, max;
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	min.z += 0.5;

	glColor3f ( 0.0, 0.0, 1.0 );
	glBegin ( GL_LINES );
	glVertex3f ( min.x, min.y, min.z );	glVertex3f ( max.x, min.y, min.z );
	glVertex3f ( min.x, max.y, min.z );	glVertex3f ( max.x, max.y, min.z );
	glVertex3f ( min.x, min.y, min.z );	glVertex3f ( min.x, max.y, min.z );
	glVertex3f ( max.x, min.y, min.z );	glVertex3f ( max.x, max.y, min.z );

	glVertex3f ( min.x, min.y, max.z );	glVertex3f ( max.x, min.y, max.z );
	glVertex3f ( min.x, max.y, max.z );	glVertex3f ( max.x, max.y, max.z );
	glVertex3f ( min.x, min.y, max.z );	glVertex3f ( min.x, max.y, max.z );
	glVertex3f ( max.x, min.y, max.z );	glVertex3f ( max.x, max.y, max.z );

    glVertex3f ( min.x, min.y, min.z ); glVertex3f ( min.x, min.y, max.z );
    glVertex3f ( min.x, max.y, min.z ); glVertex3f ( min.x, max.y, max.z );
    glVertex3f ( max.x, min.y, min.z ); glVertex3f ( max.x, min.y, max.z );
    glVertex3f ( max.x, max.y, min.z ); glVertex3f ( max.x, max.y, max.z );
	glEnd ();
}

void FluidSystem::Advance ()
{
	char *dat1, *dat1_end;
	Fluid* p;
	Vector3DF norm, z;
	Vector3DF dir, accel;
	Vector3DF vnext;
	Vector3DF min, max;
	double adj;
	float SL, SL2, ss, radius;
	float stiff, damp, speed, diff; 
	SL = m_Param[SPH_LIMIT];
	SL2 = SL*SL;
	
	stiff = m_Param[SPH_EXTSTIFF];
	damp = m_Param[SPH_EXTDAMP];
	radius = m_Param[SPH_PRADIUS];
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	ss = 0.003;//m_Param[SPH_SIMSCALE];

	dat1_end = mBuf[0].data + NumPoints()*mBuf[0].stride;
	for ( dat1 = mBuf[0].data; dat1 < dat1_end; dat1 += mBuf[0].stride ) {
		p = (Fluid*) dat1;		

		// Compute Acceleration
		accel = p->sph_force;
		accel *= m_Param[SPH_PMASS];

		// Velocity limiting 
		speed = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
		if ( speed > SL2 ) {
			accel *= SL / sqrt(speed);
		}
	
		// Boundary Conditions

		// Z-axis walls
		diff = 2 * radius - ( p->pos.z - min.z - (p->pos.x - m_Vec[SPH_VOLMIN].x) * m_Param[BOUND_ZMIN_SLOPE] )*ss;
		if (diff > EPSILON ) {			
			norm.Set ( -m_Param[BOUND_ZMIN_SLOPE], 0, 1.0 - m_Param[BOUND_ZMIN_SLOPE] );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}		

		diff = 2 * radius - ( max.z - p->pos.z )*ss;
		if (diff > EPSILON) {
			norm.Set ( 0, 0, -1 );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		
		// X-axis walls
		if ( !m_Toggle[WRAP_X] ) {
			diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1+(p->pos.y*0.025)*0.25) * m_Param[FORCE_XMIN_SIN] )*ss;	
			//diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
			if (diff > EPSILON ) {
				norm.Set ( 1.0, 0, 0 );
				adj = (m_Param[ FORCE_XMIN_SIN ] + 1) * stiff * diff - damp * norm.Dot ( p->vel_eval ) ;
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;					
			}

			diff = 2 * radius - ( max.x - p->pos.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMAX_SIN] )*ss;	
			if (diff > EPSILON) {
				norm.Set ( -1, 0, 0 );
				adj = (m_Param[ FORCE_XMAX_SIN ]+1) * stiff * diff - damp * norm.Dot ( p->vel_eval );
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
			}
		}

		// Y-axis walls
		diff = 2 * radius - ( p->pos.y - min.y )*ss;			
		if (diff > EPSILON) {
			norm.Set ( 0, 1, 0 );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		diff = 2 * radius - ( max.y - p->pos.y )*ss;
		if (diff > EPSILON) {
			norm.Set ( 0, -1, 0 );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// Wall barrier
		if ( m_Toggle[WALL_BARRIER] ) {
			diff = 2 * radius - ( p->pos.x - 0 )*ss;					
			if (diff < 2*radius && diff > EPSILON && fabs(p->pos.y) < 3 && p->pos.z < 10) {
				norm.Set ( 1.0, 0, 0 );
				adj = 2*stiff * diff - damp * norm.Dot ( p->vel_eval ) ;	
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;					
			}
		}
		
		// Levy barrier
		if ( m_Toggle[LEVY_BARRIER] ) {
			diff = 2 * radius - ( p->pos.x - 0 )*ss;					
			if (diff < 2*radius && diff > EPSILON && fabs(p->pos.y) > 5 && p->pos.z < 10) {
				norm.Set ( 1.0, 0, 0 );
				adj = 2*stiff * diff - damp * norm.Dot ( p->vel_eval ) ;	
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;					
			}
		}
		// Drain barrier
		if ( m_Toggle[DRAIN_BARRIER] ) {
			diff = 2 * radius - ( p->pos.z - min.z-15 )*ss;
			if (diff < 2*radius && diff > EPSILON && (fabs(p->pos.x)>3 || fabs(p->pos.y)>3) ) {
				norm.Set ( 0, 0, 1);
				adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
			}
		}

		// Plane gravity
		if ( m_Param[PLANE_GRAV] > 0) 
			accel += m_Vec[PLANE_GRAV_DIR];

		// Point gravity
		if ( m_Param[POINT_GRAV] > 0 ) {
			norm.x = ( p->pos.x - m_Vec[POINT_GRAV_POS].x );
			norm.y = ( p->pos.y - m_Vec[POINT_GRAV_POS].y );
			norm.z = ( p->pos.z - m_Vec[POINT_GRAV_POS].z );
			norm.Normalize ();
			norm *= m_Param[POINT_GRAV];
			accel -= norm;
		}

		// Leapfrog Integration ----------------------------
		vnext = accel;							
		vnext *= m_DT;
		vnext += p->vel;						// v(t+1/2) = v(t-1/2) + a(t) dt
		p->vel_eval = p->vel;
		p->vel_eval += vnext;
		p->vel_eval *= 0.5;					// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
		p->vel = vnext;
		vnext *= m_DT/ss;
		p->pos += vnext;						// p(t+1) = p(t) + v(t+1/2) dt

		p->temp += p->temp_eval *m_DT;
		//p->temp *= m_DT/ss;

		if ( m_Param[CLR_MODE]==1.0 ) {
			adj = fabs(vnext.x)+fabs(vnext.y)+fabs(vnext.z) / 7000.0;
			adj = (adj > 1.0) ? 1.0 : adj;
			p->clr = COLORA( 0, adj, adj, 1 );
		}
		if ( m_Param[CLR_MODE]==2.0 ) {
			float v = 0.5 + ( p->pressure / 1500.0); 
			if ( v < 0.1 ) v = 0.1;
			if ( v > 1.0 ) v = 1.0;
			p->clr = COLORA ( v, 1-v, 0, 1 );
		}

		if (m_Param[CLR_MODE] == 0.0) {
			float v = p->temp;
            float dv = MAX_T - MIN_T;
            float rgba[4] = {1.0f, 1.0f, 1.0f, 1.0f};

            if (v < MIN_T) v = MIN_T;
            if (v > MAX_T) v = MAX_T;

            if (v < (MIN_T + 0.25 * dv)) {
                rgba[0] = 0.0; //no red
                rgba[1] = (v - MIN_T) / dv; // interpolate greens
            } else if (v < (MIN_T + 0.5 * dv)) {
                rgba[0] = 0.0; //no red
                rgba[2] = (v - (MIN_T + 0.25*dv)) / dv; // interpolate blues
            } else if (v < (MIN_T + 0.75 * dv)) {
                rgba[2] = 0.0; //no blue
                rgba[1] = (v - (MIN_T + 0.5*dv)) / dv; // interpolate greens
            } else {
                rgba[2] = 0.0; //no blue
                rgba[1] = (v - (MIN_T + 0.75*dv)) / dv; // interpolate greens
            }
/*
            if (v < (MIN_T + 0.25 * dv)) {
                rgba[0] = 0.0;
                rgba[1] = 4 * (v - MIN_T) / dv;
                //rgba[1] = (v - MIN_T) / dv;
            } else if (v < (MIN_T + 0.5 * dv)) {
                rgba[0] = 0.0;
                rgba[2] = 1.0 + 4.0 * (MIN_T + 0.25 * dv - v) / dv;
                //rgba[2] = (MIN_T + 0.25 * dv - v) / dv;
            } else if (v < (MIN_T + 0.75 * dv)) {
                rgba[0] = 4.0 * (v - MIN_T - 0.5 * dv) / dv;
                //rgba[0] = (v - MIN_T - 0.5 * dv) / dv;
                rgba[2] = 0.0;
            } else {
                rgba[1] = 1.0 + 4.0 * (MIN_T + 0.75 * dv - v) / dv;
                //rgba[1] = (MIN_T + 0.75 * dv - v) / dv;
                rgba[2] = 0.0;
            }
*/
            p->clr = COLORA(rgba[0], rgba[1], rgba[2], rgba[3]);
        }

		// Euler integration -------------------------------
		/* accel += m_Gravity;
		accel *= m_DT;
		p->vel += accel;				// v(t+1) = v(t) + a(t) dt
		p->vel_eval += accel;
		p->vel_eval *= m_DT/d;
		p->pos += p->vel_eval;
		p->vel_eval = p->vel;  */	


		if ( m_Toggle[WRAP_X] ) {
			diff = p->pos.x - (m_Vec[SPH_VOLMIN].x + 2);			// -- Simulates object in center of flow
			if ( diff <= 0 ) {
				p->pos.x = (m_Vec[SPH_VOLMAX].x - 2) + diff*2;
				p->pos.z = 10;
			}
		}
	}

	m_Time += m_DT;
}

//------------------------------------------------------ SPH Setup 
//
//  Range = +/- 10.0 * 0.006 (r) =	   0.12			m (= 120 mm = 4.7 inch)
//  Container Volume (Vc) =			   0.001728		m^3
//  Rest Density (D) =				1000.0			kg / m^3
//  Particle Mass (Pm) =			   0.00020543	kg						(mass = vol * density)
//  Number of Particles (N) =		4000.0
//  Water Mass (M) =				   0.821		kg (= 821 grams)
//  Water Volume (V) =				   0.000821     m^3 (= 3.4 cups, .21 gals)
//  Smoothing Radius (R) =             0.02			m (= 20 mm = ~3/4 inch)
//  Particle Radius (Pr) =			   0.00366		m (= 4 mm  = ~1/8 inch)
//  Particle Volume (Pv) =			   2.054e-7		m^3	(= .268 milliliters)
//  Rest Distance (Pd) =			   0.0059		m
//
//  Given: D, Pm, N
//    Pv = Pm / D			0.00020543 kg / 1000 kg/m^3 = 2.054e-7 m^3	
//    Pv = 4/3*pi*Pr^3    cuberoot( 2.054e-7 m^3 * 3/(4pi) ) = 0.00366 m
//     M = Pm * N			0.00020543 kg * 4000.0 = 0.821 kg		
//     V =  M / D              0.821 kg / 1000 kg/m^3 = 0.000821 m^3
//     V = Pv * N			 2.054e-7 m^3 * 4000 = 0.000821 m^3
//    Pd = cuberoot(Pm/D)    cuberoot(0.00020543/1000) = 0.0059 m 
//
// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale)

void FluidSystem::SPH_Setup ()
{
	m_Param [ SPH_SIMSCALE ] =		0.004;			// unit size
	m_Param [ SPH_VISC ] =			0.2;			// pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)
	m_Param [ SPH_RESTDENSITY ] =	600.0;			// kg / m^3
	m_Param [ SPH_PMASS ] =			0.00020543;		// kg
	m_Param [ SPH_PRADIUS ] =		0.002;//0.004			// m
	m_Param [ SPH_PDIST ] =			0.0059;			// m
	m_Param [ SPH_SMOOTHRADIUS ] =	0.01;			// m 
	m_Param [ SPH_INTSTIFF ] =		1.00;
	m_Param [ SPH_EXTSTIFF ] =		10000.0;
	m_Param [ SPH_EXTDAMP ] =		256.0;
	m_Param [ SPH_LIMIT ] =			200.0;			// m / s

	m_Toggle [ SPH_GRID ] =		true; false;
	m_Toggle [ SPH_DEBUG ] =	true; false;

	SPH_ComputeKernels ();
}

void FluidSystem::SPH_ComputeKernels ()
{
	m_Param [ SPH_PDIST ] = pow ( m_Param[SPH_PMASS] / m_Param[SPH_RESTDENSITY], 1/3.0 );
	m_R2 = m_Param [SPH_SMOOTHRADIUS] * m_Param[SPH_SMOOTHRADIUS];
	m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow( m_Param[SPH_SMOOTHRADIUS], 9) );	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
	m_SpikyKern = -45.0f / (3.141592 * pow( m_Param[SPH_SMOOTHRADIUS], 6) );			// Laplacian of viscocity (denominator): PI h^6
	m_LapKern = 45.0f / (3.141592 * pow( m_Param[SPH_SMOOTHRADIUS], 6) );
	//std::cout << "M_spikyKern" <<m_SpikyKern << std::endl;
	//std::cout << "M_SmoothRadius" << m_Param[SPH_SMOOTHRADIUS]<< std::endl;
}

void FluidSystem::AddVolume ( Vector3DF min, Vector3DF max, float spacing,VoxelGrid* vgrid )
{
	Vector3DF pos;
	Fluid* p;
	float dx, dy, dz;
	dx = max.x-min.x;
	dy = max.y-min.y;
	dz = max.z-min.z;

	// temp counter
	int count = 0;
	/*for (float z = max.z; z >= min.z; z -= spacing ) {
		for (float y = min.y; y <= max.y; y += spacing ) {	
			for (float x = min.x; x <= max.x; x += spacing ) {
                Vector3DF index = vgrid->inVoxelGrid(x,y,z);
				if(index.x > 0 && index.y > 0 && index.z > 0){
				    std::cout << "count " <<count << std::endl;
					count++;
					p = (Fluid*)GetPoint ( AddPointReuse () );
					pos.Set ( x, y, z);
                    p->index = index;
					//pos.x += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
					//pos.y += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
					//pos.z += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
					p->pos = pos;
					p->clr = COLORA( (x-min.x)/dx, (y-min.y)/dy, (z-min.z)/dz, 1);

				}
			}
		}
	}*/
	
	for (float z = 16; z >= 0; z -= 3.2f ) {
		for (float y = 0; y <= 16; y += 3.2f ) {	
			for (float x = 0; x <= 16; x += 3.2f ) {
                Vector3DF index = vgrid->inVoxelGrid(x,y,z);
				if(index.x > 0 && index.y > 0 && index.z > 0){
				    std::cout << "count " <<count << std::endl;
					count++;
					p = (Fluid*)GetPoint ( AddPointReuse () );
					pos.Set ( x, y, z);
                    p->index = index;
					//pos.x += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
					//pos.y += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
					//pos.z += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
					p->pos = pos;
					p->clr = COLORA( (x-min.x)/dx, (y-min.y)/dy, (z-min.z)/dz, 1);

				}
			}
		}
	}
	
}

void FluidSystem::SPH_CreateExample ( int n, int nmax ) //currently creates a cube
{
	Vector3DF min, max;

	// Testing loading dragon
	m_Vec [ SPH_VOLMIN ].Set ( VOLMIN_X, VOLMIN_Y, VOLMIN_Z );
	m_Vec [ SPH_VOLMAX ].Set ( VOLMAX_X, VOLMAX_Y, VOLMAX_Z );
	m_Vec [ SPH_INITMIN ].Set ( INITMIN_X, INITMIN_Y, INITMIN_Z);
	m_Vec [ SPH_INITMAX ].Set ( INITMAX_X, INITMAX_Y, INITMAX_Z );
	//if (vgrid)
		//delete vgrid;

	switch ( n ) {
	case 0:
		// Load cube
		vgrid = new VoxelGrid("voxel/cube_4.voxels");

		break;
	case 1:
		// Load dragon
		vgrid = new VoxelGrid("voxel/dragon.voxels");
	break;
	case 2:
		vgrid = new VoxelGrid("voxel/happy.voxels");
	break;

	}
	if(!vgrid) {
		return;
	}
 	nmax = 64;// vgrid->theDim[0] * vgrid->theDim[1] * vgrid->theDim[2];
	if (nmax > 8400)
		nmax = 8400;

	ss = vgrid->voxelSize[0];
	std::cout << "nmax" << nmax  << std::endl;

	Reset ( nmax );
    // init our adjacency list
    short neighbors;
    vgrid->adj[1][1][1] = 0;
    for (int i = 0; i < vgrid->theDim[0]; i++) {
        for (int j = 0; j < vgrid->theDim[2]; j++) {
            for (int k = 0; k < vgrid->theDim[1]; k++) {
                neighbors = 0;
                if (vgrid->data[i][j][k]) { // if there is a voxel in that location
                    if (i > 0 && vgrid->data[i-1][j][k]) neighbors++;
                    if (i < vgrid->theDim[0] - 1 && vgrid->data[i+1][j][k]) neighbors++;
                    if (j > 0 && vgrid->data[i][j-1][k]) neighbors++;
                    if (j < vgrid->theDim[2] - 1 && vgrid->data[i][j+1][k]) neighbors++;
                    if (k > 0 && vgrid->data[i][j][k-1]) neighbors++;
                    if (k < vgrid->theDim[1] - 1 && vgrid->data[i][j][k+1]) neighbors++;
                } else {
                    neighbors = -1; //error state
                }
                vgrid->adj[i][j][k] = neighbors;
            }
        }
    }

	m_Param [ SPH_SIMSIZE ] = m_Param [ SPH_SIMSCALE ] * (m_Vec[SPH_VOLMAX].z - m_Vec[SPH_VOLMIN].z);
	m_Param [ SPH_PDIST ] = pow ( m_Param[SPH_PMASS] / m_Param[SPH_RESTDENSITY], 1/3.0 );
	//float ss = vgrid->voxelSize[0]*2 ;// m_Param [ SPH_PDIST ]*0.87 / m_Param[ SPH_SIMSCALE ];
	//printf ( "Spacing: %f\n", ss);

	// Hacking for now...Need to find good mapping
	AddVolume ( m_Vec[SPH_INITMIN], m_Vec[SPH_INITMAX], vgrid->voxelSize[0], vgrid);//ss, vgrid );	// Create the particles

    Fluid* f;
	Vector3DF pos;

	float cell_size = m_Param[SPH_SMOOTHRADIUS]*2.0;			// Grid cell size (2r)
	Grid_Setup ( m_Vec[SPH_VOLMIN], m_Vec[SPH_VOLMAX], m_Param[SPH_SIMSCALE], cell_size, 1.0 ); // Setup grid
	Grid_InsertParticles ();									// Insert particles

	Vector3DF vmin, vmax;
	vmin =  m_Vec[SPH_VOLMIN];
	vmin -= Vector3DF(2,2,2);
	vmax =  m_Vec[SPH_VOLMAX];
	vmax += Vector3DF(2,2,-2);
}

// Compute Pressures - Using spatial grid, and also create neighbor table
void FluidSystem::SPH_ComputePressureGrid ()
{
	char *dat1, *dat1_end;
	Fluid* p;
	Fluid* pcurr;
	int pndx;
	int i, cnt = 0;
	float dx, dy, dz, sum, dsq, c;
	float d, d2, mR, mR2;
	float radius = m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE];
	d = m_Param[SPH_SIMSCALE];
	d2 = d*d;
	mR = m_Param[SPH_SMOOTHRADIUS];
	mR2 = mR*mR;	

	dat1_end = mBuf[0].data + NumPoints()*mBuf[0].stride;
	i = 0;
	for ( dat1 = mBuf[0].data; dat1 < dat1_end; dat1 += mBuf[0].stride, i++ ) {
		p = (Fluid*) dat1;

		sum = 0.0;	
		m_NC[i] = 0;

		Grid_FindCells ( p->pos, radius );
		for (int cell=0; cell < 8; cell++) {
			if ( m_GridCell[cell] != -1 ) {
				pndx = m_Grid [ m_GridCell[cell] ];				
				while ( pndx != -1 ) {					
					pcurr = (Fluid*) (mBuf[0].data + pndx*mBuf[0].stride);					
					if ( pcurr == p ) {pndx = pcurr->next; continue; }
					dx = ( p->pos.x - pcurr->pos.x)*d;		// dist in cm
					dy = ( p->pos.y - pcurr->pos.y)*d;
					dz = ( p->pos.z - pcurr->pos.z)*d;
					dsq = (dx*dx + dy*dy + dz*dz);
					if ( mR2 > dsq ) {
						c =  m_R2 - dsq;
						sum += c * c * c;
						if ( m_NC[i] < MAX_NEIGHBOR ) {
							m_Neighbor[i][ m_NC[i] ] = pndx;
							m_NDist[i][ m_NC[i] ] = sqrt(dsq);
							m_NC[i]++;
						}
					}
					pndx = pcurr->next;
				}
			}
			m_GridCell[cell] = -1;
		}
		p->density = sum * m_Param[SPH_PMASS] * m_Poly6Kern ;	
		p->pressure = ( p->density - m_Param[SPH_RESTDENSITY] ) * m_Param[SPH_INTSTIFF];		
		p->density = 1.0f / p->density;		
	}
}

// Compute Forces - Using spatial grid with saved neighbor table. Fastest.
void FluidSystem::SPH_ComputeForceGridNC ()
{
	char *dat1, *dat1_end;	
	Fluid *p, *pcurr;
	Vector3DF force, fcurr;
	register float pterm, vterm, dterm;
	int i;
    float edge;
	float c, d;
	float dx, dy, dz;
	float mR, mR2, visc;
    float length, lap_kern;
    float neighbor_temp; // new additions to temperature by other particles
    float sa, Qi, dT; // air contribution to new temperature
    int pi, pj, pk; 

	d = m_Param[SPH_SIMSCALE];
	mR = m_Param[SPH_SMOOTHRADIUS];
	mR2 = (mR*mR);
	visc = m_Param[SPH_VISC];

	dat1_end = mBuf[0].data + NumPoints()*mBuf[0].stride;
	edge = vgrid->voxelSize[0]*2; m_Param [ SPH_PDIST ]*0.87 / m_Param[ SPH_SIMSCALE ];
    i = 0;
	for ( dat1 = mBuf[0].data; dat1 < dat1_end; dat1 += mBuf[0].stride, i++ ) {
        // reset all instance variables
		p = (Fluid*) dat1;

		force.Set (0, 0, 0);
		neighbor_temp = 0.0; sa = 0.0; Qi = 0.0;
        pi = p->index.x;
        pj = p->index.y;
        pk = p->index.z;
		
		//std::cout << "normal loop " << std::endl;
		//std::cout << "pi " << pi << " pj " << pj << " pk " << pk << std::endl;
        if (p->state == SOLID) { // hack to prevent gravity on solids for now
            force -= m_Vec[PLANE_GRAV_DIR];
            force /= m_Param[SPH_PMASS];
        } // take out after interfacial tension is implemeneted

        for (int j=0; j < m_NC[i]; j++ ) { // loop through all neighbors
            pcurr = (Fluid*) (mBuf[0].data + m_Neighbor[i][j]*mBuf[0].stride);
            dx = ( p->pos.x - pcurr->pos.x)*d; // dist in cm
            dy = ( p->pos.y - pcurr->pos.y)*d;
            dz = ( p->pos.z - pcurr->pos.z)*d;
            c = ( mR - m_NDist[i][j] ); //distance between current and neighbor?
            pterm = -0.5f * c * m_SpikyKern * ( p->pressure + pcurr->pressure) / m_NDist[i][j];
            dterm = c * p->density * pcurr->density;
            vterm = m_LapKern * visc;

            Vector3DF diff = p->pos; // distance between particle and its neighbor
            diff -= pcurr->pos;
            length = diff.Length();
            lap_kern = m_LapKern * (m_Param[SPH_SMOOTHRADIUS] - length);
            //neighbor_temp += m_Param [ SPH_PMASS ] * ((pcurr->temp - p->temp)/pcurr->density) * lap_kern; // Newtonian Heat Transfer

            if (p->state == LIQUID) {
                force.x += ( pterm * dx + vterm * (pcurr->vel_eval.x - p->vel_eval.x) ) * dterm;
                force.y += ( pterm * dy + vterm * (pcurr->vel_eval.y - p->vel_eval.y) ) * dterm;
                force.z += ( pterm * dz + vterm * (pcurr->vel_eval.z - p->vel_eval.z) ) * dterm;
                // interfacial force equations
                //force.x += (K_W + K_ICE)*(pcurr->pos.x - p->pos.x)/(dx*dx);
                //force.y += (K_W + K_ICE)*(pcurr->pos.y - p->pos.y)/(dy*dy);
                //force.z += (K_W + K_ICE)*(pcurr->pos.z - p->pos.z)/(dz*dz);
            } else { //SOLID
            }
        }
        p->sph_force = force;
        neighbor_temp *= DIFF_T;

        // air affected temperature
        if (p->state == SOLID && vgrid->data[pi][pj][pk]) { // check surface particle?
            sa = (6.0 - vgrid->adj[pi][pj][pk])/6.0;// * (edge * edge * 6.0);
            Qi = THERMAL_CONDUCTIVITY_ICE * (AMBIENT_T - p->temp) * sa;
        } else if (p->state == LIQUID) {
            Qi = THERMAL_CONDUCTIVITY_WATER * (AMBIENT_T - p->temp);
        }
        dT = Qi / (HEAT_CAP * m_Param [ SPH_PMASS ]);
        //p->temp += dT;
        //p->temp += neighbor_temp;
        //p->temp_eval = p->temp; //what?
        p->temp_eval = neighbor_temp + dT; //what?
		if (p->temp > 0.9 && p->state == SOLID) { // change state and update neighboring voxels
            vgrid->data[pi][pj][pk] = 0; // set to no particle
			std::cout << "update neighbor " << std::endl;
			std::cout << "pi " << pi << " pj " << pj << " pk " << pk << std::endl;
           // vgrid->adj[pi][pj][pk] = -1;
			std::cout << "neighbor " << vgrid->adj[pi][pj][pk] << std::endl;
            if (pi + 1 < vgrid->theDim[0]) vgrid->adj[pi+1][pj][pk]--;
            if (pi - 1 < vgrid->theDim[0]) vgrid->adj[pi-1][pj][pk]--;
            if (pj + 1 < vgrid->theDim[2]) vgrid->adj[pi][pj+1][pk]--;
            if (pj - 1 < vgrid->theDim[2]) vgrid->adj[pi][pj-1][pk]--;
            if (pk + 1 < vgrid->theDim[1]) vgrid->adj[pi][pj][pk+1]--;
            if (pk - 1 < vgrid->theDim[1]) vgrid->adj[pi][pj][pk-1]--;
            p->state = LIQUID;
		}
	}
}