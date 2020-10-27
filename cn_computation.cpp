#include <iostream>  
#include <fstream> 
#include<math.h>                                                    
#include<cstdlib>   
#include<SDL/SDL.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>


using namespace std;

struct stat st = {0};

//to compile (SDL1.2 needed)
/*g++ -O3 -o tester tester.cpp graphictools.cpp `sdl-config --cflags --libs`*/


int ripet; //repetitions of experiment

double TI; //maximum experiment time
double X_MIN,Y_MIN,X_MAX,Y_MAX; //tracking area corners, assume tracking area rectangular parallell to axes
double XL,YL,XL_2,YL_2; //full and half lenght in x and y

double grid_size=0.2;//parameters, they can be changed through parameter file. 
double time_int=2.5;//time interval for velocityfield computation
double cn_radius=3.5;//ROI radius for cn computation
double graph_mr=3;//max rotor in graphs
double graph_md=12;//max density in graphs
double DPIX;//size in pixel of a cell


class Vector2D //for operations on 2d vectors
{
public:
  double x,y,m,th;  //m magnitude, th angle
  Vector2D()
  {
    x=0;
    y=0;
    m=0;
    th=0;
  }
  Vector2D(double xp,double yp)
  {
    x=xp;
    y=yp;
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    else th=0;
    if(x<0) th=-th;
  }
  bool Uguale(Vector2D U)
  {
    bool u=true;
    u=u&&(x==U.x);
    u=u&&(y==U.y);
    u=u&&(m==U.m);
    u=u&&(th==U.th);
    return u;
  }
  void Init(double xp,double yp)
  {
    x=xp;
    y=yp;
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    else th=0;
    if(x<0) th=-th;
  }
  void Set_zero()
  {
    x=0;
    y=0;
    m=0;
    th=0;
  }
  void Magnitude()
  {
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    if(x<0) th=-th;
  }
  void Scale(double s)
  {
    x*=s;
    y*=s;
    Magnitude();
  }
  void Add(Vector2D A)
  {
    x+=A.x;
    y+=A.y;
    Magnitude();
  }
};

class State2D//position and velocity
{
public:
  Vector2D r,v;
  State2D(){}
  State2D(double rx,double ry,double vx,double vy)
  {
    r.Init(rx,ry);
    v.Init(vx,vy);
  }
  State2D(Vector2D rp,Vector2D vp)
  {
    r=rp;
    v=vp;
  }
  bool Uguale(State2D U)
  {
    bool u=true;
    u=u&&r.Uguale(U.r);
    u=u&&v.Uguale(U.v);
    return u;
  }
  void Init(double rx,double ry,double vx,double vy)
  {
    r.Init(rx,ry);
    v.Init(vx,vy);
  }
  void Init(Vector2D rp,Vector2D vp)
  {
    r=rp;
    v=vp;
  }
};

class Colore //colour for graphics
{
public:
  int r,g,b;//red green blue
  Colore()
  {
    r=0;
    g=0;
    b=0;
  }
  Colore(int r_,int g_,int b_)//initializes
  {
    r=r_;
    g=g_;
    b=b_;
  }
  void Set(int r_,int g_,int b_)//sets the values
  {
    r=r_;
    g=g_;
    b=b_;
  }
  void Copy(Colore C)//copies
  {
    r=C.r;
    g=C.g;
    b=C.b;
  }
    void Init(double rd,double gd,double bd)
  {
    r=255*rd;
    g=255*gd;
    b=255*bd;
  }  
};

Colore colorenuovo(double z,double z_max,double z_min)   //colour bar minimum limit is zero
{
  Colore nuovo[17];//some colours that make sense even when watched in b&w
  nuovo[0].Init(0.,0.,0.);
  nuovo[1].Init(0.075,0.075,0.25);
  nuovo[2].Init(0.15,0.15,0.5);
  nuovo[3].Init(0.225,0.15,0.625);
  nuovo[4].Init(0.3,0.15,0.75);
  nuovo[5].Init(0.45,0.175,0.65);
  nuovo[6].Init(0.6,0.2,0.55);
  nuovo[7].Init(0.8,0.225,0.35);
  nuovo[8].Init(1.,0.25,0.15);
  nuovo[9].Init(0.95,0.375,0.075);
  nuovo[10].Init(0.9,0.5,0.);
  nuovo[11].Init(0.9,0.625,0.05);
  nuovo[12].Init(0.9,0.75,0.1);
  nuovo[13].Init(0.9,0.825,0.3);
  nuovo[14].Init(0.9,0.9,0.5);
  nuovo[15].Init(0.95,0.95,0.75);
  nuovo[16].Init(1.,1.,1.);
  if(z>z_max) z=z_max;     //I put a maximum limit
  if(z<z_min) z=z_min;
  z=z-z_min;
  int col;
  if(z_max!=z_min) col=z/(z_max-z_min)*17;
  else col=16;
  if(col==17) col=16;
  col=16-col;
  return nuovo[col];
}


 

  
void SDL_Inizio();     //graphics 
int Init_images();
int PIXEL_X;//size of 2D images in x
int PIXEL_Y;   //size of 2D images in y

SDL_Surface *screen;      //for graphics
SDL_Surface *scene;      //for graphics 

Colore *colori;   //other functions related to graphics
void DrawPixel(SDL_Surface *screen, int x, int y, Uint8 R, Uint8 G, Uint8 B);    //graphics (draws a pixel)
void DrawIMG(SDL_Surface *img, int x, int y);     //graphics
void DrawBallnp(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 b,int rg,int x1,int y1,int xmax,int ymax);//graphics (draws a ball)


void Freccia(double vx,double vy,int ci,int cj,int size) //puts an arrow in ci, cj, colour and size as in the paper
{
  int decimo=size/5;
  double vs=sqrt(vx*vx+vy*vy);
  if(vs>1)
    {
      vx/=vs;
      vy/=vs;
      vs=1;
    }
  Colore C;
  if(vx==0)
    {
      if(vy<0) {C.r=255;C.g=0;C.b=0;}
      else {C.r=0;C.g=0;C.b=255;}
    }    
  else if(vx>0)
    {
      if(vy<0)
	{
	  C.b=0;
	  double theta=atan(vy/vx);
	  C.r=-255*sin(theta);
	  C.g=255*cos(theta);
	}
      else
	{
	  C.r=0;
	  double theta=atan(vy/vx);
	  C.b=255*sin(theta);
	  C.g=255*cos(theta);
	}
    }
  else
    {
      C.g=0;
      double theta=vy/vs;
      C.b=255*(theta+1)/2;
      C.r=255*(1-theta)/2;
    }
  if(fabs(vx)>fabs(vy))
    {
      double ix_range=(size*vx);
      double ix_start=ci-ix_range/2;
      double ix_end=ix_start+ix_range;
      double ar=fabs(vs/vx);
      double lvx=vy/vs;
      double lvy=-vx/vs;
      if(ix_end>ix_start)
	{
	  for(double i=ix_start;i<=ix_end;i+=0.5)
	    {
	      double j=cj-size/2*vy+(i-ix_start)*vy/vx;
	      int ii=int(i);
	      int ij=int(j);
	      if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(i),int(j),C.r,C.g,C.b);
	      double di=fabs((ix_end-i)*ar);
	      if(di<decimo)
		{
		  double ajyr=di*lvy;
		  double ajys=j-ajyr/2;
		  double ajye=ajys+ajyr;
		  if(ajys<ajye)
		    {
		      for(double aj=ajys;aj<=ajye;aj+=0.5)
			{
			  double ai=i-di/2*lvx+(aj-ajys)*lvx/lvy;
			  ii=int(i);
			  ij=int(j);
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }
		  else
		    {
		      for(double aj=ajys;aj>=ajye;aj-=0.5)
			{
			  double ai=i-di/2*lvx+(aj-ajys)*lvx/lvy;
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }  
		}
	    }
	}
      else
	{
	  for(double i=ix_start;i>=ix_end;i-=0.5)
	    {
	      double j=cj-size/2*vy+(i-ix_start)*vy/vx;
	      int ii=int(i);
	      int ij=int(j);
	      if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(i),int(j),C.r,C.g,C.b);
	      double di=fabs((ix_end-i)*ar);
	      if(di<decimo)
		{
		  double ajyr=di*lvy;
		  double ajys=j-ajyr/2;
		  double ajye=ajys+ajyr;   
		  if(ajys<ajye)
		    {
		      for(double aj=ajys;aj<=ajye;aj+=0.5)
			{
			  double ai=i-di/2*lvx+(aj-ajys)*lvx/lvy;
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }
		  else
		    {
		      for(double aj=ajys;aj>=ajye;aj-=0.5)
			{
			  double ai=i-di/2*lvx+(aj-ajys)*lvx/lvy;			  
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }
		}
	    }
	}	
    }
  else
    {
      double jy_range=(size*vy);
      double jy_start=cj-jy_range/2;
      double jy_end=jy_start+jy_range;
      double ar=fabs(vs/vy);
      double lvx=vy/vs;
      double lvy=-vx/vs;
      if(jy_end>jy_start)
	{
	  for(double j=jy_start;j<=jy_end;j+=0.5)
	    {
	      double i=ci-size/2*vx+(j-jy_start)*vx/vy;	  
	      int ii=int(i);
	      int ij=int(j);
	      if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(i),int(j),C.r,C.g,C.b);
	      double di=fabs((jy_end-j)*ar);
	      if(di<decimo)
		{
		  double aixr=di*lvx;
		  double aixs=i-aixr/2;
		  double aixe=aixs+aixr;
		  if(aixs<aixe)
		    {
		      for(double ai=aixs;ai<=aixe;ai+=0.5)
			{
			  double aj=j-di/2*lvy+(ai-aixs)*lvy/lvx;
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }
		  else
		    {
		      for(double ai=aixs;ai>=aixe;ai-=0.5)
			{
			  double aj=j-di/2*lvy+(ai-aixs)*lvy/lvx;
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }  
		}
	    }
	}
      else
	{
	  for(double j=jy_start;j>=jy_end;j-=0.5)
	    {
	      double i=ci-size/2*vx+(j-jy_start)*vx/vy;	  
	      int ii=int(i);
	      int ij=int(j);
	      if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(i),int(j),C.r,C.g,C.b);
	      double di=fabs((jy_end-j)*ar);
	      if(di<decimo)
		{
		  double aixr=di*lvx;
		  double aixs=i-aixr/2;
		  double aixe=aixs+aixr;
		  if(aixs<aixe)
		    {
		      for(double ai=aixs;ai<=aixe;ai+=0.5)
			{
			  double aj=j-di/2*lvy+(ai-aixs)*lvy/lvx;
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }
		  else
		    {
		      for(double ai=aixs;ai>=aixe;ai-=0.5)
			{
			  double aj=j-di/2*lvy+(ai-aixs)*lvy/lvx;
			  ii=int(i);
			  ij=int(j); 
			  if((ii>0)&&(ij>0)&&(ii<PIXEL_X)&&(ij<PIXEL_Y)) DrawPixel(scene,int(ai),int(aj),C.r,C.g,C.b);
			}
		    }  
		}
	    }
	}	
    } 
}

class Stat //for statistics
{
public:
  int conta;//counter
  double av;//average
  double sg;//standard dev
  double er;//std err
  Stat(){conta=0;av=0;sg=0;er=0;}
  void Update(double up)
  {
    av+=up;
    sg+=up*up;
    conta++;
  }
  void Final() //computes everything
  {
    if(conta) 
      {	
	av/=conta;
	sg/=conta;
	if(conta>1)
	  {
	    sg=sqrt(sg-av*av);
	    er=sg/sqrt(conta-1);
	  }
	else
	  {
	    sg=0;
	    er=0;
	  }
      }
  }
  void Azzera()
  {
    av=0;
    sg=0;
    er=0;
    conta=0;
  }
};

class DDistr// includes a vector and a matrix for statistics
{
public:
  Stat **dd;//matrix, runs over time and repetitions
  Stat *d;//vector, only repetitions
  double max;//maximum over time and repetitions
  int imr,imt;//repetition and time of max
  double maxt,delta_t,delta_t_2;//maximum time and time steps for cn computation
  int time_steps;//lenght of time vector
  int ripet;//number of repetitions
  DDistr(){}
  DDistr(int r,double mt,double dt)//intialises giving # of reps, maximum time and time step
  {
    ripet=r;
    maxt=mt;
    delta_t=dt;
    delta_t_2=delta_t/2;    
    time_steps=int(maxt/delta_t);
    d=new Stat[time_steps];//allocates
    dd=new Stat*[ripet];
    for(int i=0;i<ripet;i++) dd[i]=new Stat[time_steps];
    max=-1e10;//very negative max initialisation
  }
  void Init(int r,double mt,double dt)//intialises giving # of reps, maximum time and time step, as above
  {
    ripet=r;
    maxt=mt;
    delta_t=dt;
    delta_t_2=delta_t/2;       
    time_steps=int(maxt/delta_t);
    d=new Stat[time_steps];
    dd=new Stat*[ripet];
    for(int i=0;i<ripet;i++) dd[i]=new Stat[time_steps];
    max=-1e10;
  }
  void Update(double up,int r,double t)//adds up to the statistics
  {
    int it=int(t/delta_t);
    dd[r][it].Update(up);
  }
  void Final()//finalises
  {
    for(int t=0;t<time_steps;t++)
      {
	for(int r=0;r<ripet;r++)
	  {
	    if(dd[r][t].av>max)//finds max
	      {
		max=dd[r][t].av;
		imr=r;
		imt=t;
	      }
	    d[t].Update(dd[r][t].av);//puts in vector the average over reps
	  }
	d[t].Final();//stats over reps, now in d we have av,sg, and er over reps depending on time
      }
  }
  void Print(char name[200])//prints the vector on a file with name
  {
    ofstream out;
    out.open(name);
    for(int i=0;i<time_steps;i++) out << i*delta_t+delta_t_2 << " " << d[i].av-d[i].er << " " << d[i].av << " " << d[i].av+d[i].er << endl;
    out.close();
  }
};
 

class CN   //local field values
{
public:
  Vector2D v;   //velocity field
  double dens;  //local density
  double rot;  //local rot
  bool rotval;  //if rot could be computed
  double cn;  //local cn
  int count;  //counter
  double step2;  //square of grid size
  CN()
  {
    rot=0;
    rotval=false;
    cn=0;
    count=0;
    dens=0;
    step2=0;
  }
  void Init(double s)//initialises with grid size
  {
    rot=0;
    rotval=false;
    cn=0;
    count=0;
    dens=0;
    step2=s;
  }    
  void Final()//finalises
  {
    if(dens>0)
      {
	v.Scale(1./dens);//dens is used first to count how many peds contributed to cell, and to normalise velocity field
	dens/=count*step2;//then scaled as density (count is the number of events in that cell, basically 2.5s/(rate of tracking)
      }
  }
};



class Cells  //all the fields
{
public:
  int ripet;  //# of reps
  int time_steps;  //# of time steps
  int size_x;//# number of cells in 1D x direction
  int size_y;//# number of cells in 1D y direction
  double step;//size of grid
  double step2;//2*size of grid
  double delta_t,delta_t_2;//lenght of time cell, and half
  double x_min,x_max,y_min,y_max;//borders of tracking area
  double T_max;//maximum time
  CN ****cells;//the fields in reps,time,x,y
  bool **in_area;  //false if cell never used, used to define in area cells
  DDistr dens;//for average density stats
  DDistr max_cn;//for max cn stats stats
  DDistr av_cn;//for av cn stats stats (non zero cn only)
  DDistr av_in_cn;//for average in area  stats
  double cn_radius;//radius of ROI
  int d_cnr;//radius of ROI+1
  int dpix,dpix_2;//for graphycs, size in pixel of cell
  double graph_mr,graph_md;//for graphycs, max dens and max rotor
  Cells(){}   
  void Init(double s,double dt,double xmin,double xmax,double ymin,double ymax,double T,int r,double cnr,double gmr,double gmd,int dp)//initialises
  {
    step=s;
    step2=step*2;
    double stepq=step*step;
    delta_t=dt;
    delta_t_2=delta_t/2;
    ripet=r;
    T_max=T;
    time_steps=int(T_max/delta_t);
    x_min=xmin;
    x_max=xmax;
    y_min=ymin;
    y_max=ymax;    
    double xL=x_max-x_min;
    double yL=y_max-y_min;
    size_x=int(xL/step);
    size_y=int(yL/step);
    cells=new CN***[ripet];//allocates everything
    in_area=new bool*[size_x];
    for(int ix=0;ix<size_x;ix++)
      {
	in_area[ix]=new bool[size_y];
	for(int iy=0;iy<size_y;iy++) in_area[ix][iy]=false;	
      }
    dpix=dp;
    dpix_2=dpix/2;
    for(int ir=0;ir<ripet;ir++)
      {
	cells[ir]=new CN**[time_steps];
	for(int it=0;it<time_steps;it++)
	  {
	    cells[ir][it]=new CN*[size_x];
	    for(int ix=0;ix<size_x;ix++)
	      {	
		cells[ir][it][ix]=new CN[size_y];
		for(int iy=0;iy<size_y;iy++) cells[ir][it][ix][iy].Init(stepq);//initialises individual cell
	      }
	  }
      }
    dens.Init(ripet,T_max,delta_t);//and inits distributions
    max_cn.Init(ripet,T_max,delta_t);
    av_cn.Init(ripet,T_max,delta_t);
    av_in_cn.Init(ripet,T_max,delta_t);
    cn_radius=cnr;
    d_cnr=int(cn_radius)+1;//r=3.5-> d_cnr=4 for loop on neighs
    graph_mr=gmr;
    graph_md=gmd;    
  }
  bool Check_neighs(int s,int i,int j,int k)//checks if velocity field defined in neighs
  {
    if((cells[s][i][j-1][k].dens>0)&&(cells[s][i][j+1][k].dens>0)&&(cells[s][i][j][k-1].dens>0)&&(cells[s][i][j][k+1].dens>0)) return true;
    else return false;
  }
  bool Check_in(double t,Vector2D r)//checks if inside tracking area and time limit
  {
    if((t<T_max)&&(r.x>=x_min)&&(r.y>=y_min)&&(r.x<x_max)&&(r.y<y_max)) return true;
    else return false;
  }
  void Up_all(int s,double t)//updates number of events in each cell
  {
    if(t<T_max)
      {
	int it=int(t/delta_t);
	for(int j=0;j<size_x;j++)
	  for(int k=0;k<size_y;k++) cells[s][it][j][k].count++;
      }
  }    
  void Update(State2D S,double t,int rip)//updates velocity field value
  {
    if(Check_in(t,S.r))
      {
	int int_t=int(t/delta_t);
	int int_x=int((S.r.x-x_min)/step);
	int int_y=int((S.r.y-y_min)/step);
	cells[rip][int_t][int_x][int_y].v.Add(S.v);//here
	cells[rip][int_t][int_x][int_y].dens+=1;//and counts # of contributions to loc cell
      }
  }
  void Final_dens()//computes average v field and density in central area
  {
    for(int s=0;s<ripet;s++)
      for(int i=0;i<time_steps;i++)
	{
	  for(int j=0;j<size_x;j++)
	    for(int k=0;k<size_y;k++)
	      {
		cells[s][i][j][k].Final();//see above
		if(cells[s][i][j][k].dens>0) in_area[j][k]=true;
	      }
	}
    for(int s=0;s<ripet;s++)
      for(int i=0;i<time_steps;i++)//density is defined as the average over used cells
	{
	  for(int j=0;j<size_x;j++)
	    for(int k=0;k<size_y;k++) if(in_area[j][k]) dens.dd[s][i].Update(cells[s][i][j][k].dens);
	  dens.dd[s][i].Final();
	}
  }
  void Rotor()//computes rotor  FOR THE MOMENT I LEAVE ROTOR AND CN COMPUTATION ALSO TO CELLS NOT IN AREA
  {
    Final_dens();//first density
    for(int s=0;s<ripet;s++)
      {   
	for(int i=0;i<time_steps;i++)
	  for(int j=1;j<size_x-1;j++)
	    for(int k=1;k<size_y-1;k++)
	      {
		if(Check_neighs(s,i,j,k))//checks if defined in neighs
		  {
		    cells[s][i][j][k].rotval=true;
		    cells[s][i][j][k].rot=(cells[s][i][j+1][k].v.y-cells[s][i][j-1][k].v.y-cells[s][i][j][k+1].v.x+cells[s][i][j][k-1].v.x)/step2;//computes
		  }
	      }
      }
  }
  void Final_all()//finalises all distrs
  {
    dens.Final();
    max_cn.Final();
    av_cn.Final();
    av_in_cn.Final();
  }
  void CN_comp()//computing CN
  {
    for(int s=0;s<ripet;s++)
      {    
	for(int i=0;i<time_steps;i++)
	  for(int j=0;j<size_x;j++)
	    for(int k=0;k<size_y;k++)
	      {
		int conta_v=0;//number of non zero vel cells
		double vav=0;//average vel in ROI
		double maxr=-10000000;//inits max,min of rot
		double minr=10000000;
		for(int l=-d_cnr;l<=d_cnr;l++)
		  for(int m=-d_cnr;m<=d_cnr;m++)//scans over neighs
		    {
		      double r=sqrt(l*l+m*m);//euclidian condition
		      if(((j+l)>=0)&&((k+m)>=0)&&((j+l)<size_x)&&((k+m)<size_y)&&(r<=cn_radius))
			{
			  if(cells[s][i][j+l][k+m].rotval)//if rot defined looks for max and min
			    {
			      if(cells[s][i][j+l][k+m].rot>maxr) maxr=cells[s][i][j+l][k+m].rot;
			      if(cells[s][i][j+l][k+m].rot<minr) minr=cells[s][i][j+l][k+m].rot;
			    }
			  if(cells[s][i][j+l][k+m].v.m>0)//if vel defined updates average
			    {
			      conta_v++;
			      vav+=cells[s][i][j+l][k+m].v.m;
			    }
			}
		    }
		if(conta_v)//computes average and cn
		  {
		    vav/=conta_v;
		    if((maxr!=-10000000)&&(minr!=10000000))
		      {
			cells[s][i][j][k].cn=step*(maxr-minr)/(vav*6);
		      }
		    else cells[s][i][j][k].cn=0;//zero if v nowhere or rot nowhere
		  }
		else cells[s][i][j][k].cn=0; 
	      }
      }
    for(int s=0;s<ripet;s++)
      for(int i=0;i<time_steps;i++)//computes max and average cn 
	{
	  double loc_max=0;
	  for(int j=0;j<size_x;j++)
	    for(int k=0;k<size_y;k++)
	      {
		if(cells[s][i][j][k].cn>0)//average on non zero cn
		  {
		    av_cn.dd[s][i].Update(cells[s][i][j][k].cn);
		    if(cells[s][i][j][k].cn>loc_max) loc_max=cells[s][i][j][k].cn;//also checking max
		  }
		if(in_area[j][k]) av_in_cn.dd[s][i].Update(cells[s][i][j][k].cn);
	      }
	  av_cn.dd[s][i].Final();
	  max_cn.dd[s][i].av=loc_max;
	  av_in_cn.dd[s][i].Final();
	}
    Final_all();//finalises
  }
  void Print_v(char name_table[200],char name_graph[200],int r,int t)//prints velocity field
  {
    double vx,vy;
    for(int i=0;i<PIXEL_X;i++)//graph
      for(int j=0;j<PIXEL_Y;j++) DrawPixel(scene,i,j,255,255,255);
    for(int ix=0;ix<size_x;ix++)
      for(int iy=0;iy<size_y;iy++)
	{
	  if(cells[r][t][ix][iy].v.m>0)
	    {
	      if(cells[r][t][ix][iy].v.m<1)
		{
		  vx=cells[r][t][ix][iy].v.x;
		  vy=cells[r][t][ix][iy].v.y;
		}
	      else
		{
		  vx=cells[r][t][ix][iy].v.x/cells[r][t][ix][iy].v.m;
		  vy=cells[r][t][ix][iy].v.y/cells[r][t][ix][iy].v.m;
		}
	      Freccia(2*vx,2*vy,ix*dpix+dpix_2,iy*dpix+dpix_2,dpix);	    
	    }
	}
    for(int i=0;i<size_x;i++)
      for(int j=0;j<size_y;j++)
	{
	  int ix=i*dpix+dpix_2;
	  int iy=j*dpix+dpix_2;
	  if(!in_area[i][j]) DrawBallnp(scene,0,0,0,dpix_2,ix,iy,PIXEL_X,PIXEL_Y);
	}
    DrawIMG(scene,0,0); 
    SDL_Flip(screen);
    SDL_SaveBMP(screen,name_graph);//and table
    ofstream out(name_table);
    for(int i=0;i<size_x;i++)
      {
	for(int j=0;j<size_y;j++)
	  {
	    out << cells[r][t][i][j].v.x << " " << cells[r][t][i][j].v.y << " ";
	  }
	out << endl;
      }
    out.close();
  }
  void Print_rot(char name_table[200],char name_graph[200],int r,int t)//rotor field
  {
    Colore Loccol;
    for(int i=0;i<PIXEL_X;i++)
      for(int j=0;j<PIXEL_Y;j++)
	{
	  int ix=i/dpix;
	  int iy=j/dpix;
	  if(!in_area[ix][iy]) DrawPixel(scene,i,j,255,255,255);
	  if(!cells[r][t][ix][iy].rotval) DrawPixel(scene,i,j,255,255,255);
	  else
	    {
	      Loccol=colorepm0(cells[r][t][ix][iy].rot,graph_mr,-graph_mr);
	      DrawPixel(scene,i,j,Loccol.r,Loccol.g,Loccol.b); 
	    }
	}
    for(int i=0;i<size_x;i++)
      for(int j=0;j<size_y;j++)
	{
	  int ix=i*dpix+dpix_2;
	  int iy=j*dpix+dpix_2;
	  if(!in_area[i][j]) DrawBallnp(scene,0,0,0,dpix_2,ix,iy,PIXEL_X,PIXEL_Y);
	}
    DrawIMG(scene,0,0); 
    SDL_Flip(screen);
    SDL_SaveBMP(screen,name_graph);
    ofstream out(name_table);
    for(int i=0;i<size_x;i++)
      {
	for(int j=0;j<size_y;j++)
	  {
	    out << cells[r][t][i][j].rot << " ";
	  }
	out << endl;
      }
    out.close();
  }
  void Print_cn(char name_table[200],char name_graph[200],int r,int t)//cn field
  {
    Colore Loccol;
    for(int i=0;i<PIXEL_X;i++)
      for(int j=0;j<PIXEL_Y;j++)
	{
	  int ix=i/dpix;
	  int iy=j/dpix;
	  if(!in_area[ix][iy]) DrawPixel(scene,i,j,255,255,255);
	  if(!cells[r][t][ix][iy].cn) DrawPixel(scene,i,j,255,255,255);
	  else
	    {
	      Loccol=colorenuovo(cells[r][t][ix][iy].cn,1,0);
	      DrawPixel(scene,i,j,Loccol.r,Loccol.g,Loccol.b); 
	    }
	}
     for(int i=0;i<size_x;i++)
      for(int j=0;j<size_y;j++)
	{
	  int ix=i*dpix+dpix_2;
	  int iy=j*dpix+dpix_2;
	  if(!in_area[i][j]) DrawBallnp(scene,0,0,0,dpix_2,ix,iy,PIXEL_X,PIXEL_Y);
	}   
    DrawIMG(scene,0,0); 
    SDL_Flip(screen);
    SDL_SaveBMP(screen,name_graph);
    ofstream out(name_table);
    for(int i=0;i<size_x;i++)
      {
	for(int j=0;j<size_y;j++)
	  {
	    out << cells[r][t][i][j].cn << " ";
	  }
	out << endl;
      }
    out.close();
  }
  void Print_dens(char name_table[200],char name_graph[200],int r,int t)//density field
  {                         //plots density as an average over neighbours
    Colore Loccol;
    for(int i=0;i<PIXEL_X;i++)
      for(int j=0;j<PIXEL_Y;j++)
	{
	  int ix=i/dpix;
	  int iy=j/dpix;
	  if(!in_area[ix][iy]) DrawPixel(scene,i,j,255,255,255);
	  else
	    {
	      int lc=1;
	      double lav=0;
	      for(int ii=ix-1;ii<=ix+1;ii++)
		for(int jj=iy-1;jj<=iy+1;jj++)
		  {
		    if((ii>=0)&&(ii<size_x)&&(jj>=0)&&(jj<size_y))
		      {
			lav+=cells[r][t][ii][jj].dens;
			lc++;
		      }
		  }
	      lav/=lc;
	      if(!lav) DrawPixel(scene,i,j,255,255,255);//averages density over neighbours
	      else
		{
		  Loccol=colorenuovo(lav,graph_md,0);
		  DrawPixel(scene,i,j,Loccol.r,Loccol.g,Loccol.b);
		}
	    }
	}
    for(int i=0;i<size_x;i++)
      for(int j=0;j<size_y;j++)
	{
	  int ix=i*dpix+dpix_2;
	  int iy=j*dpix+dpix_2;
	  if(!in_area[i][j]) DrawBallnp(scene,0,0,0,dpix_2,ix,iy,PIXEL_X,PIXEL_Y);
	} 
    DrawIMG(scene,0,0); 
    SDL_Flip(screen);
    SDL_SaveBMP(screen,name_graph);
    ofstream out(name_table);
    for(int i=0;i<size_x;i++)
      {
	for(int j=0;j<size_y;j++)
	  {
	    out << cells[r][t][i][j].dens << " ";
	  }
	out << endl;
      }
    out.close();
  } 
  void Print()//prints everything to files and graphs
  {
    ofstream out;
    char name[200];
    char name_graph[200];
    sprintf(name,"data/R_%f_T_%f_L_%f/dens_time.dat",grid_size,time_int,cn_radius);
    dens.Print(name);
    sprintf(name,"data/R_%f_T_%f_L_%f/max_cn_time.dat",grid_size,time_int,cn_radius);
    max_cn.Print(name);
    sprintf(name,"data/R_%f_T_%f_L_%f/av_cn_time.dat",grid_size,time_int,cn_radius);
    av_cn.Print(name);
    sprintf(name,"data/R_%f_T_%f_L_%f/av_in_cn_time.dat",grid_size,time_int,cn_radius);
    av_in_cn.Print(name);    
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_dens_v.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_dens_v.bmp",grid_size,time_int,cn_radius);    
    Print_v(name,name_graph,dens.imr,dens.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_dens_rot.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_dens_rot.bmp",grid_size,time_int,cn_radius);
    Print_rot(name,name_graph,dens.imr,dens.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_dens_cn.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_dens_cn.bmp",grid_size,time_int,cn_radius);
    Print_cn(name,name_graph,dens.imr,dens.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_dens_dens.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_dens_dens.bmp",grid_size,time_int,cn_radius);
    Print_dens(name,name_graph,dens.imr,dens.imt);    
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_cn_v.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_cn_v.bmp",grid_size,time_int,cn_radius);    
    Print_v(name,name_graph,max_cn.imr,max_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_cn_rot.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_cn_rot.bmp",grid_size,time_int,cn_radius);
    Print_rot(name,name_graph,max_cn.imr,max_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_cn_cn.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_cn_cn.bmp",grid_size,time_int,cn_radius);
    Print_cn(name,name_graph,max_cn.imr,max_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_cn_dens.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_cn_dens.bmp",grid_size,time_int,cn_radius);
    Print_dens(name,name_graph,max_cn.imr,max_cn.imt);        
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_cn_v.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_cn_v.bmp",grid_size,time_int,cn_radius);    
    Print_v(name,name_graph,av_cn.imr,av_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_cn_rot.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_cn_rot.bmp",grid_size,time_int,cn_radius);
    Print_rot(name,name_graph,av_cn.imr,av_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_cn_cn.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_cn_cn.bmp",grid_size,time_int,cn_radius);
    Print_cn(name,name_graph,av_cn.imr,av_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_cn_dens.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_cn_dens.bmp",grid_size,time_int,cn_radius);
    Print_dens(name,name_graph,av_cn.imr,av_cn.imt);   
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_in_cn_v.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_in_cn_v.bmp",grid_size,time_int,cn_radius);    
    Print_v(name,name_graph,av_in_cn.imr,av_in_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_in_cn_rot.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_in_cn_rot.bmp",grid_size,time_int,cn_radius);
    Print_rot(name,name_graph,av_in_cn.imr,av_in_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_in_cn_cn.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_in_cn_cn.bmp",grid_size,time_int,cn_radius);
    Print_cn(name,name_graph,av_in_cn.imr,av_in_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/tables/max_av_in_cn_dens.dat",grid_size,time_int,cn_radius);
    sprintf(name_graph,"data/R_%f_T_%f_L_%f/graphs/max_av_in_cn_dens.bmp",grid_size,time_int,cn_radius);
    Print_dens(name,name_graph,av_in_cn.imr,av_in_cn.imt);
    sprintf(name,"data/R_%f_T_%f_L_%f/final.dat",grid_size,time_int,cn_radius);
    out.open(name);
    out << "Maximum CN=" << max_cn.max << " at t=" << max_cn.imt*delta_t+delta_t_2 << " ripet=" << max_cn.imr << endl;
    out << "Maximum density=" << dens.max << " at t=" << dens.imt*delta_t+delta_t_2 << " ripet=" << dens.imr << endl;
    out << "Maximum average CN=" << av_cn.max << " at t=" << av_cn.imt*delta_t+delta_t_2 << " ripet=" << av_cn.imr << endl;
    out << "Maximum average in area CN=" << av_in_cn.max << " at t=" << av_in_cn.imt*delta_t+delta_t_2 << " ripet=" << av_in_cn.imr << endl;
    out.close();
  }
};




Cells CELLS;


void Run(int sim)//reads data and computes vector field
{
  char nomenome[200];
  sprintf(nomenome,"positions/pos_%d.dat",sim);//reads data
  ifstream pos(nomenome);
  double t;
  int lped;//tracked peds
  while(1)
    {
      pos >> t;//time
      if(pos.eof()) break;
      pos >> lped;//# of pedestrains present at a given time
      CELLS.Up_all(sim,t);
      for(int i=0;i<lped;i++)
	{
	  double x,y,vx,vy;//gets position and velocity
	  pos >> x;
	  pos >> y;	  
	  pos >> vx;
	  pos >> vy;
	  State2D pedestrian(x,y,vx,vy);	  
	  CELLS.Update(pedestrian,t,sim);//updates velocity field
	}
    }
}


void Initialise()//gets environment data and computation parameters
{
  string token;
  ifstream read("parameters");
  read >> token;   
  read >> X_MIN;
  read >> token;   
  read >> Y_MIN;
  read >> token;   
  read >> X_MAX;
  read >> token;   
  read >> Y_MAX;  
  read >> token;    
  read >> ripet;//# of reps in the experiment (position data files)
  read >> token;
  read >> TI;  //total time
  read >> token;
  read >> grid_size;
  read >> token;
  read >> time_int;
  read >> token;
  read >> cn_radius;
  read >> token;
  read >> graph_mr;
  read >> token;
  read >> graph_md;
  read >> token;
  read >> DPIX;  
  read.close();
  XL=X_MAX-X_MIN;
  YL=Y_MAX-Y_MIN;
  double delta=XL-int(XL/grid_size)*grid_size;//it is preferred to have the space defined as multiple of grid size
  if(delta>0)
    {
      XL=int(XL/grid_size)*grid_size+grid_size;
      X_MAX=X_MIN+XL;
    }
  delta=YL-int(YL/grid_size)*grid_size;
  if(delta>0)
    {
      YL=int(YL/grid_size)*grid_size+grid_size;
      Y_MAX=Y_MIN+YL;
    }
  XL_2=XL/2;
  YL_2=YL/2;
  PIXEL_X=int(XL/grid_size)*DPIX;
  PIXEL_Y=int(YL/grid_size)*DPIX;
}

void Names()
{
  char name_check[200];
  sprintf(name_check,"data/R_%f_T_%f_L_%f",grid_size,time_int,cn_radius);
  if (stat(name_check, &st) == -1) {
    mkdir(name_check, 0700);
  }
  sprintf(name_check,"data/R_%f_T_%f_L_%f/graphs",grid_size,time_int,cn_radius);
  if (stat(name_check, &st) == -1) {
    mkdir(name_check, 0700);
  }
  sprintf(name_check,"data/R_%f_T_%f_L_%f/tables",grid_size,time_int,cn_radius);
  if (stat(name_check, &st) == -1) {
    mkdir(name_check, 0700);
  }  
}

int main()
{
  Initialise();//reads parameters
  Init_images();//starts graphics
  SDL_Inizio();//starts graphics
  Names();//creates directories for printing
  CELLS.Init(grid_size,time_int,X_MIN,X_MAX,Y_MIN,Y_MAX,TI,ripet,cn_radius,graph_mr,graph_md,DPIX);//initialises
  for(int sim=0;sim<ripet;sim++)
    {      
      Run(sim);//reads data
    }
  CELLS.Rotor();//computes rotor
  CELLS.CN_comp();//computes cn
  CELLS.Print(); //prints 
  return 0;
}

