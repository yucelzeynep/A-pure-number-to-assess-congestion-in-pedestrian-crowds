#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <SDL/SDL.h>
#include <iostream> 
#include <math.h>  

using namespace std;

extern SDL_Surface *screen;
extern SDL_Surface *scene;
extern int PIXEL_X;
extern int PIXEL_Y;


int min(int i,int j);//minimo
int mod(int i);//valore assoluto

double Angolo(double tx,double ty);
double Angolo_relativo(double a1,double a2);


//Slock va sempre chiamata prima di usare screen, Sulock dopo averlo usato
void Slock(SDL_Surface *screen)
{
  if ( SDL_MUSTLOCK(screen) )
  {
    if ( SDL_LockSurface(screen) < 0 )
    {
      return;
    }
  }
}

void Sulock(SDL_Surface *screen)
{
  if ( SDL_MUSTLOCK(screen) )
  {
    SDL_UnlockSurface(screen);
  }
}

// DrawPixel disegna un pixel in x,y con colore determinato da R,G,B
void DrawPixel(SDL_Surface *screen, int x, int y, Uint8 R, Uint8 G, Uint8 B)
{
  Uint32 color = SDL_MapRGB(screen->format, R, G, B);
  Uint32 *bufp;
  bufp = (Uint32 *)screen->pixels + y*screen->pitch/4 + x;
  *bufp = color;
} 

void SDL_Inizio()   //inizializza con dim X,Y
{
  if ( SDL_Init(SDL_INIT_AUDIO|SDL_INIT_VIDEO) < 0 )
  {
    printf("Unable to init SDL: %s\n", SDL_GetError());
    exit(1);
  }
  atexit(SDL_Quit);
  screen=SDL_SetVideoMode(PIXEL_X,PIXEL_Y,32,SDL_HWSURFACE|SDL_DOUBLEBUF|SDL_RESIZABLE);
  if ( screen == NULL )
  {
    printf("Unable to set video: %s\n", SDL_GetError());
    exit(1);
  }
}  




//Disegna una palla di raggio rg con centro in x1,y1 colore rgb XMAX YMAX massimi valori per i pixel
void DrawBallnp(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 b,int rg,int x1,int y1,int xmax,int ymax)
{
  Slock(screen);
  for (int x=x1-rg;x<x1+rg;x++)
  {
    for (int y=y1-rg;y<y1+rg;y++) 
    {
      if ((((x-x1)*(x-x1)+(y-y1)*(y-y1))<(rg*rg))&&(((x>=0)&&(x<xmax)&&(y>=0)&&(y<ymax)))) DrawPixel(screen,(x+xmax)%xmax,(y+ymax)%ymax,r,g,b);
    }   
  }
  Sulock(screen);
  SDL_Flip(screen);      
}

void DrawBallnp(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 b,int rg,int x1,int y1,int xmax,int ymax,double vx,double vy)
{
  Slock(screen);
  for (int x=x1-rg;x<x1+rg;x++)
  {
    for (int y=y1-rg;y<y1+rg;y++) 
    {
      double sx=x-x1;
      double sy=y1-y;
      double l=sqrt(sx*sx+sy*sy);
      double lv=sqrt(vx*vx+vy*vy);
      bool paint=true;
      if(lv!=0.) 
	{
	  double scal=(sx*vx+sy*vy)/(l*lv);
	  if(scal>sqrt(2.)/2.) paint=false;
	}
      if (paint&&(l<rg)&&(x>=0)&&(x<xmax)&&(y>=0)&&(y<ymax)) DrawPixel(screen,(x+xmax)%xmax,(y+ymax)%ymax,r,g,b);
    }   
  }
  Sulock(screen);
  SDL_Flip(screen);      
}

void DrawBallnp_hollow(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 b,int rg,int x1,int y1,int xmax,int ymax,double vx,double vy)
{
  Slock(screen);
  for (int x=x1-rg;x<x1+rg;x++)
  {
    for (int y=y1-rg;y<y1+rg;y++) 
    {
      double sx=x-x1;
      double sy=y1-y;
      double l=sqrt(sx*sx+sy*sy);
      double lv=sqrt(vx*vx+vy*vy);
      bool paint=true;
      if(lv!=0.) 
	{
	  double scal=(sx*vx+sy*vy)/(l*lv);
	  if(scal>sqrt(2.)/2.) paint=false;
	}
      if (paint&&(l<rg)&&(l>rg/2)&&(x>=0)&&(x<xmax)&&(y>=0)&&(y<ymax)) DrawPixel(screen,(x+xmax)%xmax,(y+ymax)%ymax,r,g,b);
    }   
  }
  Sulock(screen);
  SDL_Flip(screen);      
}


//Disegna una palla di raggio rg con centro in x1,y1 colore rgb XMAX YMAX massimi valori per i pixel
void DrawBallp(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 b,int rg,int x1,int y1,int xmax,int ymax)
{
  Slock(screen);
  for (int x=x1-rg;x<x1+rg;x++)
    {
      for (int y=y1-rg;y<y1+rg;y++) 
	{
	  double dx=abs(x-x1);
	  double dy=abs(y-y1);
	  if((dx*dx+dy*dy)<(rg*rg)) 
	    {
	      int hx,hy;
	      hx=(x+xmax)%xmax;
	      hy=(y+ymax)%ymax;
	      if(hx<0) hx+=xmax;
	      if(hy<0) hy+=ymax;
	      if((hx>=0)&&(hx<PIXEL_X)&&(hy>=0)&&(hy<PIXEL_Y)) DrawPixel(screen,hx,hy,r,g,b);
	    }
	}   
    }
  Sulock(screen);
  SDL_Flip(screen);      
}//usa condizioni periodiche

void DrawEllp(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 bl,double a,double b,double theta,double x0[2],double XMAX,double YMAX,int xmax,int ymax)
{
  double c=a/b;
  double ct=cos(theta);
  double st=sin(theta);
  int x1=int(x0[0]*xmax/(XMAX));
  int y1=int(x0[1]*ymax/(YMAX));
  Slock(screen);
  double ratio=b/XMAX;
  int rg=int(ratio*xmax)+1;
  if(c>1) {rg*=c;rg++;} 
  for (int x=x1-rg;x<x1+rg;x++)
    {
      for (int y=y1-rg;y<y1+rg;y++) 
	{
	    {
	      double dx=x-x1;
	      double dy=y-y1;
	      double cx=dx*ct-dy*st;
	      double cy=dx*st+dy*ct;
	      cx/=c;
	      double l=sqrt(cx*cx+cy*cy);
	      l/=xmax;
	      if(l<ratio)
		{
		  	      int hx,hy;
	      hx=(x+xmax)%xmax;
	      hy=(y+ymax)%ymax;
	      if(hx<0) hx+=xmax;
	      if(hy<0) hy+=ymax;
	      if((hx>=0)&&(hx<PIXEL_X)&&(hy>=0)&&(hy<PIXEL_Y)) DrawPixel(screen,hx,hy,r,g,bl);
		}
	    }	        
	}
    }
  Sulock(screen);
  SDL_Flip(screen);     
}

void DrawEllnp(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 bl,double a,double b,double theta,double x0[2],double XMAX,double YMAX,int xmax,int ymax)
{
  double c=a/b;
  double ct=cos(theta);
  double st=sin(theta);
  int x1=int(x0[0]*xmax/(XMAX));
  int y1=int(x0[1]*ymax/(YMAX));
  Slock(screen);
  double ratio=b/XMAX;
  int rg=int(ratio*xmax)+1;
  if(c>1) {rg*=c;rg++;} 
  for (int x=x1-rg;x<x1+rg;x++)
    {
      for (int y=y1-rg;y<y1+rg;y++) 
	{
	  if ((x>=0)&&(x<xmax)&&(y>=0)&&(y<ymax))
	    {
	      double dx=x-x1;
	      double dy=y-y1;
	      double cx=dx*ct-dy*st;
	      double cy=dx*st+dy*ct;
	      cx/=c;
	      double l=sqrt(cx*cx+cy*cy);
	      l/=xmax;
	      if(l<ratio)
		{
		  DrawPixel(screen,x,y,r,g,bl);
		}
	    }	        
	}
    }
  Sulock(screen);
  SDL_Flip(screen);     
}



//disegna intera img nel punto x,y di screen (per farla tutta usa NULL)
void DrawIMG(SDL_Surface *img, int x, int y)
{
  SDL_Rect dest;
  dest.x = x;
  dest.y = y;
  SDL_BlitSurface(img, NULL, screen, &dest);
}

//disegna img con larghezza w e altezza h a partire da x2,y2 in x,y di screen
void DrawIMG(SDL_Surface *img, int x, int y,int w, int h, int x2, int y2)
{
  SDL_Rect dest;
  dest.x = x;
  dest.y = y;
  SDL_Rect dest2;
  dest2.x = x2;
  dest2.y = y2;
  dest2.w = w;
  dest2.h = h;
  SDL_BlitSurface(img, &dest2, screen, &dest);
} 

//pittura di un colore, x y dimensioni
void PAINT(SDL_Surface *surface,Uint8 r,Uint8 g,Uint8 b,int x,int y)
{
  for (int i=0;i<x;i++) for (int j=0;j<y;j++)
  {
    DrawPixel(surface,i,j,r,g,b);
  }
}


   //aggiunge un pixel
void AddPixel(SDL_Surface *screen,Uint8 r,Uint8 g,Uint8 b,int x,int y,int XX,int YY)
{
  Slock(screen);
  if((x<XX)&&(y<YY)) DrawPixel(screen,x,y,r,g,b);
  Sulock(screen);
  SDL_Flip(screen);      
}

int Init_images() //this is for graphics
{
  Uint32 rmask, gmask, bmask, amask;
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
  rmask = 0xff000000;
  gmask = 0x00ff0000; 
  bmask = 0x0000ff00;
  amask = 0x000000ff;
#else
  rmask = 0x000000ff; 
  gmask = 0x0000ff00;
  bmask = 0x00ff0000;
  amask = 0xff000000;       
#endif 

  scene = SDL_CreateRGBSurface(SDL_SWSURFACE,PIXEL_X,PIXEL_Y,32,
			       rmask, gmask, bmask, amask);
  if(scene == NULL) {fprintf(stderr, "CreateRGBSurface failed: %s\n", SDL_GetError());exit(1);}
   
}  


void Paint(SDL_Surface *surface,int x,int y)
{
  for (int i=0;i<x;i++) 
    for (int j=0;j<y;j++) 
      DrawPixel(surface,i,j,255,255,255);
}

