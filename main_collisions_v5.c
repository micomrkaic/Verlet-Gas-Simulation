#include <SDL.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ── Layout ────────────────────────────────────────────────────────────────────
#define SIM_W       700
#define PANEL_W     230
#define WIN_W       (SIM_W + PANEL_W)
#define WIN_H       720

// ── Physics ───────────────────────────────────────────────────────────────────
#define GRAVITY      500.0
#define RADIUS         5.0
#define DIAM          (2.0 * RADIUS)
#define RESTITUTION    1.0
#define N_PARTICLES  150
#define DT_MAX       (1.0/60.0)

// ── Parabola vessel ───────────────────────────────────────────────────────────
// Bowl: y_surface(x) = PARA_BASE - PARA_A*(x - SIM_W/2)²
// Screen y increases downward; bowl opens upward (acts as floor).
#define PARA_BASE   (WIN_H - 80.0)
#define PARA_A      (0.0030)

static inline double para_y   (double x) {
    double d = x - SIM_W*0.5; return PARA_BASE - PARA_A*d*d;
}
static inline double para_dydx(double x) {
    return -2.0*PARA_A*(x - SIM_W*0.5);
}

// ── Energy / histogram history ────────────────────────────────────────────────
#define ENERGY_HIST  280    // samples in energy ring buffer
#define HIST_BINS     30    // bins in speed-distribution histogram
#define MONO_SPEED   220.0  // initial speed for monoenergetic start (px/s)

typedef struct { double x,y,vx,vy; Uint8 r,g,b; } Particle;
static Particle p[N_PARTICLES];

static double  eng_buf[ENERGY_HIST];  // stores E/E0 (normalised)
static int     eng_head=0, eng_count=0;
static double  eng_max=1.0, eng_min=0.0;
static double  eng_e0=0.0;            // initial energy for normalisation

// Smoothed histogram (exponential moving average for visual stability)
static double  hist_smooth[HIST_BINS];
static double  hist_mb[HIST_BINS];      // theoretical Maxwell-Boltzmann curve
static double  hist_ke_max = 1.0;       // current KE range for histogram x-axis

static long    collision_count = 0;     // total particle-particle collisions

static double randf(void){ return (double)rand()/((double)RAND_MAX+1.0); }

// ── Maxwell-Boltzmann distribution for 2D ideal gas ──────────────────────────
// In 2D, the speed distribution is f(v) = (m/kT)*v*exp(-mv²/2kT)
// In terms of KE = ½mv², f(KE) = (1/kT)*exp(-KE/kT)  (exponential)
// We normalise to [0,1] for display.
static void compute_mb_curve(double mean_ke) {
    if (mean_ke < 1.0) mean_ke = 1.0;
    double kT = mean_ke;   // mean KE per particle = kT in 2D
    // Fill bins: each bin centre is at ke = (b+0.5)/HIST_BINS * hist_ke_max
    double peak = 0.0;
    double raw[HIST_BINS];
    for (int b=0; b<HIST_BINS; b++) {
        double ke = (b+0.5)/HIST_BINS * hist_ke_max;
        raw[b] = exp(-ke/kT);
        if (raw[b] > peak) peak = raw[b];
    }
    for (int b=0; b<HIST_BINS; b++)
        hist_mb[b] = (peak>0) ? raw[b]/peak : 0.0;
}

// ── Particle initialisation ───────────────────────────────────────────────────
static void init_particles(int n, bool monoenergetic) {
    // Place on a rough grid in the upper part of the vessel
    int cols = (int)ceil(sqrt((double)n * SIM_W / (WIN_H*0.5)));
    if (cols < 1) cols = 1;
    double cell_w = (SIM_W - 4*RADIUS) / cols;
    double cell_h = cell_w;
    for (int i=0; i<n; i++) {
        int col = i % cols, row = i / cols;
        p[i].x = 2*RADIUS + col*cell_w + randf()*cell_w*0.5;
        p[i].y = 2*RADIUS + row*cell_h + randf()*cell_h*0.5;
        if (p[i].x > SIM_W-2*RADIUS) p[i].x = SIM_W-2*RADIUS;
        if (p[i].y > WIN_H*0.45)     p[i].y = WIN_H*0.45;

        double speed, angle;
        if (monoenergetic) {
            // All same speed, random direction — perfect for watching thermalisation
            angle = randf()*2.0*M_PI;
            speed = MONO_SPEED;
        } else {
            angle = randf()*2.0*M_PI;
            speed = 80.0 + randf()*300.0;
        }
        p[i].vx = speed*cos(angle);
        p[i].vy = speed*sin(angle);

        // Colour by angle (hue wheel) so you can track mixing
        double hue = angle/(2.0*M_PI)*360.0;
        int hi=(int)(hue/60.0)%6; double f=hue/60.0-floor(hue/60.0);
        Uint8 q=(Uint8)(255*(1-f)), t8=(Uint8)(255*f);
        switch(hi){
            case 0:p[i].r=255;p[i].g=t8; p[i].b=0;  break;
            case 1:p[i].r=q;  p[i].g=255;p[i].b=0;  break;
            case 2:p[i].r=0;  p[i].g=255;p[i].b=t8; break;
            case 3:p[i].r=0;  p[i].g=q;  p[i].b=255;break;
            case 4:p[i].r=t8; p[i].g=0;  p[i].b=255;break;
            default:p[i].r=255;p[i].g=0; p[i].b=q;  break;
        }
    }
    eng_head=eng_count=0; eng_max=1.0; eng_min=0.0; eng_e0=0.0;
    for(int b=0;b<HIST_BINS;b++){hist_smooth[b]=0.0;hist_mb[b]=0.0;}
    hist_ke_max=1.0; collision_count=0;
}

// ── Energy measurement ────────────────────────────────────────────────────────
static double compute_energy(int n, double *ke_out, double *pe_out){
    double ke=0,pe=0;
    for(int i=0;i<n;i++){
        ke+=0.5*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
        pe+=GRAVITY*(WIN_H-RADIUS-p[i].y);
    }
    if (ke_out) *ke_out = ke;
    if (pe_out) *pe_out = pe;
    return ke+pe;
}

static void push_energy(double e){
    if(eng_e0==0.0) eng_e0=(e!=0.0)?e:1.0;   // capture baseline on first push
    double norm=(eng_e0!=0.0)?e/eng_e0:1.0;
    eng_buf[eng_head]=norm; eng_head=(eng_head+1)%ENERGY_HIST;
    if(eng_count<ENERGY_HIST)eng_count++;
    eng_max=eng_min=eng_buf[0];
    for(int i=0;i<eng_count;i++){
        if(eng_buf[i]>eng_max)eng_max=eng_buf[i];
        if(eng_buf[i]<eng_min)eng_min=eng_buf[i];
    }
    // Always include 1.0 in range; keep at least ±0.001 window around it
    if(eng_max<1.0) eng_max=1.0;
    if(eng_min>1.0) eng_min=1.0;
    if(eng_max-eng_min<0.002){ eng_max=1.0+0.001; eng_min=1.0-0.001; }
}

// ── KE histogram update ───────────────────────────────────────────────────────
static void update_histogram(int n){
    // Find max KE across all particles for axis scaling
    double max_ke=1.0;
    for(int i=0;i<n;i++){
        double ke=0.5*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
        if(ke>max_ke)max_ke=ke;
    }
    // Smooth the axis scale so it doesn't jump
    hist_ke_max += (max_ke*1.1 - hist_ke_max)*0.05;

    // Build raw histogram
    double raw[HIST_BINS]={0};
    for(int i=0;i<n;i++){
        double ke=0.5*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
        int b=(int)(ke/hist_ke_max*HIST_BINS);
        if(b>=HIST_BINS)b=HIST_BINS-1;
        raw[b]+=1.0;
    }
    // Normalise raw to [0,1]
    double peak=0;
    for(int b=0;b<HIST_BINS;b++) if(raw[b]>peak)peak=raw[b];
    if(peak>0) for(int b=0;b<HIST_BINS;b++) raw[b]/=peak;
    // Exponential smoothing (α=0.15 → visual stability)
    for(int b=0;b<HIST_BINS;b++)
        hist_smooth[b]=hist_smooth[b]*0.85 + raw[b]*0.15;

    // Recompute MB theoretical curve at current mean KE
    double ke_total=0;
    for(int i=0;i<n;i++) ke_total+=0.5*(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
    compute_mb_curve(n>0 ? ke_total/n : 1.0);
}

// ── CCD: advance one particle under gravity with parabola bounce ──────────────
static void advance_particle(Particle *q, double dt){
    double t_rem=dt; int iters=0;
    while(t_rem>1e-12 && iters++<8){
        double tx  = q->x + q->vx*t_rem;
        double ty  = q->y + q->vy*t_rem + 0.5*GRAVITY*t_rem*t_rem;
        double tvy = q->vy + GRAVITY*t_rem;
        if(ty+RADIUS <= para_y(tx)){ q->x=tx; q->y=ty; q->vy=tvy; return; }
        double lo=0, hi=t_rem;
        for(int k=0;k<48;k++){
            double m=0.5*(lo+hi);
            double mx=q->x+q->vx*m, my=q->y+q->vy*m+0.5*GRAVITY*m*m;
            if(my+RADIUS>para_y(mx)) hi=m; else lo=m;
        }
        double tc=lo;
        q->x+=q->vx*tc; q->y+=q->vy*tc+0.5*GRAVITY*tc*tc; q->vy+=GRAVITY*tc;
        t_rem-=tc;
        double dydx=para_dydx(q->x), nx=-dydx, ny=1.0, len=sqrt(nx*nx+ny*ny);
        nx/=len; ny/=len;
        double vn=q->vx*nx+q->vy*ny;
        if(vn>0){ q->vx=(q->vx-2*vn*nx)*RESTITUTION; q->vy=(q->vy-2*vn*ny)*RESTITUTION; }
        q->y=para_y(q->x)-RADIUS-1e-6;
    }
    if(t_rem>1e-12){
        q->x+=q->vx*t_rem; q->y+=q->vy*t_rem+0.5*GRAVITY*t_rem*t_rem; q->vy+=GRAVITY*t_rem;
    }
}

// ── Event-driven collision times ──────────────────────────────────────────────
// Time until particle i hits the ceiling (y == RADIUS, moving upward)
// Solves: y + vy*t + 0.5*g*t^2 = RADIUS  (quadratic, g>0 so parabola opens down)
static double ceiling_time(int i, double dt){
    if(p[i].vy >= 0) return dt+1.0;   // moving down — can't hit ceiling
    double a=0.5*GRAVITY, b=p[i].vy, c=p[i].y-RADIUS;
    double disc=b*b-4.0*a*c;
    if(disc<0) return dt+1.0;
    double sq=sqrt(disc);
    double t1=(-b-sq)/(2.0*a), t2=(-b+sq)/(2.0*a);
    double best=dt+1.0;
    if(t1>=-1e-10&&t1<=dt) best=t1<0?0:t1;
    if(t2>=-1e-10&&t2<=dt&&t2<best) best=t2<0?0:t2;
    return best;
}

// Time until particle i hits left wall (x == RADIUS, moving left)
static double lwall_time(int i, double dt){
    if(p[i].vx>=0) return dt+1.0;
    double t=(RADIUS-p[i].x)/p[i].vx;
    if(t>=-1e-10&&t<=dt) return t<0?0:t;
    return dt+1.0;
}

// Time until particle i hits right wall (x == SIM_W-RADIUS, moving right)
static double rwall_time(int i, double dt){
    if(p[i].vx<=0) return dt+1.0;
    double t=(SIM_W-RADIUS-p[i].x)/p[i].vx;
    if(t>=-1e-10&&t<=dt) return t<0?0:t;
    return dt+1.0;
}

// Time until particles i,j first contact (gravity cancels in relative frame — exact)
static double contact_time(int i, int j, double dt){
    double dx=p[i].x-p[j].x, dy=p[i].y-p[j].y;
    double dvx=p[i].vx-p[j].vx, dvy=p[i].vy-p[j].vy;
    double a=dvx*dvx+dvy*dvy;
    if(a<1e-14) return dt+1.0;
    double b=2.0*(dx*dvx+dy*dvy);
    double c=dx*dx+dy*dy-DIAM*DIAM;
    if(c<0){ return (b<0)?0.0:dt+1.0; }   // already overlapping + approaching
    double disc=b*b-4.0*a*c;
    if(disc<0) return dt+1.0;
    double t=(-b-sqrt(disc))/(2.0*a);
    if(t>=-1e-10&&t<=dt) return t<0?0:t;
    return dt+1.0;
}

// ── Full simulation step: event-driven ───────────────────────────────────────
// All events (ceiling, walls, particle-particle) compete for earliest time.
// Advance all particles to that time, process the event, repeat.
// This guarantees exact collision positions regardless of timestep size.
#define EV_PP    0
#define EV_CEIL  1
#define EV_LWALL 2
#define EV_RWALL 3

static void step_sim(int n, double dt){
    double rim_y=para_y(RADIUS);
    double t_done=0;
    int max_events=n*12;    // safety cap — typically ~N events per frame

    while(t_done<dt-1e-12 && max_events-->0){
        double t_rem=dt-t_done;
        double t_min=t_rem;
        int ev_type=-1, ei=-1, ej=-1;

        // Ceiling and wall events
        for(int i=0;i<n;i++){
            double tc=ceiling_time(i,t_rem);
            if(tc<t_min){ t_min=tc; ev_type=EV_CEIL; ei=i; }
            if(p[i].y < rim_y){   // side walls only above parabola rim
                tc=lwall_time(i,t_rem);
                if(tc<t_min){ t_min=tc; ev_type=EV_LWALL; ei=i; }
                tc=rwall_time(i,t_rem);
                if(tc<t_min){ t_min=tc; ev_type=EV_RWALL; ei=i; }
            }
        }
        // Particle-particle contact events
        for(int i=0;i<n-1;i++){
            for(int j=i+1;j<n;j++){
                double tc=contact_time(i,j,t_rem);
                if(tc<t_min){ t_min=tc; ev_type=EV_PP; ei=i; ej=j; }
            }
        }

        // Advance all particles to earliest event time
        for(int i=0;i<n;i++) advance_particle(&p[i],t_min);
        t_done+=t_min;

        // Process event
        if(ev_type==EV_CEIL){
            p[ei].y=RADIUS;
            p[ei].vy=fabs(p[ei].vy)*RESTITUTION;
        } else if(ev_type==EV_LWALL){
            p[ei].x=RADIUS;
            p[ei].vx=fabs(p[ei].vx)*RESTITUTION;
        } else if(ev_type==EV_RWALL){
            p[ei].x=SIM_W-RADIUS;
            p[ei].vx=-fabs(p[ei].vx)*RESTITUTION;
        } else if(ev_type==EV_PP){
            double dx=p[ej].x-p[ei].x, dy=p[ej].y-p[ei].y;
            double dist=sqrt(dx*dx+dy*dy);
            if(dist<1e-12) dist=DIAM;
            double nx=dx/dist, ny=dy/dist;
            double dvn=(p[ej].vx-p[ei].vx)*nx+(p[ej].vy-p[ei].vy)*ny;
            if(dvn<0){
                p[ei].vx+=dvn*nx; p[ei].vy+=dvn*ny;
                p[ej].vx-=dvn*nx; p[ej].vy-=dvn*ny;
                collision_count++;
            }
            // Tiny separation to prevent immediate re-detection
            p[ei].x-=1e-6*nx; p[ei].y-=1e-6*ny;
            p[ej].x+=1e-6*nx; p[ej].y+=1e-6*ny;
        }
    }
    (void)rim_y;
}

// ── Drawing ───────────────────────────────────────────────────────────────────
static void draw_particle(SDL_Renderer *rend, double x, double y,
                           Uint8 r, Uint8 g, Uint8 b){
    int cx=(int)x,cy=(int)y,R=(int)RADIUS;
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
    for(int dy2=-R-3;dy2<=R+3;dy2++) for(int dx2=-R-3;dx2<=R+3;dx2++){
        double d=sqrt((double)(dx2*dx2+dy2*dy2));
        if(d<R+3&&d>=R){double a=1.0-(d-R)/3.0;
            SDL_SetRenderDrawColor(rend,r,g,b,(Uint8)(a*70));
            SDL_RenderDrawPoint(rend,cx+dx2,cy+dy2);}
    }
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);
    for(int dy2=-R;dy2<=R;dy2++) for(int dx2=-R;dx2<=R;dx2++)
        if(dx2*dx2+dy2*dy2<=R*R){
            double br=1.0-0.3*(dx2*dx2+dy2*dy2)/(double)(R*R);
            SDL_SetRenderDrawColor(rend,(Uint8)(r*br),(Uint8)(g*br),(Uint8)(b*br),255);
            SDL_RenderDrawPoint(rend,cx+dx2,cy+dy2);}
}

static void draw_vessel(SDL_Renderer *rend){
    SDL_SetRenderDrawColor(rend,0x55,0x77,0xAA,0xFF);
    for(int xs=0;xs<SIM_W-1;xs++){
        int y0=(int)para_y(xs),y1=(int)para_y(xs+1);
        for(int t=-1;t<=1;t++) SDL_RenderDrawLine(rend,xs,y0+t,xs+1,y1+t);
    }
    double rim_y=para_y(RADIUS);
    SDL_Rect lw={0,0,2,(int)rim_y},rw={SIM_W-2,0,2,(int)rim_y},cl={0,0,SIM_W,2};
    SDL_RenderFillRect(rend,&lw); SDL_RenderFillRect(rend,&rw);
    SDL_RenderFillRect(rend,&cl);
}

// ── Tiny 3×5 pixel font — digits 0-9, uppercase A-Z, '.', '-', '/', '(', ')' ──
// Each glyph is 5 rows of 3 bits (MSB = left column).
// Index: 0-9 = digits, 10-35 = A-Z, 36='.', 37='-', 38='/', 39='(', 40=')', 41=' '
static const Uint8 F[][5]={
    // 0-9
    {0x7,0x5,0x5,0x5,0x7},{0x2,0x6,0x2,0x2,0x7},{0x7,0x1,0x7,0x4,0x7},
    {0x7,0x1,0x7,0x1,0x7},{0x5,0x5,0x7,0x1,0x1},{0x7,0x4,0x7,0x1,0x7},
    {0x7,0x4,0x7,0x5,0x7},{0x7,0x1,0x1,0x1,0x1},{0x7,0x5,0x7,0x5,0x7},
    {0x7,0x5,0x7,0x1,0x7},
    // A-Z (index 10-35)
    {0x2,0x5,0x7,0x5,0x5}, // A
    {0x6,0x5,0x6,0x5,0x6}, // B
    {0x7,0x4,0x4,0x4,0x7}, // C
    {0x6,0x5,0x5,0x5,0x6}, // D
    {0x7,0x4,0x6,0x4,0x7}, // E
    {0x7,0x4,0x6,0x4,0x4}, // F
    {0x7,0x4,0x5,0x5,0x7}, // G
    {0x5,0x5,0x7,0x5,0x5}, // H
    {0x7,0x2,0x2,0x2,0x7}, // I
    {0x1,0x1,0x1,0x5,0x7}, // J
    {0x5,0x5,0x6,0x5,0x5}, // K
    {0x4,0x4,0x4,0x4,0x7}, // L
    {0x5,0x7,0x7,0x5,0x5}, // M
    {0x5,0x7,0x5,0x5,0x5}, // N
    {0x2,0x5,0x5,0x5,0x2}, // O
    {0x6,0x5,0x6,0x4,0x4}, // P
    {0x2,0x5,0x5,0x7,0x3}, // Q
    {0x6,0x5,0x6,0x5,0x5}, // R
    {0x7,0x4,0x7,0x1,0x7}, // S
    {0x7,0x2,0x2,0x2,0x2}, // T
    {0x5,0x5,0x5,0x5,0x7}, // U
    {0x5,0x5,0x5,0x5,0x2}, // V
    {0x5,0x5,0x7,0x7,0x5}, // W
    {0x5,0x5,0x2,0x5,0x5}, // X
    {0x5,0x5,0x2,0x2,0x2}, // Y
    {0x7,0x1,0x2,0x4,0x7}, // Z
    // punctuation (index 36+)
    {0x0,0x0,0x0,0x0,0x2}, // 36 '.'
    {0x0,0x0,0x7,0x0,0x0}, // 37 '-'
    {0x1,0x1,0x2,0x4,0x4}, // 38 '/'
    {0x3,0x4,0x4,0x4,0x3}, // 39 '('
    {0x6,0x1,0x1,0x1,0x6}, // 40 ')'
    {0x0,0x0,0x0,0x0,0x0}, // 41 ' '
};
static int gl(char c){
    if (c>='0' && c<='9') return c-'0';
    if (c>='A' && c<='Z') return 10+(c-'A');
    if (c>='a' && c<='z') return 10+(c-'a');  // lowercase maps to uppercase
    if (c=='.') return 36;
    if (c=='-') return 37;
    if (c=='/') return 38;
    if (c=='(') return 39;
    if (c==')') return 40;
    return 41; // space / unknown
}
static void ds(SDL_Renderer *rend,const char *s,int x,int y,int sc,
               Uint8 r,Uint8 g,Uint8 b){
    SDL_SetRenderDrawColor(rend,r,g,b,255);
    for(int ci=0;s[ci];ci++) for(int row=0;row<5;row++) for(int col=0;col<3;col++)
        if(F[gl(s[ci])][row]&(0x4>>col)){
            SDL_Rect px={x+ci*(3*sc+sc)+col*sc,y+row*sc,sc,sc};
            SDL_RenderFillRect(rend,&px);}
}

// Draw one stat row: label in light colour at scale 2, value in bright colour at scale 3.
// Returns new y after the row.
static int draw_stat(SDL_Renderer *rend, int px, int y,
                     const char *label, const char *val,
                     Uint8 vr, Uint8 vg, Uint8 vb){
    ds(rend, label, px, y,    2, 0xAA, 0xBB, 0xCC);
    ds(rend, val,   px, y+12, 3, vr,   vg,   vb);
    return y + 30;
}

// ── Stats panel ───────────────────────────────────────────────────────────────
static void draw_panel(SDL_Renderer *rend,int n,
                        double total_e,double ke,double pe){
    SDL_SetRenderDrawColor(rend,0x0C,0x0F,0x1C,0xFF);
    SDL_Rect bg={SIM_W,0,PANEL_W,WIN_H}; SDL_RenderFillRect(rend,&bg);
    SDL_SetRenderDrawColor(rend,0x44,0x55,0x77,0xFF);
    SDL_RenderDrawLine(rend,SIM_W,0,SIM_W,WIN_H);

    int px=SIM_W+8, ty=10;
    char buf[32];

    snprintf(buf,sizeof(buf),"%d",n);
    ty=draw_stat(rend,px,ty,"PARTICLES",buf,0xFF,0xFF,0xFF); ty+=4;

    snprintf(buf,sizeof(buf),"%.4f",(eng_e0>0)?total_e/eng_e0:1.0);
    ty=draw_stat(rend,px,ty,"TOTAL ENERGY (E/E0)",buf,0xFF,0xCC,0x44); ty+=4;

    snprintf(buf,sizeof(buf),"%.4f",(total_e!=0)?ke/total_e:0.0);
    ty=draw_stat(rend,px,ty,"KINETIC (KE/E)",buf,0x44,0xFF,0xAA); ty+=4;

    snprintf(buf,sizeof(buf),"%.4f",(total_e!=0)?pe/total_e:0.0);
    ty=draw_stat(rend,px,ty,"POTENTIAL (PE/E)",buf,0x88,0xAA,0xFF); ty+=4;

    double mean_ke=(n>0)?ke/n:0.0;
    snprintf(buf,sizeof(buf),"%.1f",mean_ke);
    ty=draw_stat(rend,px,ty,"TEMPERATURE (KE/N)",buf,0xFF,0x66,0x44); ty+=4;

    snprintf(buf,sizeof(buf),"%ld",collision_count);
    ty=draw_stat(rend,px,ty,"COLLISIONS",buf,0xCC,0x88,0xFF); ty+=10;

    // ── Energy trace ──────────────────────────────────────────────────────────
    ds(rend,"ENERGY / E0",px,ty,2,0xAA,0xBB,0xCC); ty+=13;
    int gx=SIM_W+4,gw=PANEL_W-8,gh=80,gy=ty;
    SDL_SetRenderDrawColor(rend,0x06,0x08,0x10,0xFF);
    SDL_Rect gb={gx,gy,gw,gh}; SDL_RenderFillRect(rend,&gb);

    double range=eng_max-eng_min;
    // y pixel for a normalised value v
    #define TRACE_Y(v) (gy + gh - 1 - (int)(((v)-eng_min)/range*(gh-2)))

    // Dim grid lines at 25% intervals
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend,0x22,0x33,0x55,80);
    for(int f=1;f<=3;f++) SDL_RenderDrawLine(rend,gx,gy+gh-gh*f/4,gx+gw,gy+gh-gh*f/4);
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);

    // Reference line at E/E0 = 1.0
    int ref_y=TRACE_Y(1.0);
    if(ref_y>=gy && ref_y<gy+gh){
        SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(rend,0xFF,0xFF,0xFF,60);
        SDL_RenderDrawLine(rend,gx,ref_y,gx+gw,ref_y);
        SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);
    }

    // Trace bars
    if(eng_count>1){
        int cols=(eng_count<gw)?eng_count:gw;
        for(int col=0;col<cols;col++){
            int idx=(eng_head-1-col+ENERGY_HIST)%ENERGY_HIST;
            double v=eng_buf[idx];
            double frac=(range>0)?(v-eng_min)/range:0.5;
            int bar_y=TRACE_Y(v);
            if(bar_y<gy) bar_y=gy;
            if(bar_y>=gy+gh) bar_y=gy+gh-1;
            int bx=gx+gw-1-col;
            // colour: green when near 1.0, red when drifted
            double dev=fabs(v-1.0)/(range*0.5+1e-9);
            if(dev>1.0)dev=1.0;
            Uint8 cr=(Uint8)(dev*255);
            Uint8 cg=(Uint8)((1.0-dev)*200);
            SDL_SetRenderDrawColor(rend,cr,cg,0x18,0xFF);
            SDL_RenderDrawLine(rend,bx,bar_y,bx,gy+gh-1);
            (void)frac;
        }
    }
    SDL_SetRenderDrawColor(rend,0x2A,0x3A,0x55,0xFF);
    SDL_RenderDrawRect(rend,&gb);

    // Current deviation label below graph
    if(eng_count>0){
        double cur=eng_buf[(eng_head-1+ENERGY_HIST)%ENERGY_HIST];
        double pct=(cur-1.0)*100.0;
        char devbuf[16];
        // Format as e.g. "+0.02%" or "-1.34%"
        snprintf(devbuf,sizeof(devbuf),"%+.2f PCT",pct);
        ds(rend,"DRIFT",px,ty+gh+2,2,0xAA,0xBB,0xCC);
        Uint8 dcr=(Uint8)(fabs(pct)>1.0?220:fabs(pct)>0.1?180:100);
        Uint8 dcg=(Uint8)(fabs(pct)>1.0?80 :fabs(pct)>0.1?180:220);
        ds(rend,devbuf,px+40,ty+gh+2,2,dcr,dcg,0x44);
    }
    #undef TRACE_Y
    ty+=gh+18;

    // ── KE distribution histogram + Maxwell-Boltzmann curve ───────────────────
    ds(rend,"KE DISTRIBUTION vs MB",px,ty,2,0xAA,0xBB,0xCC); ty+=13;

    int hx=SIM_W+4, hw=PANEL_W-8, hh=120, hy=ty;
    SDL_SetRenderDrawColor(rend,0x06,0x08,0x10,0xFF);
    SDL_Rect hb={hx,hy,hw,hh}; SDL_RenderFillRect(rend,&hb);

    // Grid lines
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend,0x22,0x33,0x55,80);
    for(int f=1;f<=3;f++) SDL_RenderDrawLine(rend,hx,hy+hh-hh*f/4,hx+hw,hy+hh-hh*f/4);
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);

    double bin_w=(double)hw/HIST_BINS;
    for(int b=0;b<HIST_BINS;b++){
        int bx=(int)(hx+b*bin_w);
        int bw=(int)bin_w-1;
        if (bw < 1) bw = 1;

        // Measured histogram — cyan bars
        int bar_h=(int)(hist_smooth[b]*(hh-2));
        if(bar_h>0){
            SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(rend,0x22,0xAA,0xCC,180);
            SDL_Rect br={bx,hy+hh-bar_h,bw,bar_h};
            SDL_RenderFillRect(rend,&br);
            SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);
        }

        // MB theoretical curve — bright yellow line connecting bin tops
        int mb_y=(int)(hy+hh-hist_mb[b]*(hh-2));
        if(b>0){
            int prev_mb_y=(int)(hy+hh-hist_mb[b-1]*(hh-2));
            int bx_prev=(int)(hx+(b-0.5)*bin_w);
            int bx_cur =(int)(hx+(b+0.5)*bin_w);
            SDL_SetRenderDrawColor(rend,0xFF,0xDD,0x33,0xFF);
            SDL_RenderDrawLine(rend,bx_prev,prev_mb_y,bx_cur,mb_y);
        }
        (void)mb_y;
    }
    // Draw MB curve endpoint
    {
        int prev_y=(int)(hy+hh-hist_mb[HIST_BINS-2]*(hh-2));
        int last_y=(int)(hy+hh-hist_mb[HIST_BINS-1]*(hh-2));
        SDL_SetRenderDrawColor(rend,0xFF,0xDD,0x33,0xFF);
        SDL_RenderDrawLine(rend,(int)(hx+(HIST_BINS-1.5)*bin_w),prev_y,
                                (int)(hx+(HIST_BINS-0.5)*bin_w),last_y);
    }

    SDL_SetRenderDrawColor(rend,0x2A,0x3A,0x55,0xFF);
    SDL_RenderDrawRect(rend,&hb);
    ty+=hh+6;

    // Legend
    SDL_SetRenderDrawColor(rend,0x22,0xAA,0xCC,0xFF);
    SDL_Rect lc={px,ty+2,8,8}; SDL_RenderFillRect(rend,&lc);
    ds(rend,"MEASURED",px+11,ty,2,0x44,0xCC,0xEE);
    ty+=14;
    SDL_SetRenderDrawColor(rend,0xFF,0xDD,0x33,0xFF);
    SDL_RenderDrawLine(rend,px,ty+4,px+8,ty+4);
    ds(rend,"MAXWELL-BOLTZMANN",px+11,ty,2,0xDD,0xBB,0x44);
    ty+=16;

    // Controls
    ds(rend,"R-RANDOM RESET", px,WIN_H-52,2,0x77,0x88,0xAA);
    ds(rend,"M-MONO  Q-QUIT", px,WIN_H-30,2,0x77,0x88,0xAA);
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main(int argc,char **argv){
    (void)argc;(void)argv;
    srand((unsigned)time(NULL));

    if(SDL_Init(SDL_INIT_VIDEO)!=0){fprintf(stderr,"SDL: %s\n",SDL_GetError());return 1;}
    SDL_Window *win=SDL_CreateWindow("Hard Sphere Gas Inside a Parabolic Vessel",
        SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,WIN_W,WIN_H,SDL_WINDOW_SHOWN);
    SDL_Renderer *rend=SDL_CreateRenderer(win,-1,
        SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);
    if(!win||!rend){fprintf(stderr,"SDL: %s\n",SDL_GetError());return 1;}

    int n=N_PARTICLES;
    init_particles(n,true);   // start monoenergetic to show thermalisation

    Uint32 t_prev=SDL_GetTicks(), t_last_samp=t_prev;
    bool running=true;

    while(running){
        SDL_Event ev;
        while(SDL_PollEvent(&ev)){
            if(ev.type==SDL_QUIT)running=false;
            if(ev.type==SDL_KEYDOWN) switch(ev.key.keysym.sym){
                case SDLK_ESCAPE:case SDLK_q: running=false; break;
                case SDLK_r: init_particles(n,false); break;  // random speeds
                case SDLK_m: init_particles(n,true);  break;  // monoenergetic
                case SDLK_EQUALS:case SDLK_PLUS:
                    if (n < N_PARTICLES) { n++; } break;
                case SDLK_MINUS:
                    if (n > 2) { n--; } break;
            }
        }

        Uint32 t_now=SDL_GetTicks();
        double dt=(t_now-t_prev)/1000.0;
        if(dt>DT_MAX)dt=DT_MAX;
        t_prev=t_now;

        step_sim(n,dt);

        if(t_now-t_last_samp>=80){
            double ke,pe;
            push_energy(compute_energy(n,&ke,&pe));
            update_histogram(n);
            t_last_samp=t_now;
        }

        SDL_SetRenderDrawColor(rend,0x0D,0x0F,0x1A,0xFF);
        SDL_RenderClear(rend);
        draw_vessel(rend);
        for(int i=0;i<n;i++) draw_particle(rend,p[i].x,p[i].y,p[i].r,p[i].g,p[i].b);
        double ke,pe;
        double te=compute_energy(n,&ke,&pe);
        draw_panel(rend,n,te,ke,pe);
        SDL_RenderPresent(rend);
    }

    SDL_DestroyRenderer(rend); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
