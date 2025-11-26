// compile command; gcc multibody.c -o multibodygrav.exe -I "C:\Users\Liamm\OneDrive\Desktop\C Files\SDL3\x86_64-w64-mingw32\include" -L "C:\Users\Liamm\OneDrive\Desktop\C Files\SDL3\x86_64-w64-mingw32\lib" -lSDL3

#include <SDL3/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define WINDOW_W 1000
#define WINDOW_H 800
#define MAX_TRAIL_POINTS 1000
#define MAX_BALLS 100



// Phys constants
const double ACCELDUETOGRAVITY = 0.0; // can give it downward gravity if so desired
const double FOFGRAV = 50.0f;
const double RESTITUTION = 0.8;
const double FRICTION_X = 0.995;
const double DAMPING_COEF = 0.9f;
const double FRICTION_BETW = 1.0f;

double hue = 0.0;
int ballCount = 0;



// Ball Struct
typedef struct {
    float x, y;
    float vx, vy;
    float ax, ay;
    float radius;
    float mass;
    float trail_x[MAX_TRAIL_POINTS];
    float trail_y[MAX_TRAIL_POINTS];
    int count;

    Uint8 r, g, b; // colour
} Ball;

// array of balls
Ball balls[MAX_BALLS];

Ball makeBall(){
    Ball b;

    // Position
    b.x = (rand() % WINDOW_W - 0 + 1);
    b.y = (rand() % WINDOW_H - 0 + 1);

    // Velocity initially zero
    b.vx = 0.0f;
    b.vy = 0.0f;

    // Size and mass
    b.radius = 20;
    b.mass = 100000;     // simple mass model

    // Random color
    b.r = rand() % 256;
    b.g = rand() % 256;
    b.b = rand() % 256;

    // Trail values
    float trail_x[MAX_TRAIL_POINTS] = {0};
    float trail_y[MAX_TRAIL_POINTS] = {0};
    int count;

    int p_selected_ball;

    return b;

}

// Trail Point Array
void trail_add_point(Ball *t) {
    float x = t->x;
    float y = t->y;
    if (t->count >= MAX_TRAIL_POINTS) {
        // shift everything left (like a moving window)
        memmove(&t->trail_x[0], &t->trail_x[1], sizeof(float) * (MAX_TRAIL_POINTS - 1));
        memmove(&t->trail_y[0], &t->trail_y[1], sizeof(float) * (MAX_TRAIL_POINTS - 1));
        t->count = MAX_TRAIL_POINTS - 1;
    }

    t->trail_x[t->count] = x;
    t->trail_y[t->count] = y;
    t->count++;
}


// Trail Renderer
void trail_render(SDL_Renderer *ren, Ball *t, Uint8 r, Uint8 g, Uint8 b)
{
    SDL_SetRenderDrawColor(ren, r, g, b, 255);

    for (int i = 3; i < t->count; i++) {
        SDL_RenderLine(ren,(int)t->trail_x[i - 3], (int)t->trail_y[i - 3], (int)t->trail_x[i], (int)t->trail_y[i]);
    }
}


// SDL_RenderLine(renderer, x1, y1, x2, y2)
void draw_filled_circle(SDL_Renderer *r, Ball *b) 
{
    for (int dy = -b->radius; dy <= b->radius; dy++) {
        int y = b->y + dy;
        double dx = sqrt((double)(b->radius * b->radius) - (double)(dy * dy));
        int x1 = b->x - (int)dx;
        int x2 = b->x + (int)dx;
        SDL_RenderLine(r, x1, y, x2, y);
    }
}

// yeah i steal code how did you know
// HSV to RGB converter
void HSVtoRGB(float h, float s, float v, Uint8* r, Uint8* g, Uint8* b) {
    float c = v * s;
    float x = c * (1 - fabsf(fmodf(h / 60.0f, 2) - 1));
    float m = v - c;

    float r1, g1, b1;
    if (h < 60)      { r1 = c; g1 = x; b1 = 0; }
    else if (h < 120){ r1 = x; g1 = c; b1 = 0; }
    else if (h < 180){ r1 = 0; g1 = c; b1 = x; }
    else if (h < 240){ r1 = 0; g1 = x; b1 = c; }
    else if (h < 300){ r1 = x; g1 = 0; b1 = c; }
    else             { r1 = c; g1 = 0; b1 = x; }

    *r = (Uint8)((r1 + m) * 255);
    *g = (Uint8)((g1 + m) * 255);
    *b = (Uint8)((b1 + m) * 255);
}

// force of gravity on a by b 
void applyGravity(Ball *a, Ball *b, float deltaTime)
{
    float dx = b->x - a->x;
    float dy = b->y - a->y;

    float distSq = dx*dx + dy*dy;
    if (distSq < 1.0f) distSq = 1.0f; // no div / 0 :(

    float dist = SDL_sqrtf(distSq);

    // Normal direction
    float nx = dx / dist;
    float ny = dy / dist;

    // Force magnitude
    float force = FOFGRAV * (a->mass * b->mass) / distSq;

    // Accel on ball
    float ax = nx * (force / a->mass);
    float ay = ny * (force / a->mass);

    // Apply acceleration to velocity
    a->vx += ax * deltaTime;
    a->vy += ay * deltaTime;
}

void handle_collisions(int a) {
    int i = a;
    for (int j = i + 1; j < ballCount; j++) {

        Ball *a = &balls[i];
        Ball *b = &balls[j];

        float dx = b->x - a->x;
        float dy = b->y - a->y;
        float dist = sqrtf(dx*dx + dy*dy);
        float minDist = a->radius + b->radius;

        // Are they touching?
        if (dist < minDist && dist > 0.0f) {

            // Push them apart (fix overlap)
            float overlap = minDist - dist;
            float nx = dx / dist;
            float ny = dy / dist;

            a->x -= nx * overlap * 0.5f;
            a->y -= ny * overlap * 0.5f;
            b->x += nx * overlap * 0.5f;
            b->y += ny * overlap * 0.5f;

            // Velocities: Separate tangent & normal
            float tx = -ny;
            float ty = nx;

            float vAn = a->vx * nx + a->vy * ny;
            float vBn = b->vx * nx + b->vy * ny;

            float vAt = a->vx * tx + a->vy * ty;
            float vBt = b->vx * tx + b->vy * ty;

            // 1D elastic collision along normal
            float m1 = (vAn * (a->mass - b->mass) + 2.0f * b->mass * vBn) / (a->mass + b->mass);
            float m2 = (vBn * (b->mass - a->mass) + 2.0f * a->mass * vAn) / (a->mass + b->mass);

            // Convert back to X/Y
            a->vx = (tx * vAt + nx * m1) * FRICTION_BETW;
            a->vy = (ty * vAt + ny * m1) * FRICTION_BETW;
            b->vx = (tx * vBt + nx * m2) * FRICTION_BETW;
            b->vy = (ty * vBt + ny * m2) * FRICTION_BETW;
        }
    }
}


void updatephysics(Ball *b, int i, double dt){
    b->x += b->vx * dt;
    b->y += b->vy * dt;

}

void applyGravityBetween(Ball *a, Ball *b, double dt)
{
    const float G = 500.0f; // tune for effect

    float dx = b->x - a->x;
    float dy = b->y - a->y;
    float r2 = dx*dx + dy*dy;

    if (r2 < 1.0f) return; // prevent division by zero

    float r = sqrtf(r2);

    float force = G * a->mass * b->mass / r2;

    float ax = force * dx / r / a->mass;
    float ay = force * dy / r / a->mass;

    float bx = -force * dx / r / b->mass;
    float by = -force * dy / r / b->mass;

    // Apply accelerations
    a->vx += ax * dt;
    a->vy += ay * dt;

    b->vx += bx * dt;
    b->vy += by * dt;
}

int is_point_in_circle(Ball *b, float m_x, float m_y){
    return sqrt((b->x - m_x)*(b->x - m_x) + (b->y - m_y)*(b->y - m_y)) < (b->radius * b->radius);
}

void collideWithWalls(Ball *b, float width, float height, float restitution, float frictionX) {
    // Floor
    if (b->y > height - b->radius) {
        b->y = height - b->radius;
        b->vy = -b->vy * restitution;
        b->vx *= frictionX;
    }

    // Ceiling
    if (b->y < b->radius) {
        b->y = b->radius;
        b->vy = -b->vy * restitution;
    }

    // Left wall
    if (b->x < b->radius) {
        b->x = b->radius;
        b->vx = -b->vx * restitution;
    }

    // Right wall
    if (b->x > width - b->radius) {
        b->x = width - b->radius;
        b->vx = -b->vx * restitution;
    }
}

int main(void)
{
    // Intialize SDL
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }

    // SDL3 CreateWindow(title, width, height, flags)
    SDL_Window *win = SDL_CreateWindow(
        "SDL3 Ball Gravity",
        WINDOW_W, WINDOW_H,
        SDL_WINDOW_RESIZABLE
    );

    // Did window open correctly?
    if (!win) {
        SDL_Log("Window error: %s", SDL_GetError());
        return 1;
    }

    // SDL3 CreateRenderer(window, name)
    SDL_Renderer *ren = SDL_CreateRenderer(win, NULL);
    if (!ren) {
        SDL_Log("Renderer error: %s", SDL_GetError());
        return 1;
    }

    // turns Vsync on, Framerate = Refresh rate (1 to 1)
    SDL_SetRenderVSync(ren, 1);

    // Create balls
    balls[0] = makeBall();
    balls[1] = makeBall();
    balls[2] = makeBall();
    ballCount = 3;
    

    Uint64 prev = SDL_GetTicks();   // 64 bit unsigned int to store number of ticks since library initialization, outside loop for retention

    float lastMouseX, lastMouseY;   // initialize mouse position vars outside for storage between loops
    Ball *p_selected_ball = NULL;
    bool running = true;
    while (running) {

        // TIMING
        Uint64 now = SDL_GetTicks();
        double dt = (now - prev) / 1000.0f;     // milliseconds to seconds
        if (dt <= 0.0) dt = 0.001f;             // no div / 0  :(
        if (dt > 0.1f) dt = 0.1f;               // no huge jumps in time
        prev = now;                             // shift by one for next cycle

        // Keypresses n shieet
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            
            if (e.type == SDL_EVENT_QUIT)
            running = false;

            float instantVx;
            float instantVy;
            float smoothFactor;
            float mouseX, mouseY;
            SDL_MouseButtonFlags buttonState = SDL_GetMouseState(&mouseX, &mouseY);
            
            // mouse stuff - can drag the ball, ball inherits mouses velocity  
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN){
                p_selected_ball = NULL;
                for(int i = 0; i < ballCount; i++) {
                    if(is_point_in_circle){
                        p_selected_ball = &balls[i];
                    }
                }

            p_selected_ball->x = mouseX;
            p_selected_ball->y = mouseY;

            p_selected_ball->vx = 0.0;
            p_selected_ball->vy = 0.0;

            lastMouseX = mouseX;
            lastMouseY = mouseY;
            
            }

            if (e.type == SDL_EVENT_KEY_DOWN) { // can do a little controlling
                SDL_Keycode key = e.key.key;
                
                if (key == SDLK_ESCAPE) running = false;
                if (key == SDLK_LEFT)  p_selected_ball->vx -= 50;
                if (key == SDLK_RIGHT) p_selected_ball->vx += 50;
                if (key == SDLK_UP) p_selected_ball->vy += -50;
                if (key == SDLK_DOWN) p_selected_ball->vy += 50;
            } 
        
        
            if (e.type == SDL_EVENT_MOUSE_BUTTON_UP){
                instantVx = (mouseX - lastMouseX) / dt;
                instantVy = (mouseY - lastMouseY) / dt;
                smoothFactor = 0.2f;

                p_selected_ball->vx = (p_selected_ball->vx * (1.0f - smoothFactor)) + (instantVx * smoothFactor);
                p_selected_ball->vy = (p_selected_ball->vy * (1.0f - smoothFactor)) + (instantVy * smoothFactor);
            }
        }
        

        for (int i = 0; i < ballCount; i++) {
            for (int j = i + 1; j < ballCount; j++) {
                applyGravityBetween(&balls[i], &balls[j], dt);
            }
        }

        // physics
        for(int i = 0; i < ballCount; i++) {
            updatephysics(&balls[i], i, dt);
        }
        
        for (int i = 0; i < ballCount; i++) {
            collideWithWalls(&balls[i], WINDOW_W, WINDOW_H, RESTITUTION, FRICTION_X);
            handle_collisions(i);
        }

        for (int i = 0; i < ballCount; i++){
            trail_add_point(&balls[i]);
        }

        // -------- RENDER --------
        SDL_SetRenderDrawColor(ren, 30, 30, 30, 255);
        SDL_RenderClear(ren);
        
        // colour shenanigans
        hue += 60.0f * dt;  // rotate 60° per second
        if (hue >= 360.0f) hue -= 360.0f;
        Uint8 R, G, B;
        HSVtoRGB(hue, 1.0f, 1.0f, &R, &G, &B);
        SDL_SetRenderDrawColor(ren, R, G, B, 255);

        // render balls & trails
        for (int i = 0; i < ballCount; i++) {
            draw_filled_circle(ren, &balls[i]);
            trail_render(ren, &balls[i], 255, 100, 100);
        }

        SDL_RenderPresent(ren);

        SDL_Delay(1); // slight delay
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();

    return 0;
}

