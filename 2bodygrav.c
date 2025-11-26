#include <SDL3/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define WINDOW_W 1000
#define WINDOW_H 800
#define MAX_TRAIL_POINTS 2000


// Phys constants
const double ACCELDUETOGRAVITY = 0.0; // can give it downward gravity if so desired
const double FOFGRAV = 50.0f;
const double RESTITUTION = 0.8;
const double FRICTION_X = 0.995;
const double DAMPING_COEF = 0.9f;
double hue = 0.0;

typedef struct {
    float x, y;
    float vx, vy;
    float ax, ay;
    float radius;
    float mass;

    Uint8 r, g, b; // colour
} Ball;

typedef struct {
    float x[MAX_TRAIL_POINTS];
    float y[MAX_TRAIL_POINTS];
    int count;
} Trail;


void trail_add_point(Trail *t, float x, float y)
{
    if (t->count >= MAX_TRAIL_POINTS) {
        // shift everything left (like a moving window)
        memmove(&t->x[0], &t->x[1], sizeof(float) * (MAX_TRAIL_POINTS - 1));
        memmove(&t->y[0], &t->y[1], sizeof(float) * (MAX_TRAIL_POINTS - 1));
        t->count = MAX_TRAIL_POINTS - 1;
    }

    t->x[t->count] = x;
    t->y[t->count] = y;
    t->count++;
}


void trail_render(SDL_Renderer *ren, Trail *t, Uint8 r, Uint8 g, Uint8 b)
{
    SDL_SetRenderDrawColor(ren, r, g, b, 255);

    for (int i = 3; i < t->count; i++) {
        SDL_RenderLine(ren,(int)t->x[i - 3], (int)t->y[i - 3], (int)t->x[i], (int)t->y[i]);
    }
}


// SDL_RenderLine(renderer, x1, y1, x2, y2)
void draw_filled_circle(SDL_Renderer *r, int cx, int cy, int radius) 
{
    for (int dy = -radius; dy <= radius; dy++) {
        int y = cy + dy;
        double dx = sqrt((double)(radius * radius) - (double)(dy * dy));
        int x1 = cx - (int)dx;
        int x2 = cx + (int)dx;
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

void handle_collisions(Ball *ball, Ball *centre, float RESTITUTION, float FRICTION_X) {
    // ---------- BALL WALL COLLISIONS ----------
    // Floor
    if (ball->y > WINDOW_H - ball->radius) {
        ball->y = WINDOW_H - ball->radius;
        ball->vy = -ball->vy * RESTITUTION;
        ball->vx *= FRICTION_X;
    }

    // Ceiling
    if (ball->y < ball->radius) {
        ball->y = ball->radius;
        ball->vy = -ball->vy * RESTITUTION;
    }

    // Left wall
    if (ball->x < ball->radius) {
        ball->x = ball->radius;
        ball->vx = -ball->vx * RESTITUTION;
    }

    // Right wall
    if (ball->x > WINDOW_W - ball->radius) {
        ball->x = WINDOW_W - ball->radius;
        ball->vx = -ball->vx * RESTITUTION;
    }


    // ---------- CENTRE WALL COLLISIONS ----------
    if (centre->y > WINDOW_H - centre->radius) {
        centre->y = WINDOW_H - centre->radius;
        centre->vy = -centre->vy * RESTITUTION;
        centre->vx *= FRICTION_X;
    }

    if (centre->y < centre->radius) {
        centre->y = centre->radius;
        centre->vy = -centre->vy * RESTITUTION;
    }

    if (centre->x < centre->radius) {
        centre->x = centre->radius;
        centre->vx = -centre->vx * RESTITUTION;
    }

    if (centre->x > WINDOW_W - centre->radius) {
        centre->x = WINDOW_W - centre->radius;
        centre->vx = -centre->vx * RESTITUTION;
    }


    // ---------- BALL–CENTRE COLLISION ----------
    float dx = ball->x - centre->x;
    float dy = ball->y - centre->y;
    float dist = hypotf(dx, dy);
    float minDist = ball->radius + centre->radius;

    if (dist < minDist) {

        float overlap = minDist - dist;

        // Push the ball out of the centre
        ball->x += (dx / dist) * overlap;
        ball->y += (dy / dist) * overlap;

        // Normal vector
        float nx = dx / dist;
        float ny = dy / dist;

        // Tangent vector
        float tx = -ny;
        float ty = nx;

        // Dot products: tangent
        float dpTan1 = ball->vx * tx + ball->vy * ty;
        float dpTan2 = centre->vx * tx + centre->vy * ty;

        // Dot products: normal
        float dpNorm1 = ball->vx * nx + ball->vy * ny;
        float dpNorm2 = centre->vx * nx + centre->vy * ny;

        // 1D elastic collision along normal
        float m1 = (dpNorm1 * (ball->mass - centre->mass) + 2.0f * centre->mass * dpNorm2) / (ball->mass + centre->mass);

        float m2 = (dpNorm2 * (centre->mass - ball->mass) + 2.0f * ball->mass * dpNorm1) / (ball->mass + centre->mass);

        // Apply updated velocities with slight damping
        ball->vx = (tx * dpTan1 + nx * m1) * DAMPING_COEF;
        ball->vy = (ty * dpTan1 + ny * m1) * DAMPING_COEF;

        centre->vx = (tx * dpTan2 + nx * m2) * DAMPING_COEF;
        centre->vy = (ty * dpTan2 + ny * m2) * DAMPING_COEF;
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

    SDL_SetRenderVSync(ren, 1);

    Ball centre = {
        .x = 500,
        .y = 500,
        .vx = 0,
        .vy = 0,
        .radius = 20,
        .mass = 597000,
        .r = 255, .g = 50, .b = 50
    };

    Ball ball = {
        .x = 100,
        .y = 100,
        .vx = 0,
        .vy = 0,
        .radius = 20,
        .mass = 7347.67,
        .r = 50, .g = 200, .b = 255
    };

    Trail ballTrail = {0};
    Trail centreTrail = {0};

    Uint64 prev = SDL_GetTicks();
    float lastMouseX, lastMouseY;

    bool running = true;
    while (running) {

        // TIMING
        Uint64 now = SDL_GetTicks();
        // milliseconds to seconds
        double dt = (now - prev) / 1000.0f;
        // no div / 0  :(
        if (dt <= 0.0) dt = 0.001f;
        if (dt > 0.1f) dt = 0.1f;
        prev = now;

        // Keypresses n shieet
        SDL_Event e;
        while (SDL_PollEvent(&e)) {

            if (e.type == SDL_EVENT_QUIT)
                running = false;

            if (e.type == SDL_EVENT_KEY_DOWN) {
                SDL_Keycode key = e.key.key;

                if (key == SDLK_ESCAPE) running = false;
                if (key == SDLK_LEFT)  ball.vx -= 50;
                if (key == SDLK_RIGHT) ball.vx += 50;
                if (key == SDLK_UP) ball.vy += -50;
                if (key == SDLK_DOWN) ball.vy += 50;
            } 
        float mouseX, mouseY;
        SDL_MouseButtonFlags buttonState = SDL_GetMouseState(&mouseX, &mouseY);
        float instantVx;
        float instantVy;
        float smoothFactor;

        if (buttonState & SDL_BUTTON_LMASK) { 
                ball.x = mouseX;
                ball.y = mouseY;

                ball.vx = 0.0;
                ball.vy = 0.0;

                lastMouseX = mouseX;
                lastMouseY = mouseY;
                }

        if (e.type == SDL_EVENT_MOUSE_BUTTON_UP){
            instantVx = (mouseX - lastMouseX) / dt;
            instantVy = (mouseY - lastMouseY) / dt;
            smoothFactor = 0.2f;

            ball.vx = (ball.vx * (1.0f - smoothFactor)) + (instantVx * smoothFactor);
            ball.vy = (ball.vy * (1.0f - smoothFactor)) + (instantVy * smoothFactor);
        }
        }
        

        // Apply gravity from centre to ball
        applyGravity(&ball, &centre, dt);
        
        // Apply gravity from ball to centre
        applyGravity(&centre, &ball, dt);


        // Move ball
        ball.x += ball.vx * dt;
        ball.y += ball.vy * dt;
        //move centre
        centre.x += centre.vx * dt;
        centre.y += centre.vy * dt;

        handle_collisions(&ball, &centre, RESTITUTION, FRICTION_X);

        trail_add_point(&ballTrail, ball.x, ball.y);
        trail_add_point(&centreTrail, centre.x, centre.y);

        // -------- RENDER --------
        SDL_SetRenderDrawColor(ren, 30, 30, 30, 255);
        SDL_RenderClear(ren);

        SDL_SetRenderDrawColor(ren, centre.r, centre.g, centre.b, 255);
        draw_filled_circle(ren, (int)centre.x, (int)centre.y, centre.radius);
        trail_render(ren, &centreTrail, 255, 100, 100); // red trail

        hue += 60.0f * dt;  // rotate 60° per second
        if (hue >= 360.0f) hue -= 360.0f;

        Uint8 R, G, B;
        HSVtoRGB(hue, 1.0f, 1.0f, &R, &G, &B);
        SDL_SetRenderDrawColor(ren, R, G, B, 255);
        draw_filled_circle(ren, (int)ball.x, (int)ball.y, ball.radius);
        trail_render(ren, &ballTrail,   0, 255, 255);   // cyan trail

        SDL_RenderPresent(ren);

        SDL_Delay(1); // slight delay
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();

    return 0;
}
