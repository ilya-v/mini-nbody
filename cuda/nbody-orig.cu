#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"

#define BLOCK_SIZE 256
#define SOFTENING 1e-9f

typedef struct { float x, y, z, vx, vy, vz; } Particle;

__global__
void calcForces(Particle *p, float dt, int N) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

        for (int j = 0; j < N; j++) {
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
            float invDist = rsqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        p[i].vx += dt*Fx;
        p[i].vy += dt*Fy;
        p[i].vz += dt*Fz;
    }
}

int main(const int argc, const char** argv) {

    const int
        N = 30000,
        nSteps = 100,
        nBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;

    const float dt = 0.01f; // time step

    size_t n_bytes = N*sizeof(Particle);
    Particle *particles = (Particle*)malloc(n_bytes);

    for (int i = 0; i < N; i++) {
        Particle *p = particles + i;
        p->x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p->y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p->z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p->vx = 0;
        p->vy = 0;
        p->vz = 0;
    }

    Particle *d_p;
    cudaMalloc(&d_p, n_bytes);

    StartTimer();
    for (int iter = 1; iter <= nSteps; iter++) {

        cudaMemcpy(d_p, particles, n_bytes, cudaMemcpyHostToDevice);
        calcForces <<<nBlocks, BLOCK_SIZE>>>(d_p, dt, N);
        cudaMemcpy(particles, d_p, n_bytes, cudaMemcpyDeviceToHost);

        for (int i = 0 ; i < N; i++) { // integrate position
            particles[i].x += particles[i].vx*dt;
            particles[i].y += particles[i].vy*dt;
            particles[i].z += particles[i].vz*dt;
        }
    }

    printf("N=%d, Titer=%0.3f s\n",
        N,
        GetTimer() / nSteps / 1000.0);

    free(particles);
    cudaFree(d_p);
}
