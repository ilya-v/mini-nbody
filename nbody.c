#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"

#define SOFTENING 1e-9f

typedef struct { float x, y, z, vx, vy, vz; } Particle;

void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

void bodyForce(Particle *p, float dt, int n) {
  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }

    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
}

int main(const int argc, const char** argv) {

  int nBodies = 30000;
  if (argc > 1) nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 100;  // simulation iterations

  int bytes = nBodies*sizeof(Particle);
  float *buf = (float*)malloc(bytes);
  Particle *p = (Particle*)buf;

  randomizeBodies(buf, 6*nBodies); // Init pos / vel data

  double totalTime = 0.0;
  StartTimer();

  for (int iter = 1; iter <= nIters; iter++) {

    bodyForce(p, dt, nBodies); // compute interbody forces

    for (int i = 0 ; i < nBodies; i++) { // integrate position
      p[i].x += p[i].vx*dt;
      p[i].y += p[i].vy*dt;
      p[i].z += p[i].vz*dt;
    }
  }
  totalTime = GetTimer() / 1000.0;
  double avgTime = totalTime / (double)(nIters-1); 

  printf("N=%d, Titer=%0.3f s\n", nBodies, avgTime);
  free(buf);
}
