#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"

#define BLOCK_SIZE 256
#define SOFTENING 1e-9f

typedef struct { float x, y, z, vx, vy, vz; } Body;

void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

__global__
void bodyForce(Body *p, float dt, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < n) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = rsqrtf(distSqr);
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

  int bytes = nBodies*sizeof(Body);
  float *buf = (float*)malloc(bytes);
  Body *p = (Body*)buf;

  randomizeBodies(buf, 6*nBodies); // Init pos / vel data

  float *d_buf;
  cudaMalloc(&d_buf, bytes);
  Body *d_p = (Body*)d_buf;

  int nBlocks = (nBodies + BLOCK_SIZE - 1) / BLOCK_SIZE;

  StartTimer();
  for (int iter = 1; iter <= nIters; iter++) {

    cudaMemcpy(d_buf, buf, bytes, cudaMemcpyHostToDevice);
    bodyForce<<<nBlocks, BLOCK_SIZE>>>(d_p, dt, nBodies); // compute interbody forces
    cudaMemcpy(buf, d_buf, bytes, cudaMemcpyDeviceToHost);

    for (int i = 0 ; i < nBodies; i++) { // integrate position
      p[i].x += p[i].vx*dt;
      p[i].y += p[i].vy*dt;
      p[i].z += p[i].vz*dt;
    }

  }
  double avgTime = GetTimer() / (double)(nIters-1) / 1000.0; 

  printf("N=%d, Titer=%0.3f s\n", nBodies, avgTime);
  free(buf);
  cudaFree(d_buf);
}
