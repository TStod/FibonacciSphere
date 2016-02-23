#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

float randomFloat();
float randomFloatMax(float max);
float randomFloatMinMax(float min, float max);

struct UnitVector {
  float x;
  float y;
  float z;
  float theta; // latitude in radians
  float phi; // longitude in radians
  float latitude;
  float longitude;
};


// arg[1]: number of bins to generate
// arg[2]: number of points to generate
// arg[3]: seed (optional)
int main(int argc, char **argv) {
  bool debug = true;
  int numBins = 10000;
  int numPoints = 10000000;
  time_t seed;
  if (debug) {
    seed = 123456789;
    seed = time(NULL);
    // seed = 1456190831;
  }
  else {
    if ((argc > 4) || (3 > argc)) {
      printf("Need 2 or 3 arguments:\nBins\nPoints\nSeed (optional)\n");
      return 1;
    }
    numBins = atoi(argv[1]);
    numPoints = atoi(argv[2]);
    if (3 < argc) {
      seed = atoi(argv[3]);
    }
    else {
      seed = time(NULL);
    }
  }


  printf("Bins: %d\n", numBins);
  printf("Points: %d\n", numPoints);
  printf("Seed: %d\n", seed);
  srand(seed);

  UnitVector *points;
  points = (UnitVector *) malloc(numPoints * sizeof(UnitVector));
  if (points == NULL) {
    printf("Error allocating memory for the array of points\n");
    return 1;
  }
  else {
    printf("Memory allocated for the array of points\n");
  }

  float x;
  float y;
  float z;
  float mag;
  
  UnitVector *p = NULL;
  int pointCounter = 0;
  while (pointCounter < numPoints) {
    x = randomFloatMinMax(-1, 1);
    y = randomFloatMinMax(-1, 1);
    z = randomFloatMinMax(-1, 1);
    mag = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    if ((0 < mag) && (mag < 1)) { // skip if all random numbers are 0 or if the random point is outside of the unit sphere
      p = &points[pointCounter];
      p->x = x / mag;
      p->y = y / mag;
      p->z = z / mag;
      p->theta = acos(p->z);
      p->phi = atan2(p->y, p->x);
      p->latitude = (p->theta * 180.0 / PI) - 90.0;
      p->longitude = p->phi * 180.0 / PI;
      // printf("\n");
      // printf("Point %d:\n", pointCounter);
      // printf("x: %f\n", p->x);
      // printf("y: %f\n", p->y);
      // printf("z: %f\n", p->z);
      // printf("theta: %f\n", p->theta);
      // printf("phi: %f\n", p->phi);
      // printf("latitude: %f\n", p->latitude);
      // printf("longitude: %f\n", p->longitude);
      // printf("\n");
      pointCounter++;
    }
  }



  free(points);
  return 0;
}

float randomFloat() {
  return rand() / (float) RAND_MAX;
};

float randomFloatMax(float max) {
  return max * randomFloat();
};

float randomFloatMinMax(float min, float max) {
  return randomFloatMax(max - min) + min;
};