#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <float.h>
#include <openacc.h>

#define PI 3.14159265358979323846
const double PHI = ((sqrt(5) - 1) / 2); // not to be confused with phi
const double GA = PHI * 2.0 * PI;

struct UnitVector {
  float x;
  float y;
  float z;
  float theta; // longitude in radians [-PI, PI]
  float phi; // latitude in radians [0, PI]
};

struct Bin {
  int offset;
  float theta;
};


// Declare Functions
float randomFloat();
float randomFloatMax(float max);
float randomFloatMinMax(float min, float max);
int getOffsetAfter(Bin *map, int length, float phi);
int getOffsetBefore(Bin *map, int length, float phi);
void merge(Bin *merged, Bin *left, int leftCount, Bin *right, int rightCount);
void mergeSort(Bin *array, int length);
void generateUnitVectorFromIndex(UnitVector *p, int index, int length);
void generateUnitVectorFromCartesian(UnitVector *p);
void generateUnitVectorFromSpherical(UnitVector *p);
float distanceBetween(UnitVector *p1, UnitVector *p2);


// arg[1]: number of bins to generate
// arg[2]: number of points to generate
// arg[3]: seed (optional)
int main(int argc, char **argv) {
  bool debug = false;
  int numBins = 100000;
  int numPoints = 1000000;
  time_t seed = time(NULL);
  if (debug) {
    seed = 123456789;
  }
  else {
    if ((argc > 5) || (3 > argc)) {
      printf("Need 2 or 3 arguments:\nBins\nPoints\nSeed (optional)\n");
      return 1;
    }
    numBins = atoi(argv[1]);
    numPoints = atoi(argv[2]);
    if (3 < argc) {
      seed = atoi(argv[3]);
    }
  }


  printf("Bins: %d\n", numBins);
  printf("Points: %d\n", numPoints);
  printf("Seed: %d\n", seed);
  srand(seed);

  struct timeval startGenPoints;
  struct timeval endGenPoints;
  gettimeofday(&startGenPoints, NULL);

  UnitVector *points;
  points = (UnitVector *) malloc(numPoints * sizeof(UnitVector));
  if (points == NULL) {
    printf("Error allocating memory for the array of points\n");
    return 1;
  }

  float x;
  float y;
  float z;
  float mag;

  // Generate Points
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
      generateUnitVectorFromCartesian(p);
      pointCounter++;
    }
  }

  gettimeofday(&endGenPoints, NULL);
  struct timeval startGenMaps;
  struct timeval endGenMaps;
  gettimeofday(&startGenMaps, NULL);

  float radius = pow((3 * numBins) / (4 * PI), 1.0 / 3.0);
  float exactCircumference = 2.0 * PI * radius;
  int circumference = ceil(exactCircumference);

  // Generate Unsorted Maps
  Bin **maps; 
  maps = (Bin **) calloc(circumference, sizeof(Bin *));
  if (maps == NULL) {
    printf("Error allocating memory for the maps\n");
    return 1;
  }
  maps[0] = (Bin *) malloc(circumference * sizeof(Bin));
  if (maps[0] == NULL) {
    printf("Error allocating memory for map with circumference %d\n", circumference);
    return 1;
  }
  Bin *map = maps[0];
  Bin *bin = NULL;
  int binCounter;
  for (binCounter = 0; binCounter < circumference; binCounter++) {
    bin = &map[binCounter];
    bin->theta = fmod(GA * binCounter, PI * 2.0); // [0, 2PI] only positive because we either add or subtract offsets
    bin->offset = binCounter;
  }
  int mapCounter;
  for (mapCounter = 1; mapCounter < circumference; mapCounter++) {
    maps[mapCounter] = (Bin *) malloc((circumference - mapCounter) * sizeof(Bin));
    if (maps[mapCounter] == NULL) {
      printf("Error allocating memory for map with circumference %d\n", circumference - mapCounter);
      return 1;
    }
    memcpy(maps[mapCounter], maps[mapCounter - 1], (circumference - mapCounter) * sizeof(Bin));
  }
  
  // Sort Maps
  #pragma acc parallel
  {
    #pragma acc loop
    for (mapCounter = 0; mapCounter < circumference; mapCounter++) {
      mergeSort(maps[mapCounter], circumference - mapCounter);
    }
  }
  
  gettimeofday(&endGenMaps, NULL);
  struct timeval startBinPoints;
  struct timeval endBinPoints;
  gettimeofday(&startBinPoints, NULL);

  // Bin Points
  int *results = (int *) calloc(numBins, sizeof(int));
  if (results == NULL) {
    printf("Error allocating memory for results\n");
    return 1;
  }
  int numPotentials = 4;
  int badPoints = 0;

  #pragma acc parallel
  {
    int localCircumference;
    int mapIndex;
    float upperThetaOffset; 
    float lowerThetaOffset;
    int guessIndex;
    UnitVector guess;
    int bestPointIndex;
    float bestDistance;
    float testDistance;
    UnitVector testPoint;
    bool found;
    int guessCounter;
    #pragma acc loop reduction (+ : badPoints)
    for (pointCounter = 0; pointCounter < numPoints; pointCounter++) {
      int *potentials = (int*) malloc(numPotentials * sizeof(int));
      // if (potentials == NULL) {
      //   printf("Error allocating memory for potentials\n");
      // }
      p = &points[pointCounter];
      guessIndex = floor(numBins * (p->z + 1) / 2);
      generateUnitVectorFromIndex(&guess, guessIndex, numBins);
      localCircumference = ceil(exactCircumference * sqrt(pow(p->x, 2) + pow(p->y, 2)));
      if (localCircumference > circumference) {
        localCircumference = circumference;
      }
      mapIndex = circumference - localCircumference;

      if (guess.theta < p->theta) {
        upperThetaOffset = p->theta - guess.theta;
        lowerThetaOffset = (2 * PI) - upperThetaOffset;
      }
      else {
        lowerThetaOffset = guess.theta - p->theta;
        upperThetaOffset = (2 * PI) - lowerThetaOffset;
      }

      potentials[0] = getOffsetAfter(maps[mapIndex], localCircumference, lowerThetaOffset);
      potentials[1] = getOffsetBefore(maps[mapIndex], localCircumference, lowerThetaOffset);
      potentials[2] = getOffsetBefore(maps[mapIndex], localCircumference, upperThetaOffset);
      potentials[3] = getOffsetAfter(maps[mapIndex], localCircumference, upperThetaOffset);
      
      found = false;
      bestDistance = FLT_MAX;
      for (guessCounter = 0; guessCounter < numPotentials; guessCounter++) {
        potentials[guessCounter] *= pow(-1, (guessCounter / 2) + 1);
        potentials[guessCounter] += guessIndex;
        if ((0 <= potentials[guessCounter]) && (potentials[guessCounter] < numBins)) {
          generateUnitVectorFromIndex(&testPoint, potentials[guessCounter], numBins);
          testDistance = distanceBetween(&testPoint, p);
          if (testDistance < bestDistance) {
            bestPointIndex = potentials[guessCounter];
            bestDistance = testDistance;
          }
        }
      }
      if (found) {
        #pragma acc atomic capture
        {
          results[bestPointIndex]++;
        }
      }
      else {
        badPoints++;
      }
      free(potentials);
    }
  }

  gettimeofday(&endBinPoints, NULL);

  if (badPoints > 0) {
    printf("\n");
    printf("Bad Points:%d\n", badPoints);
  }

  printf("\n");
  printf("Timing:\n");
  printf("Point Generation: %f\n", ((double)(1000000 * (endGenPoints.tv_sec - startGenPoints.tv_sec)) + (double)(endGenPoints.tv_usec - startGenPoints.tv_usec)) / (double)1000000);
  printf("Map Generation: %f\n", ((double)(1000000 * (endGenMaps.tv_sec - startGenMaps.tv_sec)) + (double)(endGenMaps.tv_usec - startGenMaps.tv_usec)) / (double)1000000);
  printf("Point Binning: %f\n", ((double)(1000000 * (endBinPoints.tv_sec - startBinPoints.tv_sec)) + (double)(endBinPoints.tv_usec - startBinPoints.tv_usec)) / (double)1000000);

  // printf("\n");
  // printf("Results:\n");
  // for (binCounter = 0; binCounter < numBins; binCounter++) {
  //   printf("%d: %d\n", binCounter, results[binCounter]);
  // }

  // Clean Up
  free(results);
  for (mapCounter = 0; mapCounter < circumference; mapCounter++) {
    free(maps[mapCounter]);
  }
  free(maps);
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

int getOffsetAfter(Bin *map, int length, float theta) {
  int offset = 0;
  int binCounter;
  for (binCounter = 0; binCounter < length; binCounter++) {
    if (theta < map[binCounter].theta) {
      return map[binCounter].offset;
    }
  }
  return map[0].offset;
}

int getOffsetBefore(Bin *map, int length, float theta) {
  int offset = 0;
  int binCounter;
  for (binCounter = length - 1; binCounter >= 0; binCounter--) {
    if (theta > map[binCounter].theta) {
      return map[binCounter].offset;
    }
  }
  return map[length - 1].offset;
}

// will overwrite the unit vector at p
void generateUnitVectorFromIndex(UnitVector *p, int index, int length) {
  p->phi = acos(-1.0 + (2.0 * ((float) (index % length) / (float) (length - 1))));
  p->theta = fmod(GA * index, PI * 2.0); // [0, 2PI]
  if (p->theta > PI) {
    p->theta -= 2 * PI; // [-PI, PI]
  }
  generateUnitVectorFromSpherical(p);
}

// expecting x, y, and z to be filled in already
void generateUnitVectorFromCartesian(UnitVector *p) {
  p->phi = acos(p->z);
  p->theta = atan2(p->y, p->x);
}

// expecting phi and theta to be filled in already
void generateUnitVectorFromSpherical(UnitVector *p) {
  p->x = sin(p->phi) * cos(p->theta);
  p->y = sin(p->phi) * sin(p->theta);
  p->z = cos(p->phi);
}

float distanceBetween(UnitVector *p1, UnitVector *p2) {
  float distance = 0;
  distance += pow(p1->x - p2->x, 2);
  distance += pow(p1->y - p2->y, 2);
  distance += pow(p1->z - p2->z, 2);
  return sqrt(distance);
}

// Merged Sort Adapted From:
// https://gist.github.com/mycodeschool/9678029#file-mergesort_c-c-L25

void merge(Bin *merged, Bin *left, int leftCount, Bin *right, int rightCount) {
  int leftIndex = 0;
  int rightIndex = 0;
  int mergedIndex = 0;

  while ((leftIndex < leftCount) && (rightIndex < rightCount)) {
    if (left[leftIndex].theta > right[rightIndex].theta) { // > ascending, < descending
      memcpy(&merged[mergedIndex++], &right[rightIndex++], sizeof(Bin));
    }
    else {
      memcpy(&merged[mergedIndex++], &left[leftIndex++], sizeof(Bin));
    }
  }
  while (leftIndex < leftCount) {
    memcpy(&merged[mergedIndex++], &left[leftIndex++], sizeof(Bin));
  }
  while (rightIndex < rightCount) {
    memcpy(&merged[mergedIndex++], &right[rightIndex++], sizeof(Bin));
  }
}

void mergeSort(Bin *array, int length) {
  int mid;
  int i;
  Bin *left;
  Bin *right;
  if (length < 2) {
    return;
  }

  mid = length / 2;

  left = (Bin*) malloc(mid * sizeof(Bin));
  right = (Bin*) malloc((length - mid) * sizeof(Bin)); 
  // TODO: check to make sure allocation worked?
  
  memcpy(left, array, mid * sizeof(Bin));
  memcpy(right, &array[mid], (length - mid) * sizeof(Bin));
  
  mergeSort(left, mid);
  mergeSort(right, length - mid);
  merge(array, left, mid, right, length - mid);

  free(left);
  free(right);
}