#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

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

  
  bool debug = true;
  int numBins = 10000;
  int numPoints = 10000;
  time_t seed;
  if (debug) {
    // seed = 1456283917;
    seed = time(NULL);
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
      generateUnitVectorFromCartesian(p);
      // printf("Point %d:\n", pointCounter);
      // printf("x: %f\n", p->x);
      // printf("y: %f\n", p->y);
      // printf("z: %f\n", p->z);
      // printf("theta: %f\n", p->theta);
      // printf("phi: %f\n", p->phi);
      // printf("\n");
      pointCounter++;
    }
  }

  float radius = pow((3 * numBins) / (4 * PI), 1.0 / 3.0);
  int circumference = ceil(2.0 * PI * radius);
  printf("radius %f\n", radius);
  printf("circumference %d\n", circumference);

  // Generate Unsorted Maps
  Bin **maps; 
  maps = (Bin **) calloc(circumference, sizeof(Bin *));
  if (maps == NULL) {
    printf("Error allocating memory for the maps\n");
    // free(points);
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
  mergeSort(maps[0], circumference);
  for (mapCounter = 0; mapCounter < circumference; mapCounter++) {
    mergeSort(maps[mapCounter], circumference - mapCounter);
    // printf("\n");
    // for (binCounter = 0; binCounter < circumference - mapCounter; binCounter++) {
      // printf("Index: %d\n", binCounter);
      // printf("offset: %d\n", maps[mapCounter][binCounter].offset);
      // printf("theta: %f\n", maps[mapCounter][binCounter].theta);
    // }
  }

  // Bin Points
  int localCircumference;
  int mapIndex;
  float upperThetaOffset; 
  float lowerThetaOffset;

  int guessIndex;
  UnitVector guess;
  
  p = &points[0];

  guessIndex = floor(numBins * (p->z + 1) / 2);
  generateUnitVectorFromIndex(&guess, guessIndex, numBins);
  localCircumference = ceil(circumference * sqrt(pow(p->x, 2) + pow(p->y, 2)));
  mapIndex = circumference - localCircumference;

  if (guess.theta < p->theta) {
    upperThetaOffset = p->theta - guess.theta;
    lowerThetaOffset = (2 * PI) - upperThetaOffset;
  }
  else {
    lowerThetaOffset = guess.theta - p->theta;
    upperThetaOffset = (2 * PI) - lowerThetaOffset;
  }

  printf("\n");
  printf("Guess Index: %d\n", guessIndex);
  
  printf("\n");
  printf("p.x:     %f\n", p->x);
  printf("p.y:     %f\n", p->y);
  printf("p.z:     %f\n", p->z);
  printf("p.phi:   %f\n", p->phi);
  printf("p.theta: %f\n", p->theta);
  printf("\n");
  printf("guess.x:     %f\n", guess.x);
  printf("guess.y:     %f\n", guess.y);
  printf("guess.z:     %f\n", guess.z);
  printf("guess.phi:   %f\n", guess.phi);
  printf("guess.theta: %f\n", guess.theta);
  printf("\n");
  printf("%f\n", p->theta);
  printf("%f\n", guess.theta);
  printf("\n");
  printf("%f\n", lowerThetaOffset);
  printf("%f\n", upperThetaOffset);

  int lowerAfter = guessIndex - getOffsetAfter(maps[mapIndex], localCircumference, lowerThetaOffset);
  int lowerBefore = guessIndex - getOffsetBefore(maps[mapIndex], localCircumference, lowerThetaOffset);
  int upperBefore = guessIndex + getOffsetBefore(maps[mapIndex], localCircumference, upperThetaOffset);
  int upperAfter = guessIndex + getOffsetAfter(maps[mapIndex], localCircumference, upperThetaOffset);

  printf("\n");
  printf("%d\n", getOffsetAfter(maps[mapIndex], localCircumference, lowerThetaOffset));
  printf("%d\n", getOffsetBefore(maps[mapIndex], localCircumference, lowerThetaOffset));
  printf("%d\n", getOffsetBefore(maps[mapIndex], localCircumference, upperThetaOffset));
  printf("%d\n", getOffsetAfter(maps[mapIndex], localCircumference, upperThetaOffset));
  printf("\n");
  printf("%d\n", lowerAfter);
  printf("%d\n", lowerBefore);
  printf("%d\n", upperBefore);
  printf("%d\n", upperAfter);
  
  printf("\n");
  for (binCounter = 0; binCounter < circumference - mapIndex; binCounter++) {
    printf("Index: %d\n", binCounter);
    printf("offset: %d\n", maps[mapIndex][binCounter].offset);
    printf("theta: %f\n", maps[mapIndex][binCounter].theta);
  }

  int bestPointIndex;
  float bestDistance;
  float testDistance;
  UnitVector testPoint;
  bool found = false;

  if ((lowerAfter <= guessIndex) && (0 <= lowerAfter)) {
    generateUnitVectorFromIndex(&testPoint, lowerAfter, numBins);
    if (!found) {
      bestPointIndex = lowerAfter;
      bestDistance = distanceBetween(&testPoint, p);
    }
    else {
      testDistance = distanceBetween(&testPoint, p);
      if (testDistance < bestDistance) {
        bestDistance = testDistance;
        bestPointIndex = lowerAfter;
      }
    }
    found = true;
  }
  if ((lowerBefore <= guessIndex) && (0 <= lowerBefore)) {
    generateUnitVectorFromIndex(&testPoint, lowerBefore, numBins);
    if (!found) {
      bestPointIndex = lowerBefore;
      bestDistance = distanceBetween(&testPoint, p);
    }
    else {
      testDistance = distanceBetween(&testPoint, p);
      if (testDistance < bestDistance) {
        bestDistance = testDistance;
        bestPointIndex = lowerBefore;
      }
    }
    found = true;
  }
  if ((upperBefore >= guessIndex) && (numBins > upperBefore)) {
    generateUnitVectorFromIndex(&testPoint, upperBefore, numBins);
    if (!found) {
      bestPointIndex = upperBefore;
      bestDistance = distanceBetween(&testPoint, p);
    }
    else {
      testDistance = distanceBetween(&testPoint, p);
      if (testDistance < bestDistance) {
        bestDistance = testDistance;
        bestPointIndex = upperBefore;
      }
    }
    found = true;
  }
  if ((upperAfter >= guessIndex) && (numBins > upperAfter)) {
    generateUnitVectorFromIndex(&testPoint, upperAfter, numBins);
    if (!found) {
      bestPointIndex = upperAfter;
      bestDistance = distanceBetween(&testPoint, p);
    }
    else {
      testDistance = distanceBetween(&testPoint, p);
      if (testDistance < bestDistance) {
        bestDistance = testDistance;
        bestPointIndex = upperAfter;
      }
    }
    found = true;
  }
  printf("Best: %d\n", bestPointIndex);
  printf("Best Distance: %f\n", bestDistance);

  // Clean Up
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
  int offset = -1;
  int binCounter;
  for (binCounter = 0; binCounter < length; binCounter++) {
    if (theta < map[binCounter].theta) {
      return map[binCounter].offset;
    }
  }
  return offset;
}

int getOffsetBefore(Bin *map, int length, float theta) {
  int offset = -1;
  int binCounter;
  for (binCounter = length - 1; binCounter >= 0; binCounter--) {
    if (theta > map[binCounter].theta) {
      return map[binCounter].offset;
    }
  }
  return offset;
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
  printf("\n");
  printf("x1: %f\n", p1->x);
  printf("x2: %f\n", p2->x);
  printf("xDiff: %f\n", p1->x - p2->x);
  printf("y1: %f\n", p1->y);
  printf("y2: %f\n", p2->y);
  printf("yDiff: %f\n", p1->y - p2->y);
  printf("z1: %f\n", p1->z);
  printf("z2: %f\n", p2->z);
  printf("zDiff: %f\n", p1->z - p2->z);
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