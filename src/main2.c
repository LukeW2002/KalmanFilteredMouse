#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include "header/linAlg.h"

#define STATES 4
#define MEASURED 2
#define DT 0.1

typedef struct mousePosition {
	int x;
	int y;
	struct mousePosition *next;
} mousePosition;

		/* MATHS STRUCTURES DEFINITIONS */
typedef struct {
	Vector *X;    // State Vector
	Vector *Xp;   // Predicted State Vector
	Matrix *P; 		// Prediction Error Covariance */
	Matrix *Q; 		// Process noise Covariance
	Matrix *R; 		// Measurement Error Covariance 
	Matrix *K; 		// Kalman Gain
	Matrix *F; 		// Jacobian of process Noise
	Matrix *H; 		// Jacobian of Measurenment Noise
	Matrix *Ht; 	// Transpose of the Jacobian
	Matrix *Ft; 	// Transpose of the process Jacobian
	Matrix *Pp; 	// Prediction Errror Covariance after prediction but before the update
	Vector *fx; 	// Output
	Vector *hx; 	// Measurement Function
}kalmanData;

mousePosition *createNode(int x, int y) 
{
	mousePosition *newNode = (mousePosition*)malloc(sizeof(mousePosition));
	if (newNode == NULL)
	{
		exit(1);
	}
	newNode->x = x;
	newNode->y = y;
	newNode->next = NULL;
	return newNode;
}

void addNode(mousePosition** head, int x, int y) 
{
	mousePosition *newNode = createNode(x, y);
	if (*head == NULL)
	{
		*head = newNode;
	}
	else
	{
		mousePosition *temp = *head;
		while (temp->next != NULL)
		{
			temp = temp->next;
		}
		temp->next = newNode;
	}
}


void freeList(mousePosition *head)
{
	mousePosition *temp;
	while (head != NULL)
	{
		temp = head;
		head = head->next;
		free(temp);
	}
}

void drawMousePos(SDL_Renderer* renderer, mousePosition *head)
{
	SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
	mousePosition *temp = head;
	while (temp != NULL)
	{
		SDL_RenderDrawPoint(renderer, temp->x, temp->y);
		temp = temp->next;
	}
}

		/* KALMAN */

double generateRandomNumber(int min, int max) 
{
	double fraction = (double)rand() / RAND_MAX;
	double randNumber = min + (fraction * (min - max));
	

	return randNumber;
}

void kalmanInitSetZero(kalmanData *kf)
{
	int mouseX, mouseY;
	SDL_PumpEvents();
	SDL_GetMouseState(&mouseX, &mouseY);
	kf->X = allocateVector(STATES);
	kf->P = allocateMatrix(STATES, STATES);
	kf->Q = allocateMatrix(STATES, STATES);
	kf->R = allocateMatrix(MEASURED, MEASURED);
	kf->K = allocateMatrix(STATES, MEASURED);
	kf->F = allocateMatrix(STATES, STATES);
	kf->H = allocateMatrix(MEASURED, STATES);
	kf->Ht = allocateMatrix(STATES, MEASURED);
	kf->Ft = allocateMatrix(STATES, STATES);
	kf->Pp = allocateMatrix(STATES, STATES);
	kf->fx = allocateVector(STATES);
	kf->hx = allocateVector(MEASURED);

	kf->X->data[0] = mouseX;
	kf->X->data[1] = 0.2;
	kf->X->data[2] = mouseY;
	kf->X->data[3] = 0.2;

	//randomly chosen
	kf->F->data[0] = 1; 
	kf->F->data[1] = 1; 
	kf->F->data[2] = 0; 
	kf->F->data[3] = 1; 
	setIdentity(kf->Q);
	for (int i = 0; i < kf->Q->rows; i++)
	{
		for (int j = 0; j < kf->Q->cols; j++)
		{
			kf->Q->data[ i * (kf->Q->cols) + j] = kf->Q->data[ i * (kf->Q->cols) + j] * 0.2;
		}
	}

	for (int i = 0; i < kf->R->rows; i++)
	{
		for (int j = 0; j < kf->R->cols; j++)
		{
			kf->R->data[ i * (kf->R->cols) + j] = kf->R->data[ i * (kf->R->cols) + j] * 0.2;
		}
	}

	for (int i = 0; i < kf->P->rows; i++)
	{
		for (int j = 0; j < kf->P->cols; j++)
		{
			kf->P->data[ i * (kf->P->cols) + j] = kf->P->data[ i * (kf->P->cols) + j] * 0.2;
		}
	}
	setZero(kf->P);
	setZero(kf->Q);
	setZero(kf->R);
	setZero(kf->K);
	setZero(kf->F);
	setZero(kf->H);
}
void freeKalman(kalmanData *kf)
{
	freeMatrix(kf->P);
	freeMatrix(kf->Q);
	freeMatrix(kf->R);
	freeMatrix(kf->K);
	freeMatrix(kf->F);
	freeMatrix(kf->H);
	freeMatrix(kf->Ht);
	freeMatrix(kf->Ft);
	freeMatrix(kf->Pp);
	freeVector(kf->X);
	freeVector(kf->fx);
	freeVector(kf->hx);
}
void kalman(kalmanData *kf, Vector *data)
{
	// Covariance Prediciton
	// P_{Kprediction} = F_{k-1}\times P_{k-1}\times F^{T}_{k-1} + Q_{T-1}
	Matrix *temp1 = allocateMatrix(STATES, STATES);

	matrixMultiply(kf->F, kf->P, temp1);

	matrixTranspose(kf->F, kf->Ft);
	matrixMultiply(temp1, kf->Ft, kf->Pp);

	matrixAdd(kf->Pp, kf->Q, kf->Pp);

	Matrix *temp2 = allocateMatrix(STATES, MEASURED);
	Matrix *temp3 = allocateMatrix(MEASURED, STATES);
	Matrix *temp4 = allocateMatrix(MEASURED, MEASURED);
	matrixTranspose(kf->H, kf->Ht);
	//printf("kalman Data: %i, %i\n", kf->H->cols, kf->H->rows);
	//printf("kalman Data: %i, %i\n", kf->Ht->cols, kf->Ht->rows);
	matrixMultiply(kf->Pp, kf->Ht, temp2);
	matrixMultiply(kf->H, kf->Pp, temp3);
	matrixMultiply(temp3, kf->Ht, temp4);
	matrixAdd(temp4, kf->R, temp4);
	matrixInverse(temp4, temp4);
	matrixMultiply(temp2, temp4, kf->K);

	Vector *temp5 = allocateVector(MEASURED);
	Vector *temp6 = allocateVector(STATES);
	vectorSubtract(data, kf->hx, temp5);
	vectorMatrixMultiply(kf->K, temp5, temp6);
	vectorAdd(kf->fx, temp6, kf->X);
	
	Matrix *id = allocateMatrix(STATES, STATES);
	setIdentity(id);
	matrixMultiply(kf->K, kf->H, temp1);
	matrixSubtract(id, temp1, temp1);
	matrixMultiply(temp1, kf->Pp, kf->P);
	freeMatrix(temp1);
	freeMatrix(temp2);
	freeMatrix(temp3);
	freeMatrix(temp4);
	freeVector(temp5);
	freeVector(temp6);
}


void model(kalmanData *kf, Vector *data)
{
	kf->fx->data[0] = kf->X->data[0] + DT * kf->X->data[1];
	kf->fx->data[1] = kf->X->data[1];
	kf->fx->data[2] = kf->X->data[2] + DT * kf->X->data[3];
	kf->fx->data[3] = kf->X->data[3];

	kf->F->data[1] = DT;
	kf->F->data[5] = DT;

	kf->hx->data[0] = kf->fx->data[0];
	kf->hx->data[1] = kf->fx->data[2];

	double dx = kf->hx->data[0] - data->data[0];
	double dy = kf->hx->data[1] - data->data[1];

	kf->H->data[0] = dx; 
	kf->H->data[3] = dy; 
}

void getMouseInput(mousePosition** head)
{
	int mouseX, mouseY;
	SDL_PumpEvents();
	SDL_GetMouseState(&mouseX, &mouseY);

	double noiseX = generateRandomNumber(-5, 5);
	double noiseY = generateRandomNumber(-5, 5);
	mouseX += (int) noiseX;
	mouseY += (int) noiseY;

	addNode(head, mouseX, mouseY);
}

void getCurrentPos(mousePosition *head, int *x, int *y)
{
	if (head == NULL)
	{
		*x = -1;
		*y = -1;
	}
	else
	{
		mousePosition *temp = head;
		while (temp->next != NULL)
		{
			temp = temp->next;
		}

		*x = temp->x;
		*y = temp->y;
	}
}

int main(void)
{
    // Initialize SDL
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
      fprintf(stderr, "SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
      return 1;
  }

  // Create window (not needed for mouse input)
  SDL_Window *window = SDL_CreateWindow("Kalman Filter with SDL2 Mouse Input",
                                        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                                        640, 480, SDL_WINDOW_SHOWN);
  if (window == NULL) {
     fprintf(stderr, "Window could not be created! SDL_Error: %s\n", SDL_GetError());
     return 1;
  }

 SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
 if (renderer == NULL) {
     fprintf(stderr, "Renderer could not be created! SDL_Error: %s\n", SDL_GetError());
     return 1;
 }
 kalmanData kf;
 kalmanInitSetZero(&kf);
 int k = 0;
 mousePosition *head = NULL;
 mousePosition *kalmanHead = NULL;
 Vector *mouse = allocateVector(MEASURED);

 while( k < 10000)
 {
	 getMouseInput(&head);
	 getMouseInput(&kalmanHead);

	 int mouseX, mouseY;
	 getCurrentPos(head, &mouseX, &mouseY);
	 mouse->data[0] = (double)mouseX;
	 mouse->data[1] = (double)mouseY;

	 kalman(&kf, mouse);
	 model(&kf, mouse);
	 
	 SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	 SDL_RenderDrawPoint(renderer, (int)kf.X->data[0], (int)kf.X->data[1]);
	 drawMousePos(renderer, head);

   SDL_RenderPresent(renderer);
   SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
   SDL_RenderClear(renderer);
	 printf("Please Work Oh god: (%f, %f)\n", kf.X->data[0],kf.X->data[1]);



   // Delay for stability (optional)
   SDL_Delay(10);
	 k++;
 }
 free(mouse);
 freeKalman(&kf);
 SDL_DestroyRenderer(renderer);
 SDL_DestroyWindow(window);
 SDL_Quit();
	return 0;
}

