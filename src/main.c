#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <unistd.h>
#include "header/linAlg.h"

#define STATES 6
#define MEASURED 2
#define DT 0.1
#define SCALAR 2

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

void drawEstimatedPosition(SDL_Renderer* renderer, mousePosition *head)
{
	SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
	mousePosition *prev = head;
	mousePosition *current = head->next;
	while (current != NULL)
	{
		SDL_RenderDrawLine(renderer, prev->x, prev->y, current->x, current->y);
		prev = current;
		current = current->next;
	}
}
void drawMousePos(SDL_Renderer* renderer, mousePosition *head)
{
	SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
	mousePosition *prev = head;
	mousePosition *current = head->next;
	while (current != NULL)
	{
		SDL_RenderDrawLine(renderer, prev->x, prev->y, current->x, current->y);
		prev = current;
		current = current->next;
	}
}

		/* KALMAN */

double generateRandomNumber(int min, int max) 
{
	double fraction = (double)rand() / RAND_MAX;
	double randNumber = min + (fraction * (min - max));
	return randNumber;
}

double generateRandomNormal()
{
	double u1 = rand() / (double)RAND_MAX;
	double u2 = rand() / (double)RAND_MAX;

	double z0 = sqrt( -2.0 * log(u1))* cos(2.0 * M_PI * u2);
	return z0;
}

void kalmanInitSetZero(kalmanData *kf)
{
	int mouseX, mouseY;
	SDL_PumpEvents();
	SDL_GetMouseState(&mouseX, &mouseY);
	kf->X = allocateVector(STATES);
	kf->fx = allocateVector(STATES);
	kf->hx = allocateVector(MEASURED);
	setZeroVector(kf->X );
	setZeroVector(kf->fx);
	setZeroVector(kf->hx);

	kf->P = allocateMatrix(STATES, STATES);
	kf->Q = allocateMatrix(STATES, STATES);
	kf->R = allocateMatrix(MEASURED, MEASURED);
	kf->K = allocateMatrix(STATES, MEASURED);
	kf->F = allocateMatrix(STATES, STATES);
	kf->H = allocateMatrix(MEASURED, STATES);
	kf->Ht = allocateMatrix(STATES, MEASURED);
	kf->Ft = allocateMatrix(STATES, STATES);
	kf->Pp = allocateMatrix(STATES, STATES);
	setZero(kf->P );
	setZero(kf->Q );
	setZero(kf->K );
	setZero(kf->F );
	setZero(kf->H );
	setZero(kf->Ht);
	setZero(kf->Ft);
	setZero(kf->Pp);

	kf->X->data[0] = mouseX;
	kf->X->data[1] = 0;
	kf->X->data[2] = 0.1;
	kf->X->data[3] = mouseY;
	kf->X->data[4] = 0;
	kf->X->data[5] = 0.1;

	//randomly chosen
	for (int i = 0; i < kf->F->rows; i++)
	{
		for (int j = 0; j < kf->F->cols; j++)
		{
			if (i == j)
			{
				kf->F->data[ i * (kf->F->cols) + j] = 1;
			}
			else
			{
				kf->F->data[ i * (kf->F->cols) + j] = 0;
			}
		}
	}
	kf->F->data[0*kf->F->cols + 1] = DT;
	kf->F->data[0*kf->F->cols + 2] = 0.5*DT*DT;
	kf->F->data[1*kf->F->cols + 2] = DT;

	kf->F->data[3*kf->F->cols + 4] = DT;
	kf->F->data[3*kf->F->cols + 5] = 0.5*DT*DT;
	kf->F->data[4*kf->F->cols + 5] = DT;

	setIdentity(kf->R);
	setIdentity(kf->P);
	setIdentity(kf->Q);
	setIdentity(kf->R);

	for (int i = 0; i < kf->P->rows; i++)
	{
		for (int j = 0; j < kf->P->cols; j++)
		{
			kf->P->data[i * kf->P->cols + i] = 500;
		}
	}

	for (int i = 0; i < kf->R->rows; i++)
	{
		for (int j = 0; j < kf->R->cols; j++)
		{
			kf->R->data[ i * (kf->R->cols) + i] = 0.2;
		}
	}
	

	double Hmat[MEASURED * STATES] =
	{
		1, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0
	};

	for (int i = 0; i < kf->H->rows; i++)
	{
		for (int j = 0; j < kf->H->cols; j++)
		{
			kf->H->data[i * kf->H->cols + j] = Hmat[i * kf->H->cols + j];
		}
	}

	double Qmat[STATES * STATES] = 
	{
		0.25*DT*DT*DT*DT, 0.5*DT*DT*DT, 0.5*DT*DT, 0 							 , 0 					 , 0        ,
		0.5*DT*DT*DT    , DT*DT       , DT       , 0 							 , 0 					 , 0        ,
		0.5*DT*DT       , DT          , 1 			 , 0 							 , 0 					 , 0        ,
		0     					, 0 					, 0 			 , 0.25*DT*DT*DT*DT, 0.5*DT*DT*DT, 0.5*DT*DT,
		0     					, 0 					, 0 			 , 0.5*DT*DT*DT    , DT*DT       , DT       ,
		0     					, 0 					, 0 			 , 0.5*DT*DT       , DT          , 1
	};
	
	for (int i = 0; i < (kf->Q->rows); i++)
	{
		for (int j = 0; j < kf->Q->cols; j++)
		{
			kf->Q->data[ i * kf->Q->cols + j] = Qmat[i * STATES + j];
		}
	}
}
void freeKalman(kalmanData *kf)
{
	freeMatrix(kf->P);
	freeMatrix(kf->Q); //Causes an Error?
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
	//printf("Temp1 of Covariance Prediction");
	//printMatrix(temp1);
	matrixTranspose(kf->F, kf->Ft);
	matrixMultiply(temp1, kf->Ft, kf->Pp);
	matrixAdd(kf->Pp, kf->Q, kf->Pp);

	Matrix *temp2 = allocateMatrix(STATES, MEASURED);
	Matrix *temp3 = allocateMatrix(MEASURED, STATES);
	Matrix *temp4 = allocateMatrix(MEASURED, MEASURED);
	matrixTranspose(kf->H, kf->Ht);
	matrixMultiply(kf->Pp, kf->Ht, temp2);
	matrixMultiply(kf->H, kf->Pp, temp3);
	matrixMultiply(temp3, kf->Ht, temp4);
	matrixAdd(temp4, kf->R, temp4);
	//printf("Temp4 of Kalman update Pre Inversion\n");
	//printMatrix(temp4);
	matrixInverse(temp4, temp4);
	//printf("Temp4 of Kalman update Post Inversion\n");
	//printMatrix(temp4);
	matrixMultiply(temp2, temp4, kf->K);

	//printf("Temp2 of Kalman update\n");
	//printMatrix(temp2);
	//printf("Temp3 of Kalman update\n");
	//printMatrix(temp3);
	//printf("Temp4 of Kalman update\n");
	//printMatrix(temp4);


	//State Update
	Vector *temp5 = allocateVector(MEASURED);
	Vector *temp6 = allocateVector(STATES);
	Vector *temp7 = allocateVector(MEASURED);
	vectorMatrixMultiply(kf->H, kf->X, temp7);
	vectorSubtract(temp7, data, temp5);
	vectorMatrixMultiply(kf->K, temp5, temp6); 
	vectorAdd(kf->fx, temp6, kf->X);

	printf("Temp5 of state update\n");
	printVector(temp5);
	printf("Temp6 of state update\n");
	printVector(temp6);
	printf("Temp7 of state update\n");
	printVector(temp7);



	
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


void model(kalmanData *kf)
{
	vectorMatrixMultiply(kf->F, kf->X, kf->fx);
}

void getMouseInput(mousePosition** head)
{
	int mouseX, mouseY;
	SDL_PumpEvents();
	SDL_GetMouseState(&mouseX, &mouseY);

	double noiseX = generateRandomNormal()*SCALAR;
	double noiseY = generateRandomNormal()*SCALAR;
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

void printStateOfKalman(kalmanData *kf)
{
	printf("user defined measurement Function Hx\n");                                                                   
	printVector(kf->hx);
	printf("Predicition error covariance P\n");                                                                      	 
	printMatrix(kf->P);                                                                                               	 
	printf("Process noise covariance Q\n");                                                                          	 
	printMatrix(kf->Q);                                                                                               	 
	printf("Measuremnent error covariance R\n");                                                                     	 
	printMatrix(kf->R);                                                                                               	 
	printf("Kalman Gain K\n");                                                                                       	 
	printMatrix(kf->K);                                                                                               	 
	printf("Jacobian of process model F\n");                                                                         	 
	printMatrix(kf->F);                                                                                               	 
	printf("Jacobian of Measurement Model H\n");                                                                     	 
	printMatrix(kf->H);                                                                                               	 
	printf("Transposes Ht\n");                                                                                       	 
	printMatrix(kf->Ht);                                                                                              	 
	printf("Transposes Ft\n");                                                                                       	 
	printMatrix(kf->Ft);                                                                                              	 
	printf("Post prediction, preupdate Covariance Pp\n");                                                            	 
	printMatrix(kf->Pp);                                                                                              	 
	printf("State Vector X\n");
	printVector(kf->X);                                                                                                  
	printf("User Defined state transition Xp\n");                                                                    	 
	printVector(kf->fx);                                                                                                 
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
 mousePosition *head = NULL;
 mousePosition *head2 = NULL;
 Vector *mouse = allocateVector(MEASURED);
 

 SDL_Event event;
 int quit = 0;
 int start = 0;
 while (!quit)
 {
	 while(SDL_PollEvent(&event))
	 {
		 if (event.type == SDL_QUIT)
		 {
			 quit = 1;
		 }
		 else if (event.type == SDL_KEYDOWN)
		 {
			 if (event.key.keysym.sym == SDLK_q)
			 {
				 quit = 1;
			 }
			 else if (event.key.keysym.sym == SDLK_y)
			 {
				 start = 1;
			 }
		 }
	 }

	 if (start)
	 {
		getMouseInput(&head);
		addNode(&head2, kf.X->data[0], kf.X->data[3]);
		int mouseX, mouseY;                                                                                                  	
		SDL_GetMouseState(&mouseX, &mouseY);                                                                                 	
		mouse->data[0] = (double)mouseX + generateRandomNormal()*SCALAR;
		mouse->data[1] = (double)mouseY + generateRandomNormal()*SCALAR;
		printf("Mouse Position\n");                                                                                      	  	
		printVector(mouse);                                                                                              	  	
		                                                                                                                  	 	
		kalman(&kf, mouse);                                                                                              	  	
		model(&kf);                                                                                                      	  	
		printStateOfKalman(&kf);                                                                                             	
		                                                                                                                 	  	
		SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);                                                            	  	
		//SDL_RenderDrawPoint(renderer, (int)kf.X->data[0], (int)kf.X->data[3]);                                         	  	
		//SDL_RenderDrawLine(renderer, (int)kf.X->data[0], (int)kf.X->data[3], (int)kf.fx->data[0], (int)kf.fx->data[1]);    	
		drawMousePos(renderer, head);                                                                                    	  	
		drawEstimatedPosition(renderer, head2);                                                                          	  	
		                                                                                                                  	 	
		SDL_RenderPresent(renderer);                                                                                         	
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);                                                                      	
		SDL_RenderClear(renderer);                                                                                           	
		// Delay for stability (optional)                                                                                   
		SDL_Delay(10);                                                                                                      
	 }
	 else
	 {
		 printf("PRESS Y TO START\n");
		 sleep(1);
	 }
 }
 free(head);
 free(head2);
 free(mouse);
 freeKalman(&kf);
 SDL_DestroyRenderer(renderer);
 SDL_DestroyWindow(window);
 SDL_Quit();
	return 0;
}

