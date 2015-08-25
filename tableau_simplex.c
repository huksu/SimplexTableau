/*
    Input Format, assumes tableau in standard form,
	meaning a maximization problem with only equalities:
        numrealvars numslackvars numconstraints 
		Constraints
        ...
        Objective Function
		
		The constraints have a basis column.  
		If the entry for the constraint is -1,
		then this constraint has no associated
		basis vector.  Usually this happens from
		a greater-than or equality constraint.
        
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // for frexp

typedef int bool;
#define true 1
#define false 0

typedef int TabErr;
#define OK 0
#define INFEASIBLE -1
#define UNBOUNDED -2

#define EPSILON 1.0e-10 // plenty of space... double epsilon is 2.220446049250313e-16
#define flt_equals(a, b) (fabs((a)-(b)) < EPSILON)
#define flt_lt(a,b) (a+EPSILON < b)

#define PRINTEXP false

typedef struct {
	int rows;
	int cols;
	int numvars;
	int numrealvars;
	int numslackvars;
	int numartificialvars;
	int numconstraints;
	int zcolumn; // todo remove z column... it is pointless?
	int constcolumn;
	int basiscolumn;
	double* data; // this is a 2d array in one dimension for better cache usage
	double* copydata; // this is a 2d array in one dimension for better cache usage
} Tableau;

double frexp10( double num, int * rExponent )
{
	// Get the proper decomposition for base 2.
	int exponent2 = 0;
	double fraction2 = frexp( num, &exponent2 );
	
	// Convert the exponent to an exponent base 10, leaving a remainder.
	double exponent10 = 0.0;
	double remainder10 = modf( exponent2 * M_LN2 / M_LN10,
							   &exponent10 );
	
	// Floor shouldn't do anything here, I think.  If it DOES do anything,
	// then the correction below should take care of it.
	double exponent10Floor = floor( exponent10 );
	*rExponent = ( int ) exponent10Floor;
	
	return fraction2 * exp(( remainder10 - exponent10Floor + exponent10 )
						   * M_LN10 );
}

void printdouble( double num )
{
	if(PRINTEXP){
		double fraction;
		int exponent;

		fraction = frexp10(num,&exponent);
		
		printf("%.2lfe%d", fraction, exponent);
	}
	else{
		printf("%.4f", num);
	}
}

double get(Tableau * t, int row, int col){
    return t->data[row * t->cols + col];
}

void set(Tableau * t, int row, int col, double value){
    t->data[row * t->cols + col] = value;
}

void addrows(Tableau * t, int x, int y, int z){
	// z = x + y
	int i;
	for(i=0; i<t->cols; i++){
		if(i!=t->zcolumn && i!=t->basiscolumn){
			set(t,z,i,get(t,x,i)+get(t,y,i));
		}
	}
}

void subtractrows(Tableau * t, int x, int y, int z){
	// z = x - y
	int i;
	for(i=0; i<t->cols; i++){
		if(i!=t->zcolumn && i!=t->basiscolumn){
			set(t,z,i,get(t,x,i)-get(t,y,i));
		}
	}
}

void addscaledrows(Tableau * t, double scalar, int x, int y, int z){
	// z = x * scalar + y
	int i;
	double value;
	for(i=0; i<t->cols; i++){
		if(i!=t->zcolumn && i!=t->basiscolumn){
			value = get(t,x,i) * scalar;
			set(t,z,i,( value + get(t,y,i) ));
		}
	}
}


void scalerow(Tableau * t, double scalar, int z){
	// z = z * scalar
	int i;
	for(i=0; i<t->cols; i++){
		if(i!=t->zcolumn && i!=t->basiscolumn){
			set(t,z,i,get(t,z,i) * scalar);
		}
	}
}
/*
void divrow(Tableau * t, double scalar, int z){
	// z = z / scalar
	int i;
	for(i=0; i<t->cols; i++){
		if(i!=t->zcolumn && i!=t->basiscolumn){
			set(t,z,i,get(t,z,i) / scalar);
		}
	}
}
*/
void addobjectiverow(Tableau * t)
{
	int i,j;
	t->rows++;
	
	t->copydata = malloc(t->rows * t->cols * sizeof *t->data);
	
	for(i=0; i<t->rows; i++){
		for(j=0; j<t->cols; j++){
			if(i<t->rows-1){
				t->copydata[i * t->cols + j] = get(t,i,j);
			}
			else{
				t->copydata[i * t->cols + j] = 0;
			}
		}
	}
	
	free(t->data);
    t->data = t->copydata;
	
	// update zcolumn
	set(t,t->rows-1,t->zcolumn,1);
}

void deleterow(Tableau * t, int targetrow)
{
	int i,j;
	int rowskip = targetrow;
	
	t->rows--;
	if (targetrow < t->numconstraints){
		t->numconstraints--;
	}
	
	t->copydata = malloc(t->rows * t->cols * sizeof *t->data);
	
	for(i=0; i<t->rows; i++){
		for(j=0; j<t->cols; j++){
			if(i<rowskip){
				t->copydata[i * t->cols + j] = get(t,i,j);
			}
			else{
				t->copydata[i * t->cols + j] = get(t,i+1,j);
			}
		}
	}
	
	free(t->data);
    t->data = t->copydata;
}

void addartificialcolumn(Tableau * t, int targetrow)
{
	int i,j;
	int colskip = t->numvars;
	
	t->cols++;
	t->numvars++;
	t->numartificialvars++;
	t->zcolumn++;
	t->constcolumn++;
	t->basiscolumn++;
	
	t->copydata = malloc(t->rows * t->cols * sizeof *t->data);

	for(i=0; i<t->rows; i++){
		for(j=0; j<t->cols; j++){
			if(j<colskip){
				t->copydata[i * t->cols + j] = t->data[i * (t->cols-1) + j];
			}
			else if(j==colskip){
				if(i==targetrow){
					t->copydata[i * t->cols + j] = 1;
				}
				else{
					t->copydata[i * t->cols + j] = 0;
				}
			}
			else{
				t->copydata[i * t->cols + j] = t->data[i * (t->cols-1) + j - 1];
			}
		}
	}
		
	free(t->data);
    t->data = t->copydata;
	
	// update basis information
	set(t,targetrow,t->basiscolumn,colskip);
}

void deletecolumn(Tableau * t, int targetcolumn)
{
	int i,j;
	int colskip = targetcolumn;
	
	t->cols--;
	t->numvars--;
	t->numartificialvars--; // assume for now it is always an artificial column
	t->zcolumn--;
	t->constcolumn--;
	t->basiscolumn--;
	
	t->copydata = malloc(t->rows * t->cols * sizeof *t->data);

	for(i=0; i<t->rows; i++){
		for(j=0; j<t->cols; j++){
			if(j<colskip){
				t->copydata[i * t->cols + j] = t->data[i * (t->cols+1) + j];
			}
			else{
				t->copydata[i * t->cols + j] = t->data[i * (t->cols+1) + j + 1];
			}
		}
	}
		
	free(t->data);
    t->data = t->copydata;
}

bool isoptimal(Tableau * t)
{
    bool ret = true;
    int i = 0;
    
    printf( "Checking Optimality...  " );
    for(i=0; i<t->numvars; i++){
		// Check to see if we have a negative value in the objective row.
        if (flt_lt(get(t,t->rows-1,i),0)){
            ret = false;
        }
    }
    
    if (!ret){
        printf( "Is not optimal.\n" );
    } else {
        printf( "Is optimal.\n" );
    }
    
    return ret;
}

int findpivotcolumn(Tableau * t)
{
    int maxnegindex = -1;
    int maxneg = 0;
    int i = 0;
    
    printf( "Finding pivot column...  \n" );
    for(i=0; i<t->numvars; i++){
        printf( "\t%f\n",get(t,t->rows-1,i) );
        if (get(t,t->rows-1,i)<maxneg){
            //printf( "*" );
            maxnegindex = i;
            maxneg = get(t,t->rows-1,i);
        }
    }
    printf( "\n" );
        
    return maxnegindex;
}

int findpivotrow(Tableau * t, int pivotcolumn)
{
    int i = 0;
    double enteringvar, constvar, currentratio, minratio;
    int minratiorow;
	bool found = false;
    
    printf( "Finding pivot row...  \n" );
    for(i=0; i<t->numconstraints; i++){
        
		enteringvar = get(t,i,pivotcolumn);
		constvar = get(t,i,t->constcolumn);
		if(enteringvar <= 0){
			// do nothing because this ratio is of no value to us
			printf( "\tCurrent Ratio: Not Applicable\n");
		}
		else{
			currentratio = constvar/enteringvar;

			printf( "\tCurrent Ratio: %f\n",currentratio);
			
			if (!found){
				minratio = currentratio;
				minratiorow = i;
			} else {
				if (currentratio < minratio){
					minratio = currentratio;
					minratiorow = i;
				}
			}
			
			found = true;
		}
    }
    printf( "\n" );
    
	if (!found){
		return UNBOUNDED;
	}
	else{
		return minratiorow;
	}
}

void printtableau(Tableau * t)
{
    int i;
    int j;
    
    printf("\nCurrent Tableau:\n");
    
	// print indices
	printf("\t");
	for(j=0; j<t->cols; j++){
        printf("%d\t",j);
    }
	printf("\n");
	
	// print headers
    printf("\t");
	for(j=0; j<t->cols; j++){
        if(j<t->numrealvars){
            printf("x%d",j);
            printf("\t");
        }
        else if (j<(t->numrealvars+t->numslackvars)){
            printf("s%d",j-t->numrealvars);
            printf("\t");
        }
		else if (j<(t->numrealvars+t->numslackvars+t->numartificialvars)){
            printf("a%d",j-(t->numrealvars+t->numslackvars));
            printf("\t");
        }
        else if (j==t->zcolumn){
            printf("z\t");
        }
        else if (j==t->constcolumn){
            printf("const\t");
        }
		else if (j==t->basiscolumn){
            printf("basic\t");
        }
    }
    printf("\n");

    // print data
    for(i=0; i<t->rows; i++){
		if(i<t->numconstraints){
			printf("c%d\t",i);
		}
		else{
			printf("obj%d\t",i-t->numconstraints);
		}
	
        for(j=0; j<t->cols; j++){
            printdouble(get(t,i,j));
			printf("\t");
			//printf("%.2f\t",get(t,i,j));
        }
        printf("\n");
    }
    printf("\n");
}

void pivot(Tableau * t, int pivotrow, int pivotcolumn)
{
    printf("Pivoting...\n");
    
    int i;
    int j;
	double scalar;
    
    printf( "Pivot Row: %d\n", pivotrow );
	printf( "Pivot Column: %d\n", pivotcolumn );
	//printf( "A[%d,%d]\n",pivotrow+1,pivotcolumn+1);
        
    double current; // the current pivot value used for each non pivot row
    
    // update the basis information
    set(t,pivotrow,t->basiscolumn,pivotcolumn);

	// scale the pivot row
	scalar = 1 / get(t,pivotrow,pivotcolumn);
	printf("scalar: %.2f\n",scalar);
	scalerow(t,scalar,pivotrow);
    
    // find the rest of the new tableau
    for(i=0; i<t->rows; i++){
        if (i != pivotrow){
			// scale pivot row by the value in current row's pivot column
			// add that result to the current row
			scalar = -get(t,i,pivotcolumn);
			addscaledrows(t, scalar, pivotrow, i, i); // z = x * scalar + y
			//divrow(t,scalar,i);
			//addrows(t,pivotrow,i,i);
			//scalerow(t,scalar,i);
        }
    }
}

void printcurrentsolution(Tableau * t){
    printf("\n");
    int i;
    int j;
    int basisvar;
    
    for(i=0; i<t->numconstraints; i++){
		basisvar = get(t,i,t->basiscolumn);
		if(basisvar<t->numrealvars){
			printf("x%d: ",basisvar);
			printdouble(get(t,i,t->constcolumn));
			printf("\n");
		}
		else if(basisvar<t->numrealvars+t->numslackvars){
			printf("x%d: ",basisvar-(t->numrealvars));
			printdouble(get(t,i,t->constcolumn));
			printf("\n");
		}
		else if(basisvar<t->numrealvars+t->numslackvars+t->numartificialvars){
			printf("a%d: ",basisvar-(t->numrealvars+t->numslackvars));
			printdouble(get(t,i,t->constcolumn));
			printf("\n");
		}
    }
    
    printf("z: ");
	printdouble(get(t,t->rows-1,t->constcolumn));
	printf("\n\n");
}

void readtableaufromfile(Tableau * t, char filename[80]){
    FILE *ifp;
    char *mode = "r";
    int numrealvars, numslackvars, numconstraints;

    int i,j;
    double value;

    printf("Reading input file...\n");
    
    ifp = fopen(filename, mode);

    if (ifp == NULL) {
        fprintf(stderr, "Can't open input file!\n");
    }

    // get the initial number of variables and constraints
    fscanf(ifp, "%d %d %d", &numrealvars, &numslackvars, &numconstraints);
    printf("Number of Real Vars: %d... Number of Slack Vars: %d... Number of Constraints: %d...\n",numrealvars,numslackvars,numconstraints);
    
    t->rows=numconstraints+1;
    t->cols=numrealvars+numslackvars+3;
    t->numvars=numrealvars+numslackvars;
    t->numrealvars=numrealvars;
    t->numslackvars=numslackvars;
	t->numartificialvars=0;
    t->numconstraints=numconstraints;
	t->zcolumn=t->cols-3;
    t->constcolumn=t->cols-2;
	t->basiscolumn=t->cols-1;

	t->data = malloc(t->rows * t->cols * sizeof *t->data);
    
	// get the tableau
	for(i=0;i<t->rows;i++){
		for(j=0;j<t->cols;j++){
			fscanf(ifp, "%lf", &value);
			//printf("%lf\t", value);
			set(t,i,j,value);
		}
		//printf("\n", value);
	}
	
	fclose(ifp);
}

bool isinsidefeasibleregion(Tableau * t)
{
	double value;
	int lookup;
	bool isfeasible = true;
	int i;

	// look for negative basis vectors
	for(i=0; i<t->numconstraints; i++){
		value = get(t,i,t->basiscolumn);
		if(value < 0) isfeasible = false;
	}
	
	return isfeasible;
}

void createartificialtableau(Tableau * t)
{
	double value;
	int i;
	
	addobjectiverow(t);

	for(i=0; i<t->numconstraints; i++){
		value = get(t,i,t->basiscolumn);
		if(value < 0){
			printf("Adding artificial variable...\n");
			addartificialcolumn(t,i);
			// add -1 to objective row
			set(t,t->rows-1,t->numvars-1,-1);
		}
	}

	for(i=0; i<t->numconstraints; i++){
		if(get(t,i,t->basiscolumn) >= (t->numrealvars+t->numslackvars)){
			addrows(t,i,t->rows-1,t->rows-1);
		}
	}
	scalerow(t,-1,t->rows-1);
}

void destroyartificialtableau(Tableau * t)
{
	// copy the relevant data in the artificial tableau
	// data back into the original tableau
	
	printf("Unpacking artificial tableau...\n");
	
	int i, k;
	
	// delete 2nd objective row
	deleterow(t,t->rows-1);
	
	// delete all artificial columns
	for(i = t->numvars-1; i >= (t->numvars-t->numartificialvars); i--){
		deletecolumn(t,i);
	}
}

TabErr phasetwo(Tableau * t){
    int pivotcolumn;
    int pivotrow;
    
	while(!isoptimal(t)){
        printf ( "Press any key to continue . . ." );
		getchar();

        pivotcolumn = findpivotcolumn(t);
        pivotrow = findpivotrow(t,pivotcolumn);
		if(pivotrow<0){
			return pivotrow; // The pivot row has been loaded with the error
		}
        pivot(t, pivotrow, pivotcolumn);
        printtableau(t);
		printcurrentsolution(t);
    }
	
	return OK;
}

TabErr phaseone(Tableau * t){
	TabErr error = OK;
	
	if(!isinsidefeasibleregion(t)){
		// Create a temporary tableau with artificial variables
		// Solve to eliminate the artificial variables from the solution
		createartificialtableau(t);
		printf("Artificial tableau created...\n");
		printtableau(t);
		error = phasetwo(t);
		if(error < 0){
			return error;
		}
		if(flt_lt(get(t,t->rows-1,t->constcolumn),0)){
			return INFEASIBLE;
		}
		destroyartificialtableau(t);
		printtableau(t);
	}
	return OK;
}

void printerror(TabErr error)
{
	if (error==UNBOUNDED){
		printf("TabErr: UNBOUNDED\n");
	}
	else if (error == INFEASIBLE){
		printf("TabErr: INFEASIBLE\n");
	}
	else if (error == OK){
		printf("TabErr: OK\n");
	}
	else{
		printf("TabErr: UNKNOWN ERROR\n");
	}
}

int main(int argc, char*argv[])
{
    int pivotcolumn;
    int pivotrow;
	char filename[80];
	Tableau * t, root;
	t = &root; // init memory
	TabErr error;
	
	if (argc <= 1){
		printf("You need to supply a file name.\n");
		return 0;
	}
    
	strcpy(filename, argv[1]);
	readtableaufromfile(t, filename);
	printtableau(t);
	

	printf("\n\n======Phase 1:=====\n\n");
	error = phaseone(t);
	if(error!=0){
		printerror(error);
		return error;
	}
	
	printf("\n\n======Phase 2:=====\n\n");
	error = phasetwo(t);
	if(error!=0){
		printerror(error);
		return error;
	}
	
	printf("\nFinal Solution:\n");
    printcurrentsolution(t);


    free(t->data);

    return 0;
}
