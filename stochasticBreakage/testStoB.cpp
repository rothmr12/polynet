#include "stoB.h"
#include "stdio.h"

void run_single_bond(double F=0.0){
	
	printf("\nTest of stoB starting\n");

	// Instantiate stoB with debug on
	stoB test(true);

	// Set Number of bonds
	test.setN(1);
	test.printAll();

	// Set temperature, and timestep
	test.setTempDt(300,10);
	test.printAll();

	// Set bond breakage rate parameters
	test.setRate(0, 1e12, 1.2, -0.1e9);

	// Set Prob at which chain breaks
	test.setPChainMax(0.99);
	test.printAll();

	// Timestep until chain breaks
	printf("Timestepping ... ");
	fflush(stdout);
	while( !test.isBroken() ){
		test.timestep(F);
	}

	// Timestep once more to see message
	test.timestep(F);

	// Self explanatory
	printf("Done\n");
}

void run_two_bond(double F=0.0){
	
	printf("\nTest of stoB starting\n");

	// Instantiate stoB with debug on
	stoB test(true);

	// Set Number of bonds
	test.setN(2);
	test.printAll();

	// Set temperature, and timestep
	test.setTempDt(300,10);
	test.printAll();

	// Set bond breakage rate parameters
	test.setRate(0, 1e12, 1.2, -0.1e9);

	// Set Prob at which chain breaks
	test.setPChainMax(0.99);
	test.printAll();

	// Timestep until chain breaks
	printf("Timestepping ... ");
	fflush(stdout);
	while( !test.isBroken() ){
		test.timestep(F);
	}

	// Timestep once more to see message
	test.timestep(F);

	// Self explanatory
	printf("Done\n");
}


int main(){

	run_single_bond(1.0e-9);
	rename("timestepping.data","single-bond-F-1nN");
	
	run_single_bond();
	rename("timestepping.data","single-bond-F-0nN");

	run_two_bond(1.0e-9);
	rename("timestepping.data","two-bond-F-1nN");

	run_two_bond();
	rename("timestepping.data","two-bond-F-0nN");

	// Force cycling on single bond

	// Multiple bonds
}
