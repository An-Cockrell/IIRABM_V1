	const int xDim=101; //grid size in x-dimension
	const int yDim=101; //grid size in y-dimension
	
	const float antibioticMultiplier=0.2; //controls proportion of bacteria killed by antibiotics; 
	//after 1 application of antibiotics, the amount of bacteria that remains is: antibioticMultiplier*(curent level of bacteria)
	const int parameterInput=0; //number of parameter sets to be read in (as opposed to generated by the below start,max,increment scheme

	const int numTimeSteps=12500; //5760=28 days
      const int numWriteSteps=12500;	//number of time steps to write out
	
	const int numReseeds=0; //how many times to reseed an individual parameter set
	const int reseedTime=505; //time point at which the simulation gets reseeded

	const int inj_numberStart=33;  //initial injury size
	const int inj_numberMax=33; //maximum injury size
	const int inj_numberIncrement=5; //injury step size
	
	const int numInfectRepeatStart=2;  //initial microbial invasiveness
	const int numInfectRepeatIncrement=1; //invasiveness step size
	const int numInfectRepeatMax=2; //maximum invasiveness

	const float oxyHealStart=0.1; //initial host resilience
	const float oxyHealIncrement=0.05; //host resilience step size
	const float oxyHealMax=0.1; //maximum host resilience
	
	const int infectSpreadStart=5; //initial toxigenesis
	const int infectSpreadIncrement=1; //toxigenesis step size
	const int infectSpreadMax=5; //maximum toxigenesis

	const int numRecurInjStart=2; //initial environmental toxicity
	const int numRecurInjIncrement=1; //environmental toxicity step size
	const int numRecurInjMax=2; //maximum environmental toxicity

	const int seedStart=2; //initial random number seed
	const int seedMax=2; //maximum random number seed
	const int seedIncrement=1; //random number seed step size

	const int injuryStep=205; //time for first recurring injury (environmental toxicity)

	const int cellCapacity=28; //maximum population per grid cell

	const int movie=0; //set to 1 to write out complete simulation state at every time step
	const int movieStep=205; //time step to begin writing out entire simulation state
	const int oxyInfectFinal=0; //set to 1 to write out the final outcomes
	const int dailyOutput=0; //set to 1 to write aggregate cytokine values daily
	const int writeEverything=1; //set to 1 to write aggregate cytokine values at every time step
	
	const string fileName="Outcome_I2_H1_T5_E2_inj33_100seed_int9_99.dat"; //output file name
	const string fileNameGA="GA_I2_H1_T5_E2_inj33_10rep_5II_length99_longrun.dat"; //file name for GA results

//// GA Parameters
	const int GA=0; //turn GA on/off with 1/0
	const int GAseed=1; //random number seed for GA
	const int reseedPostInt=1; //set to 1 to reseed training run after start of intervention
	const int interventionTimeStart=505; //time step to start intervention
	const int interventionTimeSeperation=100; //number of time steps between interventions
	const int tournamentSize=2; //tournament selection parameter
	const int ga_popSize=1000;	//initial population size
	const int maxIter=2000; //maximum number of generations
	const int genomeLength=60; //length of GA-genome, should be number of interventions*12
	const int numIndependentInterventions=5; //should equal genomeLength/12
	const int mutationRate=5; //out of 10
	const int replicateMult=10; //number of times to reseed after 1st intervention
	const int interventionLength=99; //how long each intervention lasts
	
//Intervention Test Parameters

	const int testIntervention=1; //test a specific intervention according to the below parameters
	const int testLength=36; //length of intervention vector to be tested
	const int numTests=3; //number of independent interventions to be tested
	const int testStart=205; //time step at which to start intervention test
	const string intFile="intervention8_truncatedat3.dat"; //filename of intervention to be tested
	

