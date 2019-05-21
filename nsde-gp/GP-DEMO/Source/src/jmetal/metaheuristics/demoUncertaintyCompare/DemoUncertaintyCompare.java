//  MetaDEMO.java
//
//  Author:
//       Miha Mlakar <miha.mlakar@ijs.si>
//       Dejan Petelin <dejan.petelin@ijs.si>
//
//  Copyright (c) 2012 Miha Mlakar, Dejan Petelin
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
package jmetal.metaheuristics.demoUncertaintyCompare;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.TimeUnit;

import com.sun.corba.se.spi.servicecontext.UEInfoServiceContext;

import jmetal.core.Algorithm;
import jmetal.core.Model;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.init.Params;
import jmetal.models.SPGPm;
import jmetal.models.SPGPm2;
import jmetal.qualityIndicator.Hypervolume;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.util.RankingUnderUncertanty;
import jmetal.util.RankingUnderUncertantyZaStPrimerjav;
import jmetal.util.Spea2Fitness;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.DominanceComparatorUnderUncertanty;

/**
 * This class implements the MetaDEMO algorithm. This algorithm is an adaptation
 * of the DEMO defined by T. Tusar in:
 * T. Tusar. ....
 */
public class DemoUncertaintyCompare extends Algorithm{
  
  /**
   * Stores the model used for estimation of objective values.
   */
  private Model modelLocal_;
  private Model modelGlobal_;
  
  private int logMode;
  private int genMode;
  private int evaluations;
  private int iterations;
  private int frontGen;
  private int frontFormat;
  private int frontMode;
  private int populationSize;
  private int maxEvaluations;
  private int maxIterations;
  private int selectionProcedure;

  private String genFileName;
  private String logFileName;
  private String frontFileName;
  private Operator selectionOperator;
  private Operator crossoverOperator;
  
	private int currentGenerationOfApproximateEvaluations;
	private int currentGenerationOfExactEvaluations;
	private int additionalEvaluations;
	private int exactEvaluations;
	private int exactEvaluationsDuringEvolutionProcess;
	
  private int numberOfApproximateEvaluationsOfGenetaions;
  private int numberOfExactEvaluationsOfGenetaions;
  private int afterWhatnumberOfGenerationsResetOfBaseVectors;
  private int afterWhatnumberOfGenerationsUpdateTheModel;
  private int useGenerationsForDeterminingTheUpdates;
  private int afterWhatnumberOfExactEvalSolutionsUpdateTheModel;
  //local model
  private int afterWhatnumberOfGenerationsUpdateTheLocalModel;
  private int useGenerationsForDeterminingTheLocalModelUpdates;
  private int afterWhatnumberOfExactEvalSolutionsUpdateTheLocalModel;
  int numberOfExactlyEvalSolutionsTakenForLocalModel;
  
  int steviloVsehPrimerjav = 0;
  int steviloNOBBCorrect = 0;
  int steviloNOBBIncorrect = 0;
  
  private boolean useExactEvaluation;
  
  private SolutionSet setOfExactlyEvaluatedSolutions;
  private double maximumAllowedVariance [];
  private boolean useTwoModels;
  
  private SolutionSet setOfSolutionsForLocalModel;
  
  private int limitNumberForExactEvaluations;
  
  String[] whichModelIsBetter;
  
  //private SolutionSet population;
  /**
  * Constructor
  * @param problem Problem to solve
  * @param model Model used for estimation
  */
  public DemoUncertaintyCompare(Problem problem, Model globalModel) {
    this(problem, globalModel, (Model)null);
  } // MetaDEMO
  
  /**
  * Constructor
  * @param problem Problem to solve
  * @param model Model used for estimation
  */
  public DemoUncertaintyCompare(Problem problem, Model globalModel, Model localModel) {
    super (problem);

    if(localModel != null) {
    	useTwoModels = true;
    	modelLocal_ = localModel;
    	modelGlobal_ = globalModel;
    }
    else {
    	useTwoModels = false;
    	modelGlobal_ = globalModel;
    	modelLocal_ = null;
    }
    
    initialize();
  } // MetaDEMO

  /**
   * Constructor
   * @param problem Problem to solve
   * @param model Model used for estimation
   * @param properties Properties of the algorithm
   */
   public DemoUncertaintyCompare(Problem problem, Model globalModel, Properties properties) {
     this(problem, globalModel, null, properties);
   } // MetaDEMO

 /**
  * Constructor
  * @param problem Problem to solve
  * @param model Model used for estimation
  * @param properties Properties of the algorithm
  */
  public DemoUncertaintyCompare(Problem problem, Model globalModel, Model localModel, Properties properties) {
    super (problem, properties);
    
    if(localModel != null) {	
    	useTwoModels = true;
    	modelLocal_ = localModel;
    	modelGlobal_ = globalModel;
    }
    else {
    	useTwoModels = false;
    	modelGlobal_ = globalModel;
    	modelLocal_ = null;
    }
    
    initialize();
  } // MetaDEMO
  
  private void initialize() {
    frontFileName      = getParameter("Output.frontFileName" );
    logFileName        = getParameter("Output.logFileName");
    genFileName        = getParameter("Output.genFileName");
    genMode            = Integer.parseInt(getParameter("Output.genMode"));
    logMode            = Integer.parseInt(getParameter("Output.logMode"));    
    
    frontGen           = Integer.parseInt(getParameter("Output.frontGen"));
    frontFormat        = Integer.parseInt(getParameter("Output.frontFormat"));
    frontMode          = Integer.parseInt(getParameter("Output.frontMode"));

    populationSize     = Integer.parseInt(getParameter("Algorithm.populationSize"));
    selectionProcedure = Integer.parseInt(getParameter("Algorithm.selectionProcedure"));
    maxEvaluations     = Integer.parseInt(getParameter("Algorithm.maxEvaluations"));
    maxIterations      = maxEvaluations / populationSize + (int)Math.ceil(maxEvaluations%populationSize);
    
    selectionOperator  = Params.getSelectionOperator();		//operators_.get("selection");   
    crossoverOperator  = Params.getCrossoverOperator();			//operators_.get("crossover");
    
  	currentGenerationOfApproximateEvaluations = 0;
  	currentGenerationOfExactEvaluations = 0;
  	additionalEvaluations = 0;
  	exactEvaluations = 0;
  	exactEvaluationsDuringEvolutionProcess = 0;
  	
    numberOfExactEvaluationsOfGenetaions = Params.getNumberOfExactEvaluations();
    numberOfApproximateEvaluationsOfGenetaions = Params.getNumberOfApproximateEvaluations();
    afterWhatnumberOfGenerationsResetOfBaseVectors     = Integer.parseInt(getParameter("Model.afterWhatnumberOfGenerationsResetOfBaseVectors")); //Integer.valueOf(properties.getProperty("Model.windowSize"))
    afterWhatnumberOfGenerationsUpdateTheModel = Integer.parseInt(getParameter("Model.afterWhatnumberOfGenerationsUpdateTheModel"));
    afterWhatnumberOfGenerationsUpdateTheModel = 100000000;
    useGenerationsForDeterminingTheUpdates = Integer.parseInt(getParameter("Model.useGenerationsForDeterminingTheUpdates")); 
    afterWhatnumberOfExactEvalSolutionsUpdateTheModel = Integer.parseInt(getParameter("Model.afterWhatnumberOfExactEvalSolutionsUpdateTheModel")); 
    afterWhatnumberOfExactEvalSolutionsUpdateTheModel = 1000000000;
    
    //local model
    afterWhatnumberOfGenerationsUpdateTheLocalModel = Integer.parseInt(getParameter("ModelLocal.afterWhatnumberOfGenerationsUpdateTheLocalModel"));
    useGenerationsForDeterminingTheLocalModelUpdates = Integer.parseInt(getParameter("ModelLocal.useGenerationsForDeterminingTheLocalModelUpdates")); 
    afterWhatnumberOfExactEvalSolutionsUpdateTheLocalModel = Integer.parseInt(getParameter("ModelLocal.afterWhatnumberOfExactEvalSolutionsUpdateTheLocalModel")); 
    numberOfExactlyEvalSolutionsTakenForLocalModel = 0;
  	
    limitNumberForExactEvaluations= Integer.parseInt(getParameter("Algorithm.limitNumberForExactEvaluations"));
    
    useExactEvaluation = true;        
     
    setOfSolutionsForLocalModel = new SolutionSet();
    setOfSolutionsForLocalModel.setMaxSize(Integer.MAX_VALUE);
        
  }
  
  /**   
   * Runs of the MetaDEMO algorithm.
   * @return a <code>SolutionSet</code> that is a set of non dominated solutions
   * as a result of the algorithm execution.
   * @throws JMException 
   */  
  public SolutionSet execute() throws JMException, ClassNotFoundException {    
    int[] evolutionStat;
    int[] state;
    SolutionSet population;
    setOfExactlyEvaluatedSolutions = new SolutionSet(Integer.MAX_VALUE);
    maximumAllowedVariance = new double[Params.getNumberOfObjectives()];	
    for(int i=0; i<Params.getNumberOfObjectives();i++)
    	maximumAllowedVariance[i] = -1;	//???
    
    //delete content of "dejan.txt" file
    deleteDejan();    
    
    // set variable values
    iterations = 0;
    evaluations = 0;
    population = new SolutionSet(populationSize * 2);        
    currentGenerationOfApproximateEvaluations = 0;
    currentGenerationOfExactEvaluations = 0;
    additionalEvaluations = 0;
    
    createInitialPopulation(population);
    exactEvaluations = populationSize;
   
    // first iteration is initial population
    iterations++;             

    printFront(population);
    if (genMode != 0) {
    	state = new int[]{-1, -1, -1};
    	printGeneration(population, state, exactEvaluations, exactEvaluations);
    }
         
    //setOfExactlyEvaluatedSolutions.printVariablesAndObjectivesToFile(logFileName + "dejan.txt");
    
    modelGlobal_.update(setOfExactlyEvaluatedSolutions);
    if(useTwoModels)
    {
    	modelLocal_.update(setOfExactlyEvaluatedSolutions);
    	
    	for(int i =0; i<setOfExactlyEvaluatedSolutions.size();i++)
    		setOfSolutionsForLocalModel.add(setOfExactlyEvaluatedSolutions.get(i));
    	
    	numberOfExactlyEvalSolutionsTakenForLocalModel = setOfExactlyEvaluatedSolutions.size();
    }
    
    // OPTIMIZATION PROCESS
    while (iterations < maxIterations) {
    	
    	if(iterations>1)
    	{
    	try {
	      updateSurrogateModels(population);
      } catch (InterruptedException e) {  e.printStackTrace(); } catch (ExecutionException e) { e.printStackTrace();  }
    	}
    	
    	evolutionStat = doEvolution(population);
    	exactEvaluationsDuringEvolutionProcess = setOfExactlyEvaluatedSolutions.size();    	    	
    	
    	population = doSelection(population);
    	    
    	/*if(useTwoModels)
    	{
    		setOfSolutionsForLocalModel.clear();
    		for(int i=0; i<population.size(); i++)
    		{
    			if(population.get(i).isExactllyEvaluated() == true)
    				setOfSolutionsForLocalModel.add(population.get(i));
    		}
    		modelLocal_.update(setOfExactlyEvaluatedSolutions);
    	} */
    	
    	System.out.println("All evaluations: " + evaluations);
    	System.out.println("Exact evaluations: " + exactEvaluations);
    	
    	setOfExactlyEvaluatedSolutions.printVariablesAndObjectivesToFile(logFileName + "dejan.txt");   	
    	
    	iterations++;
    	
      printFront(population);
      if (genMode != 0) {
      	state = new int[]{evolutionStat[0], evolutionStat[1], evolutionStat[2]};
      	printGeneration(population, state, exactEvaluations, exactEvaluationsDuringEvolutionProcess);
      }
      
      
    } // while
  
    try{
    	if(modelLocal_ != null)
    		modelLocal_.closeRConnection();
    	
    	modelGlobal_.closeRConnection();
    }
    catch (Exception e) {
		}
    
    
    
    
    population.printToFile("D:\\Doktorat\\Metka\\app\\jMetalCI\\primerjava.txt", "Vse skupaj: " + String.valueOf(steviloVsehPrimerjav)+"   NOBB correct: "+ steviloNOBBCorrect + ",   NOBB incorrect " + steviloNOBBIncorrect);
    
    
    
    
    //recalculate last front so that all solutions are exactly evaluated
    additionalEvaluations = 0;
    Ranking ranking = new Ranking(population);
    for(int i=0; i< ranking.getSubfront(0).size(); i++)
    {
    	if(ranking.getSubfront(0).get(i).isExactllyEvaluated() == false)
    	{
    		problem_.evaluate(ranking.getSubfront(0).get(i));
    		ranking.getSubfront(0).get(i).setStandardDevianceToZero();
    		problem_.evaluateConstraints(ranking.getSubfront(0).get(i));    		
    		ranking.getSubfront(0).get(i).setExactllyEvaluated(true);
    		additionalEvaluations++;    		
    	}
    }    	
    
    Hypervolume hyp = new Hypervolume();
    population.printToFile(genFileName, "Additional evaluations: " + additionalEvaluations);
    population.printToFile(genFileName, "Final hypervolume: " + hyp.hypervolume(ranking.getSubfront(0).writeObjectivesToMatrix(), Params.getFrontForHypervolume(), problem_.getNumberOfObjectives()));
    
    population.printToFile(frontFileName, "---------------------------------------------------------------------------------------------------");
    
    if (frontFormat == 1) // decision and objective vector
      ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
    else if (frontFormat == 0) // only objective vector
      ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);  
    
    // Return the first non-dominated front
    System.out.println("Dodatni izracuni " + additionalEvaluations);
    //return new RankingUnderUncertanty(population, problem_).getSubfront(0);
    return new Ranking(population).getSubfront(0);
  } // execute

  
  private void updateSurrogateModels(SolutionSet population) throws ClassNotFoundException, JMException, InterruptedException, ExecutionException {

  	
  	if(setOfExactlyEvaluatedSolutions.size() != 0)
  	{   
  		int paralel = 0;
  		
  		if(paralel ==1)
  		{
  		List<Future> futuresList = new ArrayList(); 		
  		String taskResult;
  		ForkJoinPool fjPool = new ForkJoinPool(2);
  		Future global = null;
  		Future local = null;
  		
  		try 
  		{
  				if(iterations%1 == 0)
  				{
  					global = (fjPool.submit(new UpdateGlobalModel(modelGlobal_, setOfExactlyEvaluatedSolutions, iterations%afterWhatnumberOfGenerationsResetOfBaseVectors == 0)));  					
  					setOfExactlyEvaluatedSolutions.clear();
  	  		}
  				
  				
  		    if(useTwoModels)
  		  	{
  		  		setOfSolutionsForLocalModel.clear();
  		  		for(int i=0; i<population.size(); i++)
  		  		{
  		  			if(population.get(i).isExactllyEvaluated() == true)
  		  				setOfSolutionsForLocalModel.add(population.get(i));
  		  		}  		  		
  		  		if(setOfSolutionsForLocalModel.size() != 0)
  		  		{
  		  			local = (fjPool.submit(new UpdateLocalModel(modelLocal_, setOfSolutionsForLocalModel)));    					
  		  		}
  		  	}
  		    
  		    taskResult = (String)local.get();
          System.out.println("result ForkJoin "+taskResult);
  		    taskResult = (String)global.get();
          System.out.println("result ForkJoin "+taskResult);
  		}
  		finally
  		{}

  		System.out.println("Modeli posodobljeni!");
  		
  	}
  	
 		else if(paralel == 2)
 		{
  	if(setOfExactlyEvaluatedSolutions.size() != 0)
  	{   
  		 		
  		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
  		try 
  		{
  				if(iterations%1 == 0)
  				{
  					exec.submit(new UpdateModelParalel(modelGlobal_, modelGlobal_.getName(), setOfExactlyEvaluatedSolutions, iterations%afterWhatnumberOfGenerationsResetOfBaseVectors == 0));
  					setOfExactlyEvaluatedSolutions.clear();
  	  		}
  				
  		    if(useTwoModels)
  		  	{
  		  		setOfSolutionsForLocalModel.clear();
  		  		for(int i=0; i<population.size(); i++)
  		  		{
  		  			if(population.get(i).isExactllyEvaluated() == true)
  		  				setOfSolutionsForLocalModel.add(population.get(i));
  		  		}
  		  		
  		  		if(setOfSolutionsForLocalModel.size() != 0)
  		  			exec.submit(new UpdateModelParalel(modelLocal_, "localModel", setOfSolutionsForLocalModel, false));
  		  	}
  		    
  		} 
  		finally {
  		    exec.shutdown();
  		}
  		
  		while (! exec.isTerminated()) {
  		  try { 
  		  	exec.awaitTermination(1000, TimeUnit.MILLISECONDS);
  		  }
  		  catch (InterruptedException e) {
  		     
  		  }
  		}
  		System.out.println("Modeli posodobljeni!");
  		
  	}
 		}
  		
  	else	//without paralelism
  	{
  		boolean firstCondition = (((iterations-1)%afterWhatnumberOfGenerationsUpdateTheModel == 0) && useGenerationsForDeterminingTheUpdates == 1); //-1 becouse after initial (first) generation we always do the update
  		boolean secondCondition = ((setOfExactlyEvaluatedSolutions.size() >= afterWhatnumberOfExactEvalSolutionsUpdateTheModel) && useGenerationsForDeterminingTheUpdates == 0);
  		
  		if(firstCondition || secondCondition)   
  		{
	  		if(modelGlobal_.getName() == "SPGP")
	  			((SPGPm)(modelGlobal_)).update(setOfExactlyEvaluatedSolutions, iterations%afterWhatnumberOfGenerationsResetOfBaseVectors == 0);
	  		else
	  			modelGlobal_.update(setOfExactlyEvaluatedSolutions); 
	  		
	    	setOfExactlyEvaluatedSolutions.clear();
  		}
  		
  		
		if(useTwoModels)
		{		  	
					
			for(int i=0; i<population.size(); i++)
  		{
  			if(population.get(i).isExactllyEvaluated() == true)
  			{
  				if(!isSolutionAlreadyInTheDataSet(population.get(i), setOfSolutionsForLocalModel, Integer.parseInt(getParameter("ModelLocal.activeSet"))))
  				{
  					setOfSolutionsForLocalModel.add(population.get(i));  					
  					//if(setOfSolutionsForLocalModel.size()> Integer.parseInt(getParameter("ModelLocal.activeSet")))
  					//	setOfSolutionsForLocalModel.remove(0);
  				}
  			}
  		}
		
			boolean firstConditionLocalModel = (((iterations-1)%afterWhatnumberOfGenerationsUpdateTheLocalModel == 0) && useGenerationsForDeterminingTheLocalModelUpdates == 1); //-1 becouse after initial (first) generation we always do the update
  		boolean secondConditionLocalModel = (((setOfSolutionsForLocalModel.size() - numberOfExactlyEvalSolutionsTakenForLocalModel) >= afterWhatnumberOfExactEvalSolutionsUpdateTheLocalModel) && useGenerationsForDeterminingTheLocalModelUpdates == 0);
  		
  		if(firstConditionLocalModel || secondConditionLocalModel)   
  		{
  			System.out.println("Local zacetek!");
  			   			
  			SolutionSet updateLocalModelDataSet = new SolutionSet(Integer.parseInt(getParameter("ModelLocal.activeSet")));
  			
  			for(int i=numberOfExactlyEvalSolutionsTakenForLocalModel; i<setOfSolutionsForLocalModel.size();i++)
  				updateLocalModelDataSet.add(setOfSolutionsForLocalModel.get(i));
  				
  			modelLocal_.update(updateLocalModelDataSet);
  			
  			//at the beginning this number is equal to the population size (first all exactli eval solutions)
  			numberOfExactlyEvalSolutionsTakenForLocalModel = setOfSolutionsForLocalModel.size();	// - numberOfExactlyEvalSolutionsTakenForLocalModel;
  			//setOfSolutionsForLocalModel.clear();
  			
  			System.out.println("Local posodobljeni!");
  		}
 				  		
		  } 
		
  	}
		  
  	}
  	
  	
  }

	private boolean isSolutionAlreadyInTheDataSet(Solution solution,  SolutionSet setOfSolutions, int activeSet) throws JMException {

		boolean alreadyInside = false;
		int indexDeterminingTheBeginingOfTheComparison;		//to achieve that we do not compare with solutions that are not "active" any more becouse of the activeSet limitation
		
		if(setOfSolutions.size() > activeSet)
			indexDeterminingTheBeginingOfTheComparison = setOfSolutions.size() - activeSet;
		
		else
			indexDeterminingTheBeginingOfTheComparison = 0;
		
		for(int i=indexDeterminingTheBeginingOfTheComparison; i< setOfSolutions.size(); i++)
			if(areSolutionsTheSame(solution, setOfSolutions.get(i)))
			{
				alreadyInside = true;
				break;
			}
		
		return alreadyInside;
  }

	private void deleteDejan() {
  	try {
    	BufferedWriter bw;    
	    bw = new BufferedWriter(new FileWriter(logFileName + "dejan.txt"));    
    bw.write("");
    bw.flush();
    bw.close();
    } catch (IOException e) {
	    e.printStackTrace();
    }
  }

	private void printFront(SolutionSet population) throws JMException {
  	Ranking ranking;
  	
    if (iterations % frontGen == 0){
      if (frontMode != 0) {
        //ranking = new RankingUnderUncertanty(population, problem_);
      	ranking = new Ranking(population);
        if (frontFormat == 1) // decision and objective vector
          ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
        else if (frontFormat == 0) // only objective vector
          ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);  
      }
    }
  }
  
  
  private void printGeneration(SolutionSet population, int[] state, int exactEvaluation, int exactEvaluationsDuringEvolutionProcess) {
    int numberOfParentBetter;
    int numberOfIncomparable;
    int numberOfChildBetter;

  	numberOfChildBetter = state[0];
    numberOfParentBetter = state[1];
    numberOfIncomparable = state[2];
  	
    population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, numberOfChildBetter, numberOfParentBetter, numberOfIncomparable, exactEvaluation, exactEvaluationsDuringEvolutionProcess);    
  }
  
  private void createInitialPopulation(SolutionSet population) throws ClassNotFoundException, JMException {
  	Solution newSolution;
  	
  	// Create the initial solutionSet
    for (int i = 0; i < populationSize; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      newSolution.setStandardDevianceToZero();
      problem_.evaluateConstraints(newSolution);      
      newSolution.setExactllyEvaluated(true);
      setOfExactlyEvaluatedSolutions.add(newSolution);
      evaluations++;
      if (logMode != 0)
        newSolution.printVariablesandObjectivesAndStandardDeviationToFile(logFileName, evaluations, newSolution.isExactllyEvaluated(), newSolution.getstandardDeviances(), whichModelIsBetter, null, null, null);
      
      population.add(newSolution);
    } // for
  }
  
  private int reevaluateFrontSolutions(SolutionSet population) throws JMException {
    int tempAdditionalEvaluations;
    int additionalEvaluations;
    Ranking ranking;
    
    additionalEvaluations = 0;
    do {
      tempAdditionalEvaluations = 0;
      ranking = new Ranking(population);
      for (int i = 0; i < ranking.getSubfront(0).size(); i++) {
          if (ranking.getSubfront(0).get(i).isExactllyEvaluated() == false) { // this means that it was calculated by NN
            additionalEvaluations++;
            tempAdditionalEvaluations++;
            problem_.evaluate(ranking.getSubfront(0).get(i));
            ranking.getSubfront(0).get(i).setStandardDevianceToZero();
            problem_.evaluateConstraints(ranking.getSubfront(0).get(i));            
            ranking.getSubfront(0).get(i).setExactllyEvaluated(true);
          }                    
        }  
      
    } while (tempAdditionalEvaluations > 0);
    return additionalEvaluations;
  }
  
  private SolutionSet doSelection(SolutionSet population) throws JMException, ClassNotFoundException {  	
  	
  	//additionalEvaluations += reevaluateFrontSolutions(population);

    if (selectionProcedure == 0) { // NSGA 2
      // Ranking the offspring population      
      population = environmentalSelectionNSGAIIWithUncertanty(population);
    }
    else if (selectionProcedure == 1) { // IBEA
      while (population.size() > populationSize) { // offspringPopulation
        //removeWorst(offspringPopulation);
      }
    }
    else if (selectionProcedure == 2) { // SPEA 2
      population=environmentalSelectionSPEA2(population, populationSize);  // offspringPopulation
    }
                  
    // we get population of reduced (half) size after the environmental selection               
    // so we double it back. 
    population.setMaxSize(populationSize*2);
    
    return population;
  }

  private int[] doEvolution(SolutionSet population) throws JMException, ClassNotFoundException {
  	int numberOfChildBetter;
    int numberOfParentBetter;
    int numberOfIncomparable;
    Solution parents[];    
    int result;

    numberOfChildBetter = 0;
    numberOfParentBetter = 0;
    numberOfIncomparable = 0;
    double[] solutionApproximatedValue = new double[Params.getNumberOfObjectives()];
    double[] solutionConfidenceInterval = new double[Params.getNumberOfObjectives()];
    
    String[] allObjectivesAndConfidences;
    
      
    
    
    try
    {

    
    for (int i = 0; i < populationSize; i++) {
      // Obtain parents. Two parameters are required: the population and the 
      //                 index of the current individual
      parents = (Solution [])selectionOperator.execute(new Object[] { population, i} ); // GENERIRANJE KANDIDATA

      Solution child ;
      // Crossover. Two parameters are required: the current individual and the 
      //            array of parents
      child = (Solution)crossoverOperator.execute(new Object[]{population.get(i), parents});   
      
      // Candidate evaluation  popravek
      /*allObjectivesAndConfidences = approximateSolution(child);
      
      for(int j=0; j<Params.getNumberOfObjectives();j++)
      {
      	solutionApproximatedValue[j] = child.getObjective(j);
      	solutionConfidenceInterval[j] = child.getstandardDeviance(j);
      }
      //model_.evaluate(child);
      problem_.evaluateConstraints(child);
      */
      
      useExactEvaluation = false;
      //if the variance is too high we need to calculate to explore search space
      for(int j=0; j< Params.getNumberOfObjectives(); j++)     
      	if(child.getstandardDeviance(j)> maximumAllowedVariance[j])
      		useExactEvaluation = true;
          
      if (useExactEvaluation) {          		
        problem_.evaluate(child) ;
        child.setStandardDevianceToZero();
        problem_.evaluateConstraints(child);        
        child.setExactllyEvaluated(true);
        setOfExactlyEvaluatedSolutions.add(child);
        exactEvaluations ++;
      }
      
      
      evaluations++ ;
      //if (logMode != 0)
      //  child.printVariablesandObjectivesAndStandardDeviationToFile(logFileName, evaluations, child.isExactllyEvaluated(), child.getstandardDeviances());
                         
      // Dominance test
      if (areSolutionsTheSame(population.get(i), child)) {
      	steviloVsehPrimerjav ++;
      	if(population.get(i).isExactllyEvaluated() == false)
      	{
      		if(child.isExactllyEvaluated() == true)
      		{
      			population.replace(i, child);
      			numberOfIncomparable++;
      		}     	
      		else
      		{
      			population.replace(i, child);
      			numberOfIncomparable++;
      		}
      	}
      	else
      		numberOfIncomparable++;
        
      }
      
      else // the two solutions are not same 
      {
      	//pravilen rezultat, ker sta obe rešitvi ekzaktni
      	result = compareAndEvaluateSolutionsUnderDominance(population.get(i), child, population);
      	steviloVsehPrimerjav ++ ;
      	
      	int rezultatAproksimacije = primerjajApproksimiraniResitviNOBB(population.get(i), child);
      	
      	if(result == rezultatAproksimacije)
      		steviloNOBBCorrect ++;
      	else
      		steviloNOBBIncorrect++;
      	
      	      	
        if (result == -1) { // Solution i dominates child
          numberOfParentBetter++;                    
        } // if
        else if (result == 1) { // child dominates
          population.replace(i, child);
          numberOfChildBetter++;                                            
        } // else if
        else if (result == 0){ // the two solutions are non-dominated
          population.add(child);
          numberOfIncomparable++;
        }
        	       
      }
      
      //if (logMode != 0)
        //child.printVariablesandObjectivesAndStandardDeviationToFile(logFileName, evaluations, child.isExactllyEvaluated(), child.getstandardDeviances(), whichModelIsBetter, solutionApproximatedValue, solutionConfidenceInterval, allObjectivesAndConfidences);
      
    } // for

  	numberOfParentBetter = 0;
  
    
    }

    catch (Exception e) {
    	numberOfParentBetter = 0;
		}

    return new int[]{numberOfChildBetter, numberOfParentBetter, numberOfIncomparable};
  }
  
  private int primerjajApproksimiraniResitviNOBB(Solution solution, Solution child) throws ClassNotFoundException, JMException {
	  
  	Solution s1 = new Solution(solution);
  	Solution s2 = new Solution(child);
  	
  	approximateSolution(s1);
  	s1.setStandardDevianceToZero();
  	s1.setExactllyEvaluated(true);
  	approximateSolution(s2);
  	s2.setStandardDevianceToZero();
  	s2.setExactllyEvaluated(true);
  	
  	return compareAndEvaluateSolutionsUnderDominance(s1, s2, new SolutionSet());  	
  }

	private String[] approximateSolution(Solution child) throws ClassNotFoundException, JMException {

  	double distanceToFisibility;
  	double standardDevianceOfGlobalModel;
  	double[] confidenceIntervalsOfGlobalModel = new double[child.numberOfObjectives()];
  	double[] confidenceIntervalsOfLocalModel = new double[child.numberOfObjectives()];
  	double[] approximatedObjectivesOfGlobalModel = new double[child.numberOfObjectives()];
  	double[] approximatedObjectivesOfLocalModel = new double[child.numberOfObjectives()];
  	String[] allObjectivesAndConfidences = new String[child.numberOfObjectives() * 2 * 2  + 4];	//obj + confidence, 2 models, 4 labels
  	
  	Solution globalModelSolution;

  	whichModelIsBetter = new String[child.numberOfObjectives()]; 
  	for(int i=0; i<whichModelIsBetter.length;i++)
  		whichModelIsBetter[i] = "OneModel";
  	
  	modelGlobal_.evaluate(child);
  	problem_.evaluateConstraints(child);
  	//write data for log file
  	for(int i=0; i<child.numberOfObjectives(); i++)
		{
  		confidenceIntervalsOfGlobalModel[i] = child.getstandardDeviance(i);
  		approximatedObjectivesOfGlobalModel[i] = child.getObjective(i);
		}
  	
  	child.setExactllyEvaluated(false);
  	
  	if(useTwoModels)
  	{
  		globalModelSolution = new Solution(child);
  		
  		modelLocal_.evaluate(child);
  		problem_.evaluateConstraints(child);
  		//write data for log file
  		for(int i=0; i<child.numberOfObjectives(); i++)
  		{
    		confidenceIntervalsOfLocalModel[i] = child.getstandardDeviance(i);
    		approximatedObjectivesOfLocalModel[i] = child.getObjective(i);
  		}
  		
  		
  		for(int i=0; i<child.numberOfObjectives(); i++)
  		{	
  			if(child.getstandardDeviance(i) == Double.NaN)
  			{
  				i++;
  				i--;
  			}
  			if(globalModelSolution.getstandardDeviance(i) == Double.NaN)
  			{
  				i++;
  				i--;
  			}
  			if(child.getstandardDeviance(i) > confidenceIntervalsOfGlobalModel[i])
  			{
  				child.setObjective(i, globalModelSolution.getObjective(i));
  				if(globalModelSolution.numberOfFeatures()>0)
  					child.setFeature(i, globalModelSolution.getFeature(i));
  				child.setStandardDeviance(i, globalModelSolution.getstandardDeviance(i));
  				whichModelIsBetter[i] = " |G|";
  			}
  			else
  				whichModelIsBetter[i] = " |L|";
  		}
  		  	
  	}
  	
  	if(problem_.getName() == "StoreSteel" || problem_.getName() == "ExternEvaluator")
  	{
  		
	  	if(problem_.getNumberOfFeatures() > 0)
	  	{	
	  		//calculate objectives from features
	  		for(int i=0; i<problem_.getNumberOfFeatures(); i++)
	  			child.setObjective(i, getObjective(i, child.getFeature(i)));  		
	  	}  	
  	
	  	distanceToFisibility =  calculateDistanceToFeasibility(child);
	  	if(distanceToFisibility != 0)
	  	{
	  		child.setNumberOfViolatedConstraint(1);
	  		child.setOverallConstraintViolation(distanceToFisibility);
	  	}
	  	
	  	else
	  	{
	  		child.setNumberOfViolatedConstraint(0);
	  		child.setOverallConstraintViolation(0);
	  	}
  	}
  	
  	  int counter = 0;
  		allObjectivesAndConfidences[counter] = " G|Obj. ";
  		counter++;
  		for(int j=0; j< child.numberOfObjectives(); j++)  	
  		{
  			allObjectivesAndConfidences[counter] = String.valueOf(approximatedObjectivesOfGlobalModel[j]) + " ";
  			counter++;
  		}

  		allObjectivesAndConfidences[counter] = " G|Conf. ";
  		counter++;
  		for(int j=0; j< child.numberOfObjectives(); j++)  	
  		{
  			allObjectivesAndConfidences[counter] = String.valueOf(confidenceIntervalsOfGlobalModel[j] + " ");
  			counter++;
  		}
  	
  		allObjectivesAndConfidences[counter] = "  L|Obj. ";
  		counter++;
  		for(int j=0; j< child.numberOfObjectives(); j++)  	
  		{
  			allObjectivesAndConfidences[counter] = String.valueOf(approximatedObjectivesOfLocalModel[j] + " ");
  			counter++;
  		}

  		allObjectivesAndConfidences[counter] = " L|Conf. ";
  		counter++;
  		for(int j=0; j< child.numberOfObjectives(); j++)  	
  		{
  			allObjectivesAndConfidences[counter] = String.valueOf(confidenceIntervalsOfLocalModel[j] + " ");
  			counter++;
  		}
  		
  		return allObjectivesAndConfidences;
  }

  
	private double calculateDistanceToFeasibility(Solution solution) {
double distance = 0;
		
		if(solution.getFeature(0) < 10)
			distance = distance + 10 - solution.getFeature(0);
		else if(solution.getFeature(0) > 11)
			distance = distance + solution.getFeature(0) - 11;
		
		if(solution.getFeature(1) < 11)
			distance = distance + 11 - solution.getFeature(1);
		else if(solution.getFeature(1) > 15)
			distance = distance + (solution.getFeature(1) - 15)/4;
		
		if(solution.getFeature(2) < 1115)
			distance = distance + 1115 - solution.getFeature(2);
		else if(solution.getFeature(2) > 1130)
			distance = distance + (solution.getFeature(2) - 1130)/15;
		
	  return (-1 * distance);
  }

	private double getObjective(int i, double feature) {

		double objective;
	  
	  if(i== 0)
	  	objective = Math.abs(feature-10);
	  else if(i== 1)
	  	objective = Math.abs(feature-13);
	  else
	  	objective = Math.abs(feature-1122.5);
	  
	  return objective;
  }

	private int compareAndEvaluateSolutionsUnderDominance(Solution parent, Solution child, SolutionSet population) throws JMException, ClassNotFoundException {

  	Comparator dominance;
  	dominance = new DominanceComparatorUnderUncertanty();
  	int result;
  	
  	result = dominance.compare(parent, child); 
    if (result == -1) { // Solution i dominates child
      return -1;                    
    } // if
    else if (result == 1) { // child dominates
      return 1;                                            
    } // else if
    else if (result == 0){ // the two solutions are non-dominated    	
      return 0;
    }
    else if (result == 2 || result == 10)		//child could be better
    {
     	if(exactEvaluations > limitNumberForExactEvaluations)
    		finishOptimization(population);
     	
    	problem_.evaluate(child) ;
    	child.setStandardDevianceToZero();
      problem_.evaluateConstraints(child);      
      child.setExactllyEvaluated(true);
      setOfExactlyEvaluatedSolutions.add(child);
      exactEvaluations ++;
      
      return compareAndEvaluateSolutionsUnderDominance(parent, child, population);
    }
    else // (result == -2)   parent could be better
    {
    	problem_.evaluate(parent) ;
    	parent.setStandardDevianceToZero();
      problem_.evaluateConstraints(parent);      
      parent.setExactllyEvaluated(true);
      setOfExactlyEvaluatedSolutions.add(parent);
      exactEvaluations ++;
      
      return compareAndEvaluateSolutionsUnderDominance(parent, child, population);
    }
	  
  }

	public SolutionSet environmentalSelectionSPEA2(SolutionSet union, int archiveSize) {
    SolutionSet res;
    /*
     * Racunanje razdalj med osebki. Razdalje so shranjene v spremenljivki distance objekta Spa2Fitness
     */
    Spea2Fitness spea = new Spea2Fitness(union);
      
    /*
     * Calculating:   strength(moc)
     *                R(i) - rawFitness(groba uspesnost)
     *                D(i) - distance (gostota)
     *                F(i) = R(i) + D(i) - (prava uspesnost)
     */
    spea.fitnessAssign();
      
    /*
     * Ce je nedominiranih osebkov vec kot je velikost novega arhiva, izloci osebke, ki so blizu drugim osebkom (4.3)
     * Sicer pa, ce nedominirani osebki ne zapolnijo arhiva, arhiv dopolni z najboljsimi dominiranimi osebki iz P(t-1) in A(t-1) (4.4)
     * ce je zaustavitveni kritrij izpolnjen, koncaj (4.5)
     */    
    res = spea.environmentalSelection(archiveSize);  
      
    return res;
  }

  public SolutionSet environmentalSelectionNSGAIIWithUncertanty(SolutionSet population) throws JMException, ClassNotFoundException{
    SolutionSet finalPop;
    Distance distance;
    RankingUnderUncertantyZaStPrimerjav ranking;
  	distance = new Distance();     
    int remain = populationSize;
    int index  = 0;
    SolutionSet front = null;
    finalPop = new SolutionSet(remain);
    int tempSizeOfSet = setOfExactlyEvaluatedSolutions.size();
        
    ranking = new RankingUnderUncertantyZaStPrimerjav(population, problem_, modelGlobal_, problem_ );
    
    //steviloVsehPrimerjav = steviloVsehPrimerjav + ranking.steviloPrimerjav;
    //steviloNOBBCorrect = steviloNOBBCorrect + ranking.NOBBSteviloPravilnihPrimerjav;
    //steviloNOBBIncorrect = steviloNOBBIncorrect + ranking.NOBBSteviloNepravilnihPrimerjav;
    
    exactEvaluations = exactEvaluations + setOfExactlyEvaluatedSolutions.size() - tempSizeOfSet;
    
    if(exactEvaluations > limitNumberForExactEvaluations)
  		finishOptimization(population);
    
    front = ranking.getSubfront(index);

    while ((remain > 0) && (remain >= front.size())) { 
      // Assign crowding distance to individuals
      distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
      
      for (int k = 0; k < front.size(); k++)
        finalPop.add(front.get(k));

      // Decrement remain
      remain = remain - front.size();

      // Obtain the next front
      index++;
      if (remain > 0)
        front = ranking.getSubfront(index);
    } // while

    // remain is less than front(index).size, insert only the best one
    if (remain > 0) { // front contains individuals to insert                        
      while (front.size() > remain) {
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        front.remove(front.indexWorst(new CrowdingComparator()));
      }
      for (int k = 0; k < front.size(); k++)
        finalPop.add(front.get(k));
      remain = 0; 
    } // if    
    return finalPop;      
  }
  
  private void finishOptimization(SolutionSet population) throws JMException {
	  
  	//System.out.println("Zgoden konec! Exact: "+exactEvaluations);
  	population.printToFile(genFileName, "Exceeden number of exact evaluations! Exactly evaluated solutions: "+exactEvaluations);
  	
  	//Reduce number of solutions because of exeded number of exact eval. we can have more than 100 solutions in the final front
  	if(population.size() > populationSize)
  	{
  		population = environmentalSelectionNSGAII(new Ranking(population), new Distance(), new SolutionSet(), populationSize, 0);	//copied from DEMO
  	}
  	
  	additionalEvaluations = 0;
    Ranking ranking = new Ranking(population);
    for(int i=0; i< ranking.getSubfront(0).size(); i++)
    {
    	if(ranking.getSubfront(0).get(i).isExactllyEvaluated() == false)
    	{
    		problem_.evaluate(ranking.getSubfront(0).get(i));
    		ranking.getSubfront(0).get(i).setStandardDevianceToZero();
    		problem_.evaluateConstraints(ranking.getSubfront(0).get(i));    		
    		ranking.getSubfront(0).get(i).setExactllyEvaluated(true);
    		additionalEvaluations++;    		
    	}
    }    	
   
    
    Hypervolume hyp = new Hypervolume();
    population.printToFile(genFileName, "Additional evaluations: " + additionalEvaluations);
    population.printToFile(genFileName, "Final hypervolume: " + hyp.hypervolume(ranking.getSubfront(0).writeObjectivesToMatrix(), Params.getFrontForHypervolume(), problem_.getNumberOfObjectives()));
    population.printToFile(genFileName, "Nondominated solutions: " + ranking.getSubfront(0).size());
    
    population.printToFile(frontFileName, "---------------------------------------------------------------------------------------------------");
    
    if (frontFormat == 1) // decision and objective vector
      ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
    else if (frontFormat == 0) // only objective vector
      ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);  
    
    // Return the first non-dominated front
    System.out.println("Dodatni izracuni " + additionalEvaluations);
	  
    System.exit(0);
  }

  public SolutionSet environmentalSelectionNSGAII(Ranking ranking, Distance distance, SolutionSet front, int remain, int index){
    SolutionSet finalPop=new SolutionSet(remain);

    // Obtain the next front
    front=ranking.getSubfront(index);

    while ((remain > 0) && (remain >= front.size())){ 
	    //Assign crowding distance to individuals
	    distance.crowdingDistanceAssignment(front,problem_.getNumberOfObjectives());
	
	    for (int k = 0; k < front.size(); k++ ) {
        finalPop.add(front.get(k));
	    } // for
	
	    //Decrement remain
	    remain = remain - front.size();
	
	    //Obtain the next front
	    index++;
	    if (remain > 0) {
        front = ranking.getSubfront(index);
	    } // if
    }//end while

     // remain is less than front(index).size, insert only the best one
    if (remain > 0) {  // front contains individuals to insert                        
	    while (front.size() > remain) {
	    	distance.crowdingDistanceAssignment(front,problem_.getNumberOfObjectives());
	    	front.remove(front.indexWorst(new CrowdingComparator()));
	    }
	    for (int k = 0; k < front.size(); k++) {
	    	finalPop.add(front.get(k));
	    }
	    remain = 0; 
	  } // if    
	  	  
	  return finalPop;    
  }
  
	private boolean areSolutionsTheSame(Solution parent, Solution child) throws JMException {
    for (int i = 0; i < Params.getNumberOfVariables(); i++) {
      if (parent.getDecisionVariables()[i].getValue() != child.getDecisionVariables()[i].getValue())
        return false;
      
    }
    return true;
  }

} // MetaDEMO

class UpdateGlobalModel extends RecursiveTask
{
	Model model;
	String modelName;
	SolutionSet setOfExactlyEvaluatedSolutions;
	boolean b;

	public UpdateGlobalModel(Model model,  SolutionSet setOfExactlyEvaluatedSolutions, boolean b) throws ClassNotFoundException, JMException 
	{
		this.model = model;
		this.setOfExactlyEvaluatedSolutions = setOfExactlyEvaluatedSolutions;
		this.b = b;
				
  }

	@Override
  public String compute() {

			System.out.println("Zacetek global");
			if(model.getName() == "SPGP")
	      try 
				{
	        ((SPGPm)(model)).update(setOfExactlyEvaluatedSolutions, b);
        } 
				catch (ClassNotFoundException | JMException e) 
				{
	        e.printStackTrace();
        }
			
      else
      {
	      try 
	      {
	        model.update(setOfExactlyEvaluatedSolutions);
        } 
	      catch (ClassNotFoundException | JMException e) 
	      {
	        e.printStackTrace();
        }
      }
			System.out.println("Konec global");
			return "Global konc";
		}
  
} 

class UpdateLocalModel extends RecursiveTask
{
	Model model;
	String modelName;
	SolutionSet setOfExactlyEvaluatedSolutions;
	boolean b;

	public UpdateLocalModel(Model model,  SolutionSet setOfExactlyEvaluatedSolutions) throws ClassNotFoundException, JMException 
	{
		this.model = model;
		this.setOfExactlyEvaluatedSolutions = setOfExactlyEvaluatedSolutions;			
  }

	@Override
  public String compute() {
		System.out.println("Zacetek local");
		try {
      model.update(setOfExactlyEvaluatedSolutions);
    } catch (ClassNotFoundException | JMException e) {
      e.printStackTrace();
    }	
		System.out.println("Konec local");
		return "Local konc";
  }
} 

class UpdateModelParalel implements Runnable
{
	Model model;
	String modelName;
	SolutionSet setOfExactlyEvaluatedSolutions;
	boolean b;
	
	public UpdateModelParalel(Model model, String modelName,  SolutionSet setOfExactlyEvaluatedSolutions, boolean b) throws ClassNotFoundException, JMException 
	{
		this.model = model;
		this.modelName = modelName;
		this.setOfExactlyEvaluatedSolutions = setOfExactlyEvaluatedSolutions;
		this.b = b;
				
  }

	@Override
  public void run() {
		if(modelName == "localModel")
		{
System.out.println("Zacetek local");
			try {
	      model.update(setOfExactlyEvaluatedSolutions);
      } catch (ClassNotFoundException | JMException e) {
	      e.printStackTrace();
      }	
			System.out.println("Konec local");
		}
		
		else  //global model
		{
			System.out.println("Zacetek global");
			if(model.getName() == "SPGP")
	      try 
				{
	        ((SPGPm)(model)).update(setOfExactlyEvaluatedSolutions, b);
        } 
				catch (ClassNotFoundException | JMException e) 
				{
	        e.printStackTrace();
        }
			
      else
      {
	      try 
	      {
	        model.update(setOfExactlyEvaluatedSolutions);
        } 
	      catch (ClassNotFoundException | JMException e) 
	      {
	        e.printStackTrace();
        }
      }
			System.out.println("Konec global");
		}
  }
  
}